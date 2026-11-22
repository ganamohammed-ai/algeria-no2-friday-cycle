"""Portable replay of publication analyses from the archived daily tables.

Default mode covers central Friday results and the reference-range export.
Extended mode also covers ventilation, seasonal contrasts, sensitivity and figures.
Full mode additionally reruns expensive containment/ML, injection and calibration.
All writes go to --out; inputs and reference results remain unchanged.
"""
from pathlib import Path
import argparse,subprocess,sys,os,csv,json,math,shutil,time

ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-root',type=Path,default=ROOT.parent/'Zenodo/processed')
    ap.add_argument('--out',type=Path,default=ROOT/'reproduced_results')
    ap.add_argument('--mode',choices=['quick','extended','full'],default='quick')
    args=ap.parse_args();data=args.data_root.resolve();out=args.out.resolve();ref=data/'outputs'
    if out==data or out==ref or data in out.parents:raise SystemExit('--out must be outside the input-data tree.')
    out.mkdir(parents=True,exist_ok=True);logs=out/'logs';logs.mkdir(exist_ok=True)
    study=ref/'no2_ventilation_daily_v2.csv';control=ref/'control/control_no2_ventilation_daily.csv'
    comparisons=[];commands=[]
    env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONIOENCODING='utf-8',MPLCONFIGDIR=str(out/'matplotlib_cache'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    def run(label,script,*arguments):
        cmd=[sys.executable,str(ROOT/'scripts/analysis'/script),*[str(x) for x in arguments]];start=time.time()
        p=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8')
        (logs/(label+'.txt')).write_text(p.stdout+'\n'+p.stderr,encoding='utf-8')
        commands.append({'label':label,'command':cmd,'exit_code':p.returncode,'seconds':time.time()-start})
        print(label,'OK' if p.returncode==0 else 'FAILED',flush=True)
        if p.returncode:raise RuntimeError('See '+str(logs/(label+'.txt')))
    def compare(label,actual,reference,tolerance=1e-7):
        def read(p):
            with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
        aa,bb=read(actual),read(reference);errors=[]
        if len(aa)!=len(bb):errors.append('row count')
        for i,(a,b) in enumerate(zip(aa,bb)):
            if set(a)!=set(b):errors.append('columns');break
            for k in b:
                try:
                    x,y=float(a[k]),float(b[k])
                    if math.isnan(x) and math.isnan(y):continue
                    if not math.isclose(x,y,rel_tol=tolerance,abs_tol=tolerance):errors.append([i,k,a[k],b[k]])
                except (ValueError,TypeError):
                    if a[k]!=b[k]:errors.append([i,k,a[k],b[k]])
        comparisons.append({'check':label,'rows':len(aa),'matches':not errors,'tolerance':tolerance,'differences':errors[:20]})
        if errors:print('REFERENCE DIFFERENCE:',label,errors[:2],flush=True)
    run('friday_original','friday_cycle.py','--table',study,'--controls',control,'--out',out)
    compare('Friday original',out/'friday_cycle.csv',ref/'friday_cycle.csv')
    for label,script,table,reference in [('holdout','detect_friday_holdout.py',ref/'study_2026/no2_ventilation_daily.csv',ref/'study_2026/detect_friday_holdout_2026.csv'),('cities','detect_friday_cities.py',ref/'cities/no2_ventilation_daily.csv',ref/'cities/detect_friday_cities.csv')]:
        a=[table]+([ref/'friday_cycle.csv'] if label=='holdout' else [])
        run(label,script,*a,'--out',out);compare(label,out/reference.name,reference)
    for label,inputs,reference in [('2019_2025',[study],ref/'detect_friday_se_check_2019_2025.csv'),('2026',[ref/'study_2026/no2_ventilation_daily.csv'],ref/'detect_friday_se_check_2026.csv'),('pooled',[study,ref/'study_2026/no2_ventilation_daily.csv',ref/'cities/no2_ventilation_daily.csv'],ref/'detect_friday_se_check.csv')]:
        dest=out/('se_'+label);run('se_'+label,'detect_friday_se_check.py',dest,*inputs);compare('SE '+label,dest/'detect_friday_se_check.csv',reference)
    run('table8','export_table8.py','--results',ref,'--out',out)
    run('friday_subperiods','export_friday_subperiods.py','--table',study,'--out',out)
    if args.mode in ['extended','full']:
        run('lockdown_ols','lockdown_normalisation.py','--table',study,'--controls',control,'--out',out,'--skip-ml')
        for name in ['lockdown_effect.csv','lockdown_placebo.csv','trend_version_steps.csv']:
            compare(name,out/name,ref/name,1e-5)
        run('seasonal','detect_seasonal.py','--table',study,'--controls',control,'--out',out)
        compare('Seasonal windows',out/'detect_seasonal.csv',ref/'detect_seasonal.csv',1e-5)
        for script,name,extra in [('decompose_blh_wind.py','sensitivity_blh_wind_split.csv',[]),('sensitivity_version_steps.py','sensitivity_version_steps.csv',['--controls',control]),('l2_ak_analysis.py','l2_ak_elasticity.csv',['--l2',ref/'l2_ak_daily.csv'])]:
            run(script,script,'--table',study,*extra,'--out',out);compare(name,out/name,ref/name,1e-5)
        # Model-stage CLI reads this basename. This is a copy in the output tree.
        for tag,table in [('study',study),('control',control)]:
            dest=out/('model_'+tag);dest.mkdir(exist_ok=True);shutil.copy2(table,dest/'no2_ventilation_daily.csv')
            sites='algiers,oran,arzew,constantine,setif,hassi_messaoud' if tag=='study' else 'control_tell,control_highland'
            run('ventilation_'+tag,'analyze_no2_ventilation.py','--stage','model','--sites',sites,'--out',dest)
        dustdir=data/'data/04_dust_aai'
        if dustdir.exists():
            dest=out/'model_study'
            run('dust_and_seasons','sensitivity_analyses.py','--out',dest,'--aai-dir',dustdir)
            for name in ['sensitivity_outliers.csv','sensitivity_dust.csv','sensitivity_season_ring.csv']:
                compare(name,dest/name,ref/name,1e-4)
        run('figure1','fig_no2_vs_ventilation.py','--table',study,'--out',out/'fig1_no2_vs_ventilation.png')
        run('figure2','fig_elasticity_controls.py','--study',out/'model_study/regression_results.csv','--controls',out/'model_control/regression_results.csv','--out',out/'fig_elasticity_with_controls.png')
        run('figure3','fig_friday.py','--effects',out/'friday_cycle.csv','--out',out/'fig_friday_effect.png')
        run('figure4','fig_lockdown.py','--effects',ref/'lockdown_effect.csv','--out',out/'fig_lockdown_effect.png')
        run('figureS1','fig_wind_rotation.py','--input-dir',ref/'wind_rotation','--out',out/'fig_wind_rotation.png')
    if args.mode=='full':
        run('lockdown_full','lockdown_normalisation.py','--table',study,'--controls',control,'--out',out)
        for name in ['lockdown_effect.csv','lockdown_placebo.csv','lockdown_ml.csv','trend_version_steps.csv']:compare(name,out/name,ref/name,1e-4)
        run('injection_full','detect_injection.py','--table',study,'--controls',control,'--out',out,'--gb')
        extra=out/'control_extra';extra.mkdir(exist_ok=True)
        run('extra_injection','run_extra_detect.py','--table',ref/'control_extra/no2_ventilation_daily.csv','--out',extra,'--gb')
        run('calibration','detect_calibrate.py',out,study,control,ref/'control_extra/no2_ventilation_daily.csv')
        run('table8_recalculated','export_table8.py','--results',out,'--out',out)
    summary={'mode':args.mode,'python':sys.version,'commands':commands,'comparisons':comparisons,'all_compared_results_match':all(x['matches'] for x in comparisons),'scope':'This replay starts from processed daily tables. Raw satellite extraction is documented separately; archived binary files were not all reprocessed.'}
    (out/'validation_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    if not summary['all_compared_results_match']:raise SystemExit('Reference differences found; see validation_summary.json')
    print('Replay passed:',len(comparisons),'numeric comparisons.',flush=True)
if __name__=='__main__':main()
