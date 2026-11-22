"""Figure 2 from newly fitted study and background-site regression tables."""
import argparse,csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--study',type=Path,required=True);ap.add_argument('--controls',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    def read(p):
        with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
    rows=[r for r in read(a.study)+read(a.controls) if r['ventilation']=='vc10']
    order=['algiers','arzew','oran','setif','constantine','hassi_messaoud','control_tell','control_highland'];lookup={r['site']:r for r in rows}
    fig,ax=plt.subplots(figsize=(8.4,5.5),dpi=200)
    for i,s in enumerate(order):
        r=lookup[s];b=float(r['slope_lnNO2_per_lnVC']);lo=float(r['ci95_low']);hi=float(r['ci95_high']);color='#b5421f' if i<6 else '#364b59'
        ax.errorbar(b,i,xerr=[[b-lo],[hi-b]],fmt='o',color=color,capsize=3)
    ax.axvline(0,color='0.4',ls='--');ax.set_yticks(range(len(order)),[s.replace('_',' ').title() for s in order]);ax.invert_yaxis();ax.set_xlabel('Elasticity of NO₂ to ventilation (VC10), 95% CI');fig.tight_layout();a.out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.out);plt.close(fig)
if __name__=='__main__':main()
