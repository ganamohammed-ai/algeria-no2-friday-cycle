"""Reconstruct the manuscript's descriptive reference-window ranges on the log scale.

This is an export step added for the publication package, not a new statistical test.
The ranges are not coverage-validated confidence intervals.
"""
import argparse,csv,math
from pathlib import Path
import numpy as np

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    def read(name):
        with (args.results/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
    windows=read('detect_calibration_null_windows.csv');estimates=read('detect_calibration_2020.csv')
    rows=[]
    for r in estimates:
        if r['model']!='adjusted' or r['se_method']!='nw7rows':continue
        p=np.array([float(w['est_pct']) for w in windows if w['site']==r['site'] and w['model']=='adjusted'])
        beta=math.log1p(float(r['est_pct'])/100);sd=np.log1p(p/100).std(ddof=1)
        rows.append(dict(site=r['site'],estimate_pct=float(r['est_pct']),nw_lo_pct=float(r['lo_pct']),nw_hi_pct=float(r['hi_pct']),
                         descriptive_lo_pct=float(100*np.expm1(beta-1.96*sd)),descriptive_hi_pct=float(100*np.expm1(beta+1.96*sd)),
                         sd_log_coefficients_times100=float(100*sd),mean_reference_shift_pct=float(p.mean()),
                         reference_windows_lower_pct=float(100*np.mean(p<float(r['est_pct']))),n_reference_windows=len(p)))
    with (args.out/'table8_reference_window_ranges.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('Exported Table 8:',len(rows),'sites')
if __name__=='__main__':main()
