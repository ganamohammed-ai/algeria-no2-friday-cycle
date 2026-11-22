"""Export the four within-site Friday subperiod estimates described in the article."""
import argparse,csv
from pathlib import Path
import friday_cycle as F

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--table',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    rows=F.A.read_csv(a.table);results=[]
    for start,end in [('2019-08-07','2020-02-29'),('2020-03-01','2021-06-30'),('2021-07-01','2023-12-31'),('2024-01-01','2025-12-31')]:
        for site in F.STUDY:
            sub=[r for r in rows if r['site']==site and start<=r['date']<=end]
            n,o=F.fit(sub,2,'fs')
            results.append(dict(site=site,start=start,end=end,n=n,friday_pct=o['Fri'][0],friday_lo=o['Fri'][1],friday_hi=o['Fri'][2],saturday_pct=o['Sat'][0],saturday_lo=o['Sat'][1],saturday_hi=o['Sat'][2]))
    with (a.out/'friday_subperiods.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    print('Exported',len(results),'site-period estimates')
if __name__=='__main__':main()
