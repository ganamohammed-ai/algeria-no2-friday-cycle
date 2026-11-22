"""Added-control injection experiments; preserves the original model and options."""
import argparse,sys,os
import detect_injection as DI
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--table',required=True);ap.add_argument('--out',required=True);ap.add_argument('--gb',action='store_true');a=ap.parse_args()
    os.makedirs(a.out,exist_ok=True)
    DI.L.STUDY=[];DI.L.CONTROLS=['ctl_x01','ctl_x02','ctl_x03','ctl_x04']
    sys.argv=['detect_injection','--table',a.table,'--controls',a.table,'--out',a.out]+(['--gb'] if a.gb else [])
    DI.main()
if __name__=='__main__':main()
