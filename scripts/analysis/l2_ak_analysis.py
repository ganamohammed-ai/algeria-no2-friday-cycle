#!/usr/bin/env python3
"""Elasticity of NO2 to ventilation (VC10) with L2 columns: as retrieved, and with replaced profile shapes.
Usage: python l2_ak_analysis.py --table no2_ventilation_daily_v2.csv --l2 l2_ak_daily.csv --out outputs
"""
import argparse, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_no2_ventilation as A

SITES = ["algiers", "oran", "arzew", "constantine", "setif", "hassi_messaoud"]
MIN_PIX = 5

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True); ap.add_argument("--l2", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    main_rows = A.read_csv(a.table)
    l2 = {(r["site"], r["date"]): r for r in A.read_csv(a.l2)}
    out, desc = [], []
    for s in SITES:
        rows = []
        for r in main_rows:
            if r["site"] != s or r["usable"] != 1:
                continue
            q = l2.get((s, r["date"]))
            if q is None or q["n_pix"] < MIN_PIX:
                continue
            r = dict(r)
            r.update(l2_vcd=q["vcd"], l2_bl=q["vcd_bl"], l2_un=q["vcd_un"], f_bl=q["f_bl"], f_un=q["f_un"],
                     amf=q["amf_trop"], n_pix=q["n_pix"])
            rows.append(r)
        row = {"site": s}
        for key, lab in (("no2_inner", "L3"), ("l2_vcd", "L2"), ("l2_bl", "L2_bl"), ("l2_un", "L2_un"),
                         ("f_bl", "ln_f_bl"), ("f_un", "ln_f_un"), ("amf", "ln_amf")):
            f = A.fit_site(rows, "vc10", y_key=key)
            row["n"] = f["n"]
            row[lab + "_slope"] = f["slope_lnNO2_per_lnVC"]
            row[lab + "_lo"] = f["ci95_low"]; row[lab + "_hi"] = f["ci95_high"]
        row["median_f_bl"] = float(np.median([r["f_bl"] for r in rows]))
        row["median_f_un"] = float(np.median([r["f_un"] for r in rows]))
        row["corr_L3_L2"] = float(np.corrcoef(np.log([r["no2_inner"] for r in rows if r["no2_inner"] > 0 and r["l2_vcd"] > 0]),
                                               np.log([r["l2_vcd"] for r in rows if r["no2_inner"] > 0 and r["l2_vcd"] > 0]))[0, 1])
        out.append(row)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()})
    path = os.path.join(a.out, "l2_ak_elasticity.csv")
    A.write_csv(path, out)
    print("wrote", path)

if __name__ == "__main__":
    main()
