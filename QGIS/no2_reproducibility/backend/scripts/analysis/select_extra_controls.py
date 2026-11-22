#!/usr/bin/env python3
"""Analysis D site selection (plan v2, section 6). Static layers only; NO NO2 is read.
Eligibility = v1 base rule + the two v1 amendments (lat >= 34.5, dist to nearest study site <= 150 km):
  built10<=0.5%, built50<=2%, pop10<=5000, nl10<=1.0, water10<=5%, dist>=100 km.
Greedy by lowest night light (ties within 0.05: lowest built10, then elev sd); a candidate is kept only if
>=100 km from every already selected site, INCLUDING the two existing controls. Max 10 new sites."""
import csv, math, sys
def f(x):
    try: return float(x)
    except: return float("nan")
def hav(a,b):
    R=6371.0; p1,p2=math.radians(a[0]),math.radians(b[0]); dl=math.radians(b[1]-a[1]); dp=p2-p1
    h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(math.sqrt(h))
rows=[]
for r in csv.DictReader(open(sys.argv[1],encoding="utf-8")):
    d=dict(lat=f(r["lat"]),lon=f(r["lon"]),b10=f(r["built_10km"]),b50=f(r["built_50km"]),water=f(r.get("water_10km") or 0),
           pop=f(r["pop_10km"]),nl=f(r["nightlight_10km"]),elev=f(r["elev_mean_10km"]),sd=f(r["elev_sd_10km"]),dist=f(r["dist_nearest_site_km"]))
    if d["water"]!=d["water"]: d["water"]=0.0
    if any(v!=v for v in d.values()): continue
    rows.append(d)
el=[r for r in rows if r["lat"]>=34.5 and r["dist"]<=150 and r["dist"]>=100 and r["b10"]<=0.005 and r["b50"]<=0.02
    and r["pop"]<=5000 and r["nl"]<=1.0 and r["water"]<=0.05]
print("candidates",len(rows),"eligible",len(el))
el.sort(key=lambda r:(round(r["nl"]/0.05),r["b10"],r["sd"],r["nl"]))
sel=[(35.0,5.25),(34.5,-1.0)]; new=[]
for r in el:
    if all(hav((r["lat"],r["lon"]),s)>=100 for s in sel):
        new.append(r); sel.append((r["lat"],r["lon"]))
    if len(new)==10: break
for i,r in enumerate(new,1): r["name"]=f"ctl_x{i:02d}"
w=csv.DictWriter(open("control_sites_extra_selected.csv","w",newline=""),fieldnames=["name"]+[k for k in new[0] if k!="name"] if new else ["name"]);w.writeheader();w.writerows(new)
for r in new: print(r["name"],r["lat"],r["lon"],"elev",round(r["elev"]),"nl",round(r["nl"],2),"dist",round(r["dist"]))
print("n new",len(new))
