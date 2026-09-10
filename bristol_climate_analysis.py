#!/usr/bin/env python3
"""Export overlapping raw monthly values from the Bristol 5.1 JSON.

This helper does not reproduce the interactive baseline/lag inference.
Use the dashboard CSV for that selection. No confidence claim is made here.

No network access is required. Example:
  python bristol_climate_analysis.py ../data/bristol-climate.json \
      --a swep_rain_mm:value --b frome_monthly_flow_m3s:value --out comparison.csv
"""
import argparse,csv,json,math,statistics
def mean(x): return statistics.fmean(x)
def corr(a,b):
 if len(a)<3 or len(a)!=len(b) or len(set(a))<2 or len(set(b))<2:return None
 ma,mb=mean(a),mean(b); n=sum((x-ma)*(y-mb) for x,y in zip(a,b))
 return n/math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
def lag1(x): return corr(x[:-1],x[1:]) if len(x)>3 else 0
def extract(doc,spec):
 dataset,key=spec.split(':',1)
 return {r['date'][:7]:float(r[key]) for r in doc['series'][dataset] if r.get(key) is not None}
def main():
 p=argparse.ArgumentParser();p.add_argument('json');p.add_argument('--a',required=True);p.add_argument('--b',required=True);p.add_argument('--out',default='comparison.csv');a=p.parse_args()
 d=json.load(open(a.json,encoding='utf-8')); x,y=extract(d,a.a),extract(d,a.b); dates=sorted(set(x)&set(y)); av=[x[k] for k in dates];bv=[y[k] for k in dates]
 r=corr(av,bv)
 with open(a.out,'w',newline='',encoding='utf-8') as f:
  w=csv.writer(f);w.writerow(['date',a.a,a.b]);w.writerows((k,x[k],y[k]) for k in dates)
 print(json.dumps({'version':'5.1.0','mode':'raw monthly values; seasonality retained; no statistical inference','n':len(dates),'pearson_r':r,'retrieved_at':d['updated_at']},indent=2,allow_nan=False))
if __name__=='__main__': main()
