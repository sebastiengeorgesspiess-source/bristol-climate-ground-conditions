"""Reproduce Bristol 5.1 annual indices and inference offline.
Usage: python bristol_extremes_reproducible.py path/to/bristol-climate-daily.json
Install numpy and scipy. The snapshot date fixes the complete-year cutoff.
Regional temperature and precipitation are never treated as Bristol station data.
"""
import json, sys, calendar, statistics
from pathlib import Path
from datetime import date as RealDate
import numpy as np
from scipy.stats import genextreme, kendalltau, theilslopes, linregress

def annual_extremes(rain,mean_t,max_t,min_t):
 def group(rows):
  out={}
  for r in rows:out.setdefault(r['date'][:4],[]).append(r)
  return out
 rg,tg,xg,ng=map(group,(rain,mean_t,max_t,min_t)); out=[]
 years=sorted(set(rg)|set(tg)|set(xg)|set(ng))
 wet_base=[r['value'] for r in rain if '1961-01-01'<=r['date']<='1990-12-31' and r['value']>=1]
 p95=float(np.percentile(wet_base,95)) if wet_base else None
 for y in years:
  rr=sorted(rg.get(y,[]),key=lambda z:z['date']); rv=[r['value'] for r in rr]; xv=[r['value'] for r in xg.get(y,[])]; nv=[r['value'] for r in ng.get(y,[])]; tv=[r['value'] for r in tg.get(y,[])]
  expected=366 if calendar.isleap(int(y)) else 365
  if int(y)>=date.today().year:continue
  rv=rv if len(rv)==expected else [];xv=xv if len(xv)==expected else [];nv=nv if len(nv)==expected else [];tv=tv if len(tv)==expected else []
  if rv:
   rolls=[sum(rv[i-4:i+1]) for i in range(4,len(rv))]; runs=[]; run=0
   for v in rv:
    run=run+1 if v<1 else 0;runs.append(run)
  else:rolls=[];runs=[]
  row={'year':int(y),'rain_days':len(rv),'temp_days':len(tv)}
  if rv:row.update({'rain_total_mm':round(sum(rv),1),'rx1day_mm':round(max(rv),1),'rx5day_mm':round(max(rolls),1),'wet_days':sum(v>=1 for v in rv),'r10_days':sum(v>=10 for v in rv),'r20_days':sum(v>=20 for v in rv),'cdd_days':max(runs),'r95p_mm':round(sum(v for v in rv if p95 is not None and v>p95),1)})
  if tv:row['tmean_c']=round(statistics.mean(tv),2)
  if xv:row.update({'txx_c':round(max(xv),1),'summer_days':sum(v>25 for v in xv),'ice_days':sum(v<0 for v in xv)})
  if nv:row.update({'tnn_c':round(min(nv),1),'frost_days':sum(v<0 for v in nv),'tropical_nights':sum(v>20 for v in nv)})
  if xv and nv:
   n=min(len(xv),len(nv));row['dtr_c']=round(statistics.mean(xv[i]-nv[i] for i in range(n)),2)
  out.append(row)
 return out,p95

def pettitt(values):
 x=np.asarray(values,float);n=len(x);u=[]
 for t in range(n-1):u.append(sum(np.sign(x[t+1:]-x[:t+1,None]).ravel()))
 k=int(np.argmax(np.abs(u)));K=abs(u[k]);p=min(1.0,2*np.exp((-6*K*K)/(n**3+n**2)))
 return k,float(K),float(p)

def trend_test(annual,key):
 z=[(r['year'],r[key]) for r in annual if r.get(key) is not None];years=np.array([x[0] for x in z],float);v=np.array([x[1] for x in z],float)
 tau,p=kendalltau(years,v);sen=theilslopes(v,years,.95);ols=linregress(years,v);k,K,pp=pettitt(v)
 return {'variable':key,'n':len(v),'mann_kendall_tau':round(float(tau),4),'mann_kendall_p':round(float(p),6),'theil_sen_per_decade':round(float(sen.slope*10),4),'theil_sen_95ci_per_decade':[round(float(sen.low_slope*10),4),round(float(sen.high_slope*10),4)],'ols_per_decade':round(float(ols.slope*10),4),'ols_p':round(float(ols.pvalue),6),'pettitt_year':int(years[k]),'pettitt_p_approx':round(pp,6),'warning':'Pettitt identifies one candidate shift; metadata and multiple-testing context are required before interpretation.'}

def gev_bootstrap(annual,key,nboot=1000,seed=20260831):
 z=np.asarray([r[key] for r in annual if r.get(key) is not None],float);c,loc,scale=genextreme.fit(z);periods=(20,50,100);base=np.array([genextreme.ppf(1-1/t,c,loc=loc,scale=scale) for t in periods]);rng=np.random.default_rng(seed);pars=[];levels=[]
 for _ in range(nboot):
  try:
   cc,ll,ss=genextreme.fit(rng.choice(z,len(z),replace=True));pars.append((cc,ll,ss));levels.append([genextreme.ppf(1-1/t,cc,loc=ll,scale=ss) for t in periods])
  except:pass
 def ci(a,i):return [round(float(np.percentile(a[:,i],2.5)),4),round(float(np.percentile(a[:,i],97.5)),4)]
 pa=np.asarray(pars);la=np.asarray(levels)
 return {'variable':key,'method':'GEV MLE with nonparametric year bootstrap','seed':seed,'requested_bootstrap':nboot,'successful_bootstrap':len(pa),'parameters':{'shape':{'estimate':round(float(c),5),'ci95':ci(pa,0)},'location':{'estimate':round(float(loc),4),'ci95':ci(pa,1)},'scale':{'estimate':round(float(scale),4),'ci95':ci(pa,2)}},'return_levels':{str(t):{'estimate':round(float(base[i]),2),'ci95':ci(la,i)} for i,t in enumerate(periods)}}

if __name__ == '__main__':
    source=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'data/bristol-climate-daily.json'
    snapshot=json.loads(source.read_text())
    class date(RealDate):
        @classmethod
        def today(cls): return RealDate.fromisoformat(snapshot['updated_at'][:10])
    d=snapshot['series']
    annual,p95=annual_extremes(d['swep_daily'],d['cet_daily_mean'],d['cet_daily_max'],d['cet_daily_min'])
    result={'version':'5.1.0','cutoff':str(date.today()),'annual':annual,'wet_day_p95':p95,
            'trends':[trend_test(annual,k) for k in ('rx1day_mm','rx5day_mm','cdd_days','wet_days','summer_days','frost_days','txx_c','tnn_c')],
            'bootstrap':[gev_bootstrap(annual,k) for k in ('rx1day_mm','rx5day_mm','txx_c')]}
    print(json.dumps(result,indent=2,allow_nan=False))
