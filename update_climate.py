#!/usr/bin/env python3
from urllib.request import Request,urlopen
from datetime import date,datetime,timezone,timedelta
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,hashlib,calendar,time,statistics,urllib.parse,gzip,io
import numpy as np
from scipy.stats import genextreme,kendalltau,theilslopes,linregress
from scipy.optimize import minimize
from shapely.geometry import shape,mapping
BASE=Path(__file__).resolve().parent; OUT=BASE/'data'; OUT.mkdir(exist_ok=True); HEAD={'User-Agent':'Sebastien-Spiess-Climate-Commission/1.0'}
def text(u):
 for attempt in range(6):
  try:return urlopen(Request(u,headers=HEAD),timeout=120).read().decode()
  except Exception as exc:
   if getattr(exc,'code',None)!=429 or attempt==5:raise
   time.sleep(5*(attempt+1))
def js(u): return json.loads(text(u))
def payload(u):
 b=urlopen(Request(u,headers=HEAD),timeout=180).read()
 return b,json.loads(b)
def matrix(raw,start=0):
 out=[]
 for line in raw.splitlines():
  p=line.split()
  if len(p)<2 or not p[0].isdigit() or len(p[0])!=4: continue
  y=int(p[0])
  for m,v in enumerate(p[1:13],1):
   try:v=float(v)
   except:continue
   if np.isfinite(v) and v>-90 and y>=start and date(y,m,1)<date.today().replace(day=1):out.append({'date':f'{y:04d}-{m:02d}','value':v})
 return out
def era_chunk(years):
 a,b=years; end=min(date.today()-timedelta(days=7),date(b,12,31))
 daily='temperature_2m_mean,relative_humidity_2m_mean,precipitation_sum,surface_pressure_mean,soil_moisture_0_to_7cm_mean,soil_moisture_7_to_28cm_mean,soil_moisture_28_to_100cm_mean,wind_speed_10m_mean,et0_fao_evapotranspiration_sum'
 u=f'https://archive-api.open-meteo.com/v1/archive?latitude=51.4545&longitude=-2.5879&start_date={a}-01-01&end_date={end.isoformat()}&daily={daily}&timezone=UTC&models=era5'
 return js(u)['daily']
def monthly_era(chunks):
 sums={}; means={'temperature_2m_mean','relative_humidity_2m_mean','surface_pressure_mean','soil_moisture_0_to_7cm_mean','soil_moisture_7_to_28cm_mean','soil_moisture_28_to_100cm_mean','wind_speed_10m_mean'}
 for d in chunks:
  for i,day in enumerate(d['time']):
   key=day[:7]; b=sums.setdefault(key,{})
   for var,arr in d.items():
    if var=='time' or i>=len(arr) or arr[i] is None:continue
    q=b.setdefault(var,[]); q.append(float(arr[i]))
 out=[]
 for key,b in sorted(sums.items()):
  row={'date':key}
  expected=calendar.monthrange(int(key[:4]),int(key[5:]))[1]
  for var,v in b.items():row[var]=(round(sum(v)/len(v),3) if var in means else round(sum(v),3)) if len(v)==expected else None
  row['valid_days']={var:len(v) for var,v in b.items()};row['expected_days']=expected
  out.append(row)
 return out
def bgs_features(collection,limit=1000):
 u=f'https://ogcapi.bgs.ac.uk/collections/{collection}/items?bbox=-3.1,51.2,-2.2,51.75&limit={limit}&f=json'
 d=js(u)
 return {'type':'FeatureCollection','features':d.get('features',[]),'numberMatched':d.get('numberMatched'),'source':u}
def bgs_geology(collection):
 u=f'https://ogcapi.bgs.ac.uk/collections/{collection}/items?bbox=-3.1,51.2,-2.2,51.75&limit=1000&f=json'
 d=js(u); out=[]
 for f in d.get('features',[]):
  try: geom=mapping(shape(f['geometry']).simplify(.002,preserve_topology=True))
  except Exception: continue
  p=f.get('properties',{})
  keep={k:p.get(k) for k in ('lex_d','lex_rcs_d','rcs_d','rock_d','max_period','min_period','feature_d','version','nom_scale') if p.get(k) not in (None,' ')}
  out.append({'type':'Feature','geometry':geom,'properties':keep})
 return {'type':'FeatureCollection','features':out,'numberMatched':d.get('numberMatched'),'source':u}
def nrfa_station(station=53006):
 meta=js(f'https://nrfaapps.ceh.ac.uk/nrfa/ws/station-info?station={station}&format=json-object&fields=all')['data'][0]
 flow=js(f'https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series?station={station}&data-type=gmf&format=json-object')
 stream=flow.get('data-stream',[]); monthly=[]
 for i in range(0,len(stream)-1,2):
  if stream[i+1] is not None: monthly.append({'date':stream[i],'value':stream[i+1]})
 keep=['id','name','river','location','latitude','longitude','catchment-area','opened','closed','gdf-start-date','gdf-end-date','gdf-mean-flow','gdf-q95-flow','gdf-q70-flow','gdf-q10-flow','bfihost','bfihost19','saar','urban-extent','factors-affecting-runoff']
 return {'station':{k:meta.get(k) for k in keep},'monthly_flow_m3s':monthly,'source':f'https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series?station={station}&data-type=gmf'}
def hydro_indicators(era):
 bymonth={m:[] for m in range(1,13)}
 for r in era:
  y=int(r['date'][:4]); m=int(r['date'][5:]);
  if 1991<=y<=2020: bymonth[m].append(r)
 out=[]; rolling=[]
 for r in era:
  m=int(r['date'][5:]); base=bymonth[m]; p=r.get('precipitation_sum'); sm=r.get('soil_moisture_28_to_100cm_mean')
  pv=[x['precipitation_sum'] for x in base if x.get('precipitation_sum') is not None]; sv=[x['soil_moisture_28_to_100cm_mean'] for x in base if x.get('soil_moisture_28_to_100cm_mean') is not None]
  pz=(p-statistics.mean(pv))/statistics.stdev(pv) if p is not None and len(pv)>1 else None
  smpct=sum(v<=sm for v in sv)/len(sv)*100 if sm is not None and sv else None
  rolling.append(p); rolling=rolling[-12:]
  out.append({'date':r['date'],'precip_z':round(pz,3) if pz is not None else None,'deep_soil_percentile':round(smpct,1) if smpct is not None else None,'rain_12m_mm':round(sum(rolling),2) if len(rolling)==12 and all(v is not None for v in rolling) else None})
 return out
def ea_readings(measure,start,limit=100000):
 mid=urllib.parse.quote('http://environment.data.gov.uk/hydrology/id/measures/'+measure,safe='')
 u=f'https://environment.data.gov.uk/hydrology/data/readings.json?measure={mid}&mineq-date={start}&_limit={limit}'
 b,d=payload(u); items=d.get('items',[])
 return items,{'url':u,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b),'api_version':d.get('meta',{}).get('version','2.1.1')}
def daily_summary(items,mode='mean'):
 by={}
 for x in items:
  if x.get('value') is not None: by.setdefault(x.get('date') or x.get('dateTime','')[:10],[]).append(float(x['value']))
 out=[]
 for day,v in sorted(by.items()):
  value=sum(v)/len(v) if mode=='mean' else max(v) if mode=='max' else sum(v)
  out.append({'date':day,'value':round(value,4),'n':len(v)})
 return out
def monthly_summary(rows,mode='mean'):
 by={}
 for r in rows:
  if r.get('value') is not None: by.setdefault(r['date'][:7],[]).append(float(r['value']))
 return [{'date':k,'value':round(sum(v)/len(v) if mode=='mean' else sum(v),4) if len(v)==calendar.monthrange(int(k[:4]),int(k[5:]))[1] else None,'valid_days':len(v),'expected_days':calendar.monthrange(int(k[:4]),int(k[5:]))[1]} for k,v in sorted(by.items())]
def station_monthly(rows,key,mode='mean'):
 by={}
 for r in rows:
  if r.get(key) is not None: by.setdefault(r['date'][:7],[]).append(float(r[key]))
 return monthly_summary([{'date':r['date'],'value':r[key]} for r in rows if r.get(key) is not None],mode)
def ghcn_station():
 u='https://www.ncei.noaa.gov/access/services/data/v1?dataset=daily-summaries&stations=UKE00105675&startDate=1850-01-01&endDate='+date.today().isoformat()+'&format=json&units=metric&includeAttributes=true'
 b,d=payload(u); rows=[]
 for x in d:
  r={'date':x.get('DATE')}
  for src,dst in [('TMAX','tmax_c'),('TMIN','tmin_c'),('PRCP','precip_mm')]:
   if x.get(src) not in (None,''): r[dst]=float(x[src])
  if len(r)>1: rows.append(r)
 return rows,{'url':u,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b),'version':'GHCN-Daily access snapshot'}
def daily_rows(raw):
 out=[]
 for line in raw.splitlines():
  p=line.split()
  if len(p)<2 or len(p[0])!=10 or p[0][4]!='-': continue
  try:v=float(p[1])
  except:continue
  if v>-90:out.append({'date':p[0],'value':v})
 return out
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
def return_levels(annual,key):
 values=[r[key] for r in annual if r.get(key) is not None]
 if len(values)<30:return None
 c,loc,scale=genextreme.fit(values)
 return {'method':'GEV maximum-likelihood fit','n_years':len(values),'shape':round(float(c),5),'location':round(float(loc),4),'scale':round(float(scale),4),'levels':{str(t):round(float(genextreme.ppf(1-1/t,c,loc=loc,scale=scale)),2) for t in (20,50,100)},'warning':'Exploratory stationary return levels; parameter uncertainty and non-stationarity are not represented.'}
def fetch_bytes(u):return urlopen(Request(u,headers=HEAD),timeout=120).read()
def isd_lite_station(usaf,start,end,name,lat,lon):
 def one(y):
  u=f'https://www.ncei.noaa.gov/pub/data/noaa/isd-lite/{y}/{usaf}-99999-{y}.gz'
  try:return y,u,fetch_bytes(u)
  except:return y,u,None
 with ThreadPoolExecutor(max_workers=10) as p:parts=list(p.map(one,range(start,end+1)))
 rows=[]; artifacts=[]
 for y,u,b in parts:
  if not b:continue
  artifacts.append({'id':f'isd_{usaf}_{y}','url':u,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b),'version':'NOAA ISD-Lite yearly artifact'})
  for line in gzip.decompress(b).decode(errors='replace').splitlines():
   q=line.split()
   if len(q)<8:continue
   try:
    yy,mm,dd,hh=map(int,q[:4]);t=int(q[4]);dew=int(q[5]);p=int(q[6]);wd=int(q[7]);ws=int(q[8]) if len(q)>8 else -9999
   except:continue
   r={'dateTime':f'{yy:04d}-{mm:02d}-{dd:02d}T{hh:02d}:00:00Z'}
   if t>-9990:r['temp_c']=t/10
   if dew>-9990:r['dewpoint_c']=dew/10
   if p>-9990:r['sea_level_pressure_hpa']=p/10
   if wd>-9990:r['wind_direction_deg']=wd
   if ws>-9990:r['wind_speed_ms']=ws/10
   if len(r)>1:rows.append(r)
 by={}
 for r in rows:by.setdefault(r['dateTime'][:10],[]).append(r)
 daily=[]
 for d,rr in sorted(by.items()):
  tv=[x['temp_c'] for x in rr if x.get('temp_c') is not None];dv=[x['dewpoint_c'] for x in rr if x.get('dewpoint_c') is not None];pv=[x['sea_level_pressure_hpa'] for x in rr if x.get('sea_level_pressure_hpa') is not None];wv=[x['wind_speed_ms'] for x in rr if x.get('wind_speed_ms') is not None]
  o={'date':d,'n_hours':len(rr)}
  if tv:o.update({'tmean_c':round(statistics.mean(tv),2),'tmax_c':max(tv),'tmin_c':min(tv)})
  if dv:o['dewpoint_c']=round(statistics.mean(dv),2)
  if pv:o['pressure_hpa']=round(statistics.mean(pv),2)
  if wv:o['wind_ms']=round(statistics.mean(wv),2)
  daily.append(o)
 return {'station':{'name':name,'usaf':usaf,'latitude':lat,'longitude':lon,'requested_coverage':f'{start}–{end}','classification':'OBSERVED · NOAA ISD-LITE · UNHOMOGENISED POINT STATION'},'daily':daily,'artifacts':artifacts}
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
def gev_nonstationary(annual,key):
 z=np.asarray([(r['year'],r[key]) for r in annual if r.get(key) is not None],float);yr=z[:,0];x=z[:,1];yc=(yr-yr.mean())/10;c0,l0,s0=genextreme.fit(x)
 def nll(q):
  c,a,b,logs=q;scale=np.exp(logs);return -np.sum(genextreme.logpdf(x,c,loc=a+b*yc,scale=scale))
 fit=minimize(nll,[c0,l0,0,np.log(s0)],method='Nelder-Mead');ll1=-fit.fun;ll0=np.sum(genextreme.logpdf(x,c0,loc=l0,scale=s0));aic0=2*3-2*ll0;aic1=2*4-2*ll1
 return {'variable':key,'stationary_aic':round(float(aic0),2),'linear_location_aic':round(float(aic1),2),'delta_aic_stationary_minus_linear':round(float(aic0-aic1),2),'location_change_per_decade':round(float(fit.x[2]),4),'preferred_by_aic':'linear-location' if aic1+2<aic0 else 'no material preference','warning':'Exploratory model comparison, not proof of a physical cause.'}
sources={
 'cet':'https://www.metoffice.gov.uk/hadobs/hadcet/data/meantemp_monthly_totals.txt',
 'cetmax':'https://www.metoffice.gov.uk/hadobs/hadcet/data/maxtemp_monthly_totals.txt',
 'cetmin':'https://www.metoffice.gov.uk/hadobs/hadcet/data/mintemp_monthly_totals.txt',
 'swep':'https://www.metoffice.gov.uk/hadobs/hadukp/data/monthly/HadSWEP_monthly_totals.txt',
 'ewp':'https://www.metoffice.gov.uk/hadobs/hadukp/data/monthly/HadEWP_monthly_totals.txt',
 'nao':'https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table'}
with ThreadPoolExecutor(max_workers=2) as p:
 raw=dict(zip(sources,p.map(text,sources.values())))
 old_path=OUT/'bristol-climate.json'
 old=json.loads(old_path.read_text()) if old_path.exists() else {}
 cached=old.get('series',{}).get('bristol_reanalysis',[])
 # Refresh every run; an old cache must never permanently freeze reanalysis.
 years=[(y,min(y+9,date.today().year)) for y in range(1950,date.today().year+1,10)]
 chunks=list(p.map(era_chunk,years,chunksize=1)); era=monthly_era(chunks)
with ThreadPoolExecutor(max_workers=6) as p:
 boreholes,landslides,hydrology,bedrock,superficial,faults=list(p.map(lambda x:x(),[lambda:bgs_features('onshoreboreholeindex'),lambda:bgs_features('landslideindex'),nrfa_station,lambda:bgs_geology('bgsgeology625kbedrock'),lambda:bgs_geology('bgsgeology625ksuperficial'),lambda:bgs_geology('bgsgeology625kfaults')]))
doc={'schema_version':'2.0.0','updated_at':datetime.now(timezone.utc).isoformat(),'scope':{'local':'Bristol grid point 51.4545N, 2.5879W','regional':'South West England & Wales','national':'England & Wales','long_record':'Central England Temperature region includes Bristol as a boundary vertex','geology_bbox':'-3.1, 51.2, -2.2, 51.75'},'series':{'cet_mean_c':matrix(raw['cet']),'swep_rain_mm':matrix(raw['swep']),'ewp_rain_mm':matrix(raw['ewp']),'nao_index':matrix(raw['nao'],1950),'bristol_reanalysis':era,'derived_hydroclimate':hydro_indicators(era),'frome_monthly_flow_m3s':hydrology['monthly_flow_m3s']},'geology':{'boreholes':boreholes,'landslides':landslides},'hydrology':hydrology,'sources':[{'publisher':'Met Office Hadley Centre','dataset':'HadCET v2 monthly mean temperature','coverage':'1659–present','classification':'OBSERVED · homogenised instrumental regional series','url':sources['cet'],'licence':'UK Open Government Licence','caveat':'Regional Central England series; not a Bristol station record.'},{'publisher':'Met Office Hadley Centre','dataset':'HadUKP South West England & Wales precipitation','coverage':'1873–present','classification':'OBSERVED · homogenised regional instrumental series','url':sources['swep'],'licence':'UK Open Government Licence','caveat':'Area average; recent months may be revised after quality control.'},{'publisher':'Met Office Hadley Centre','dataset':'HadUKP England & Wales precipitation','coverage':'1766–present','classification':'OBSERVED · homogenised national instrumental series','url':sources['ewp'],'licence':'UK Open Government Licence','caveat':'Area average, not local rainfall.'},{'publisher':'NOAA Climate Prediction Center','dataset':'Monthly NAO index','coverage':'1950–present','classification':'OBSERVED/DERIVED · standardised teleconnection index','url':sources['nao'],'licence':'US Government public data','caveat':'Statistical circulation index; association is not deterministic attribution.'},{'publisher':'ECMWF Copernicus via Open-Meteo','dataset':'ERA5 / ERA5-Land historical reanalysis for Bristol grid point','coverage':'1950–present','classification':'REANALYSIS · model–observation synthesis','url':'https://open-meteo.com/en/docs/historical-weather-api','licence':'CC BY 4.0 / Copernicus terms','caveat':'Grid-cell estimate, not a station observation; early-period uncertainty is larger.'},{'publisher':'UKCEH National River Flow Archive','dataset':'Frome (Bristol) at Frenchay · station 53006','coverage':hydrology['station'].get('gdf-start-date','1961')+'–present','classification':'OBSERVED · quality-controlled gauged flow','url':'https://nrfaapps.ceh.ac.uk/nrfa/nrfa-api.html','licence':'UKCEH terms','caveat':'Catchment outlet measurement; not groundwater level or city-wide runoff.'},{'publisher':'British Geological Survey','dataset':'Onshore Borehole Index and National Landslide Database index','coverage':'current registry snapshot','classification':'OBSERVED INDEX RECORDS','url':'https://www.bgs.ac.uk/technologies/web-services/api/','licence':'Open Government Licence','caveat':'Presence records and locations; not complete subsurface logs, hazard probabilities or absence evidence.'},{'publisher':'Zenodo / UK research community','dataset':'Weather types and large-scale climate drivers (D2.3)','coverage':'ERA5 1959–2021','classification':'RESEARCH DERIVATIVE · not merged','url':'https://zenodo.org/records/16984885','licence':'Record-specific','caveat':'Relevant methodological context, but not a Bristol ground observation. Kept separate from operational series.'}],'methods':{'baseline':'1991–2020','aggregation':'Daily reanalysis values aggregated to calendar months; means for state variables, sums for precipitation and reference evapotranspiration.','recharge_proxy':'max(precipitation − FAO-56 reference evapotranspiration, 0); screening proxy only, not groundwater recharge.','derived':'Monthly precipitation z-score and deep-soil percentile use the 1991–2020 calendar-month distribution. They are derived indicators, not measurements.','nao':'Pearson correlations are descriptive associations using overlapping monthly or winter observations; no causal attribution.'},'limitations':['CET, HadUKP and the Bristol reanalysis grid cell represent different spatial supports and are never merged into one record.','BGS borehole and landslide points are index records; completeness varies and absence does not mean no feature or hazard.','ERA5-Land soil moisture is modelled grid-cell water content, not an in-situ measurement.','Monthly means can conceal short-duration rainfall and geotechnical extremes.','No reconstruction is currently spliced into the observed record. Any future reconstruction must remain visibly separate.']}
felton_measure='12011f07-3dcc-47d8-97cc-ed3f2727b184-gw-logged-i-subdaily-mAOD-qualified'
doc['series']['cet_max_c']=matrix(raw['cetmax'])
doc['series']['cet_min_c']=matrix(raw['cetmin'])
rain_measure='e584b095-fa71-4106-b877-bfdd9a248d29-rainfall-t-900-mm-qualified'
rain_daily_measure='e584b095-fa71-4106-b877-bfdd9a248d29-rainfall-t-86400-mm-qualified'
flow_measure='f4a37857-9d23-4ab2-af0d-81171ad0d52a-flow-i-900-m3s-qualified'
with ThreadPoolExecutor(max_workers=4) as p:
 futures=[p.submit(ea_readings,felton_measure,'2025-01-01'),p.submit(ea_readings,rain_measure,'2025-01-01'),p.submit(ea_readings,rain_daily_measure,'2000-01-01'),p.submit(ea_readings,flow_measure,'2025-01-01'),p.submit(ghcn_station)]
 groundwater,storm_rain,rain_daily,flow_15m,station_daily=[f.result() for f in futures]
doc['schema_version']='3.1.0'
doc['geology'].update({'bedrock':bedrock,'superficial':superficial,'faults':faults})
gw_daily=daily_summary(groundwater[0],'mean'); rain_daily_rows=daily_summary(rain_daily[0],'sum')
doc['series'].update({
 'long_ashton_station_daily':station_daily[0],
 'long_ashton_tmean_monthly_c':station_monthly([{'date':r['date'],'tmean_c':(r['tmax_c']+r['tmin_c'])/2} for r in station_daily[0] if r.get('tmax_c') is not None and r.get('tmin_c') is not None],'tmean_c','mean'),
 'long_ashton_rain_monthly_mm':station_monthly(station_daily[0],'precip_mm','sum'),
 'felton_groundwater_daily_mAOD':gw_daily,
 'felton_groundwater_monthly_mAOD':monthly_summary(gw_daily,'mean'),
 'clifton_rain_daily_mm':rain_daily_rows,
 'clifton_rain_monthly_mm':monthly_summary(rain_daily_rows,'sum'),
 'clifton_rain_15m_extremes_mm':sorted([{'dateTime':x.get('dateTime'),'value':float(x['value']),'quality':x.get('quality')} for x in storm_rain[0] if x.get('value') is not None],key=lambda x:x['value'],reverse=True)[:500],
 'frenchay_flow_15m_extremes_m3s':sorted([{'dateTime':x.get('dateTime'),'value':float(x['value']),'quality':x.get('quality')} for x in flow_15m[0] if x.get('value') is not None],key=lambda x:x['value'],reverse=True)[:500]
})
doc['hydrology']['groundwater_station']={'name':'Felton','latitude':51.388347,'longitude':-2.688753,'opened':'1982-06-04','borehole_depth_m':157.6,'aquifer':'Carboniferous Limestone','measure':felton_measure,'unit':'mAOD','sampling':'15-minute qualified observations; displayed as daily means'}
doc['hydrology']['rain_station']={'name':'Clifton Oakfield','latitude':51.460499,'longitude':-2.610283,'measure':rain_measure,'unit':'mm','sampling':'15-minute totals plus qualified daily totals'}
doc['station_observations']={'station':'Long Ashton','station_id':'UKE00105675','latitude':51.43,'longitude':-2.67,'elevation_m':51,'coverage':'2000-12-31–2002-06-29','warning':'Genuine near-Bristol station observations, but too short for long-term trend inference. Never spliced into reanalysis.'}
doc['variables']={
 'cet_mean_c':{'label':'Central England monthly mean temperature','unit':'°C','support':'regional','frequency':'monthly','missing':'source sentinel values below −90 excluded'},
 'swep_rain_mm':{'label':'South West England & Wales precipitation','unit':'mm/month','support':'area average','frequency':'monthly','missing':'source sentinel values below −90 excluded'},
 'ewp_rain_mm':{'label':'England & Wales precipitation','unit':'mm/month','support':'area average','frequency':'monthly','missing':'source sentinel values below −90 excluded'},
 'nao_index':{'label':'North Atlantic Oscillation index','unit':'standardised index','support':'North Atlantic circulation','frequency':'monthly','missing':'unavailable months omitted'},
 'bristol_reanalysis':{'label':'ERA5/ERA5-Land Bristol grid cell','unit':'variable-specific: °C, %, hPa, mm, m³/m³, km/h','support':'model grid cell','frequency':'daily aggregated to monthly','missing':'null API values omitted; aggregation n varies'},
 'frome_monthly_flow_m3s':{'label':'Frome at Frenchay gauged mean flow','unit':'m³/s','support':'Frome catchment outlet','frequency':'monthly','missing':'NRFA missing values omitted'},
 'long_ashton_station_daily':{'label':'Long Ashton GHCN-Daily station','unit':'°C and mm/day','support':'point station','frequency':'daily','missing':'variables absent from a day remain absent; record is short'},
 'felton_groundwater_daily_mAOD':{'label':'Felton groundwater level','unit':'metres above Ordnance Datum','support':'157.6 m borehole, Carboniferous Limestone','frequency':'15-minute qualified values aggregated to daily mean','missing':'days without readings omitted; n retained'},
 'clifton_rain_15m_extremes_mm':{'label':'Clifton Oakfield rainfall','unit':'mm/15 min','support':'point gauge','frequency':'15-minute','missing':'qualified API values only; dashboard retains top 500 events'},
 'frenchay_flow_15m_extremes_m3s':{'label':'Frenchay river flow','unit':'m³/s instantaneous','support':'gauging station','frequency':'15-minute','missing':'qualified API values only; dashboard retains top 500 events'}
}
doc['variables'].update({
 'cet_max_c':{'label':'Central England monthly mean daily maximum temperature','unit':'°C','support':'regional Central England series','frequency':'monthly from 1878','missing':'source sentinel values below −90 excluded'},
 'cet_min_c':{'label':'Central England monthly mean daily minimum temperature','unit':'°C','support':'regional Central England series','frequency':'monthly from 1878','missing':'source sentinel values below −90 excluded'}
})
artifact_meta={'long_ashton_ghcn':station_daily[1],'felton_groundwater':groundwater[1],'clifton_rain_15m':storm_rain[1],'clifton_rain_daily':rain_daily[1],'frenchay_flow_15m':flow_15m[1]}
for k,v in raw.items(): artifact_meta[k]={'url':sources[k],'sha256':hashlib.sha256(v.encode()).hexdigest(),'bytes':len(v.encode()),'version':'publisher current snapshot'}
artifact_meta['bristol_reanalysis']={'url':'https://archive-api.open-meteo.com/v1/archive','sha256':hashlib.sha256(json.dumps(chunks,separators=(',',':')).encode()).hexdigest(),'bytes':len(json.dumps(chunks,separators=(',',':')).encode()),'version':'ERA5/ERA5-Land API snapshot'}
doc['artifacts']=[{'id':k,**v} for k,v in artifact_meta.items()]
doc['geology']['detailed_wms']={'title':'BGS Geology 50K','scale':'1:50,000','service':'https://map.bgs.ac.uk/arcgis/services/BGS_Detailed_Geology/MapServer/WMSServer','access':'view service; underlying vector data are not redistributed','licence':'BGS WMS terms / personal non-commercial viewing','limitation':'Map interpretation varies with source sheet age and scale; consult BGS for licensed analytical data.'}
doc['evidence_gaps']=[{'item':'Time-resolved slope reaction','status':'NOT AVAILABLE AS AN OPEN MATCHED SERIES','reason':'The BGS landslide index provides occurrence locations, not continuous displacement or pore-pressure monitoring. No synthetic response is shown.'},{'item':'Long Bristol station climatology','status':'PARTIAL','reason':'Long Ashton GHCN observations are genuine but only cover 2000–2002. Long trends therefore use clearly separated regional records and reanalysis.'},{'item':'Detailed 1:50k geology','status':'VIEW SERVICE','reason':'Official BGS WMS is displayed; analytical vector redistribution requires an appropriate BGS licence.'}]
doc['methods'].update({'inference':'Pearson and Spearman coefficients are accompanied by 95% intervals and p-values. Effective sample size uses lag-1 autocorrelation (Bretherton-type adjustment).','seasonality':'Monthly climatological anomalies use 1991–2020 calendar-month means; this removes the mean seasonal cycle before inference.','multiple_testing':'Sensitivity grid results are controlled with Benjamini–Hochberg false-discovery rate q=0.05.','sensitivity':'The workbench repeats comparisons across 1/3/12-month smoothing and 0/1/3/6/12-month lags; coefficients are descriptive, not causal.','high_frequency':'Qualified 15-minute rain and flow readings support event screening. Groundwater is aggregated to daily means with observation count retained.'})
doc['limitations'].extend(['The Long Ashton station series is genuine but short and cannot support a secular Bristol trend.','The detailed BGS 1:50,000 layer is a viewing service, not redistributed analytical geometry.','No openly matched slope-displacement series is available; rainfall and flow peaks are potential forcing indicators, not observed slope reactions.'])
doc['sources'].extend([{'publisher':'NOAA NCEI','dataset':'GHCN-Daily · Long Ashton UKE00105675','coverage':'2000–2002','classification':'OBSERVED · POINT STATION','url':station_daily[1]['url'],'licence':'NOAA public data','caveat':'Short and incomplete; shown independently and not used as a long climatology.'},{'publisher':'Environment Agency','dataset':'Felton qualified groundwater level · Hydrology API 2.1.1','coverage':'1982–present; current dashboard window from 2025','classification':'OBSERVED · BOREHOLE','url':groundwater[1]['url'],'licence':'Open Government Licence 3.0','caveat':'A single Carboniferous Limestone borehole south-west of Bristol; not a regional water table.'},{'publisher':'Environment Agency','dataset':'Clifton Oakfield rainfall and Frenchay flow · qualified 15-minute series','coverage':'current high-frequency monitoring window','classification':'OBSERVED · HIGH FREQUENCY','url':storm_rain[1]['url'],'licence':'Open Government Licence 3.0','caveat':'Event screening only; a rainfall peak is not evidence of a landslide response.'},{'publisher':'British Geological Survey','dataset':'BGS Geology 50K WMS','coverage':'Great Britain','classification':'AUTHORITATIVE MAP VIEW','url':'https://www.bgs.ac.uk/technologies/web-map-services-wms/web-map-services-geology-50k/','licence':'BGS WMS terms','caveat':'View-only local detail; vector redistribution and analysis require appropriate licensing.'}])
doc['sources'].append({'publisher':'Met Office Hadley Centre','dataset':'HadCET monthly maximum and minimum temperature','coverage':'1878–present','classification':'OBSERVED · HOMOGENISED REGIONAL SERIES','url':'https://www.metoffice.gov.uk/hadobs/hadcet/','licence':'UK Open Government Licence','caveat':'Regional Central England series; not a Bristol point station. Mean temperature extends further back to 1659.'})
daily_sources={
 'cet_daily_mean':'https://www.metoffice.gov.uk/hadobs/hadcet/data/meantemp_daily_totals.txt',
 'cet_daily_max':'https://www.metoffice.gov.uk/hadobs/hadcet/data/maxtemp_daily_totals.txt',
 'cet_daily_min':'https://www.metoffice.gov.uk/hadobs/hadcet/data/mintemp_daily_totals.txt',
 'ewp_daily':'https://www.metoffice.gov.uk/hadobs/hadukp/data/daily/HadEWP_daily_totals.txt',
 'swep_daily':'https://www.metoffice.gov.uk/hadobs/hadukp/data/daily/HadSWEP_daily_totals.txt'}
with ThreadPoolExecutor(max_workers=5) as p:daily_raw=dict(zip(daily_sources,p.map(text,daily_sources.values())))
daily={k:daily_rows(v) for k,v in daily_raw.items()}
annual_swep,p95_swep=annual_extremes(daily['swep_daily'],daily['cet_daily_mean'],daily['cet_daily_max'],daily['cet_daily_min'])
annual_ewp,p95_ewp=annual_extremes(daily['ewp_daily'],daily['cet_daily_mean'],daily['cet_daily_max'],daily['cet_daily_min'])
doc['schema_version']='5.1.0'
doc['extremes']={'baseline':'1961–1990 wet-day distribution','swep_annual':annual_swep,'ewp_annual':annual_ewp,'swep_wet_day_p95_mm':round(p95_swep,2),'ewp_wet_day_p95_mm':round(p95_ewp,2),'return_levels':{'swep_rx1day':return_levels(annual_swep,'rx1day_mm'),'swep_rx5day':return_levels(annual_swep,'rx5day_mm'),'cet_txx':return_levels(annual_swep,'txx_c')},'definitions':{'RX1day':'Annual maximum one-day precipitation','RX5day':'Annual maximum consecutive five-day precipitation','CDD':'Maximum consecutive days with precipitation <1 mm','R10/R20':'Days with precipitation ≥10/20 mm','R95p':'Annual precipitation above the 1961–1990 wet-day 95th percentile','SU25':'Days with daily maximum temperature >25 °C','FD0':'Days with daily minimum temperature <0 °C','ID0':'Days with daily maximum temperature <0 °C','TR20':'Nights with daily minimum temperature >20 °C','TXx/TNn':'Annual maximum Tmax / minimum Tmin','DTR':'Mean daily Tmax − Tmin'}}
daily_doc={'schema_version':'5.0.0','updated_at':doc['updated_at'],'classification':'OBSERVED · DAILY REGIONAL SERIES; spatial supports remain separate','series':daily,'variables':{'cet_daily_mean':{'unit':'°C','coverage':'1772–present','support':'Central England region'},'cet_daily_max':{'unit':'°C','coverage':'1878–present','support':'Central England region'},'cet_daily_min':{'unit':'°C','coverage':'1878–present','support':'Central England region'},'ewp_daily':{'unit':'mm/day','coverage':'1931–present','support':'England & Wales area average'},'swep_daily':{'unit':'mm/day','coverage':'1931–present','support':'South West England & Wales area average'}},'sources':daily_sources}
daily_doc['schema_version']='5.1.0'
for k,v in daily_raw.items():doc['artifacts'].append({'id':k,'url':daily_sources[k],'sha256':hashlib.sha256(v.encode()).hexdigest(),'bytes':len(v.encode()),'version':'publisher current snapshot'})
doc['sources'].extend([{'publisher':'Met Office Hadley Centre','dataset':'HadCET daily mean / maximum / minimum temperature','coverage':'1772–present (mean); 1878–present (maximum/minimum)','classification':'OBSERVED · HOMOGENISED REGIONAL DAILY SERIES','url':'https://www.metoffice.gov.uk/hadobs/hadcet/data/download.html','licence':'UK Open Government Licence','caveat':'Central England regional evidence; not a Bristol point station.'},{'publisher':'Met Office Hadley Centre','dataset':'HadUKP daily precipitation · South West England & Wales / England & Wales','coverage':'1931–present','classification':'OBSERVED · HOMOGENISED REGIONAL DAILY SERIES','url':'https://www.metoffice.gov.uk/hadobs/hadukp/data/download.html','licence':'UK Open Government Licence','caveat':'Area averages; suitable for regional extremes context, not a local gauge.'},{'publisher':'Met Office','dataset':'HadUK-Grid v1.3.2.0','coverage':'variable-dependent; monthly precipitation from 1836, daily precipitation from 1891','classification':'AUTHORITATIVE GRIDDED OBSERVATIONS · ACCESS-GATED, NOT IMPORTED','url':'https://www.metoffice.gov.uk/research/climate/maps-and-data/data/haduk-grid/datasets','licence':'UK Open Government Licence','caveat':'The CEDA archive requires authenticated access. It is registered here but no grid values are shown as imported.'},{'publisher':'Copernicus Climate Change Service','dataset':'E-OBS gridded observations for Europe','coverage':'1950–present','classification':'AUTHORITATIVE GRIDDED OBSERVATIONS · NOT IMPORTED','url':'https://cds.climate.copernicus.eu/datasets/insitu-gridded-observations-europe','licence':'Copernicus licence','caveat':'CDS account/API access is required; DOI 10.24381/cds.151d3ec6.'},{'publisher':'Met Office','dataset':'UKCP18 climate projections','coverage':'baseline and projections to 2100','classification':'SCENARIO-CONDITIONED PROJECTIONS · NOT IMPORTED','url':'https://www.metoffice.gov.uk/research/approach/collaboration/ukcp','licence':'UKCP terms','caveat':'Account-mediated datasets. No projection is presented as an observation or deterministic forecast.'}])
doc['publications']=[
 {'doi':'10.1029/2023EF004073','year':2024,'scope':'Bristol pluvial flooding','evidence':'Peer-reviewed event and urban flood analysis','role':'Supports interpretation of intense rainfall and urban surface-water processes; no values merged into observations.'},
 {'doi':'10.1029/2023WR035533','year':2024,'scope':'Bristol and Bath future urban flooding','evidence':'Peer-reviewed UKCP Local modelling','role':'Scenario context for future pluvial flooding; model output remains separate.'},
 {'doi':'10.1080/02626667.2022.2038791','year':2022,'scope':'NAO, rainfall and streamflow','evidence':'Peer-reviewed hydroclimatology','role':'Supports the NAO comparison design; correlation is not causal attribution.'},
 {'doi':'10.1093/biohorizons/hzp009','year':2009,'scope':'Bristol Avon hydrology and ecology','evidence':'Peer-reviewed regional study','role':'Catchment context only; no digitised values imported.'},
 {'doi':'10.1007/s10584-014-1101-8','year':2014,'scope':'South West England climate variability','evidence':'Peer-reviewed regional climate study','role':'Regional interpretation and variable selection.'},
 {'doi':'10.3390/su12083233','year':2020,'scope':'Bristol resilience under multiple models','evidence':'Peer-reviewed systems study','role':'Resilience context; not a substitute for observed climate data.'}]
with ThreadPoolExecutor(max_workers=3) as p:
 local_segments=list(p.map(lambda f:f(),[
  lambda:isd_lite_station('036280',1937,2017,'Filton',51.517,-2.583),
  lambda:isd_lite_station('037260',1981,1997,'Bristol Weather Centre',51.467,-2.600),
  lambda:isd_lite_station('037243',1988,2025,'Bristol Airport / Lulsgate',51.383,-2.719)]))
doc['local_stations']=[{'station':s['station'],'daily':s['daily']} for s in local_segments]
for s in local_segments:doc['artifacts'].extend(s['artifacts'])
doc['sources'].append({'publisher':'NOAA National Centers for Environmental Information','dataset':'Integrated Surface Database Lite · Filton, Bristol Weather Centre and Bristol Airport','coverage':'station-dependent, 1937–2025','classification':'OBSERVED · SUBDAILY POINT STATIONS · UNHOMOGENISED','url':'https://www.ncei.noaa.gov/products/land-based-station/integrated-surface-database','licence':'NOAA public data','caveat':'Separate station segments with different locations, instruments and completeness. They are not spliced into a homogenised Bristol climate series.'})
def local_annual(segment):
 by={}
 for r in segment['daily']:
  if r.get('tmean_c') is not None:by.setdefault(int(r['date'][:4]),[]).append(r)
 return [{'year':y,'tmean_c':round(statistics.mean(x['tmean_c'] for x in rr),2),'txx_c':max(x['tmax_c'] for x in rr if x.get('tmax_c') is not None),'tnn_c':min(x['tmin_c'] for x in rr if x.get('tmin_c') is not None),'days':len(rr)} for y,rr in sorted(by.items()) if len(rr)>=300]
doc['local_station_annual']={s['station']['usaf']:local_annual(s) for s in local_segments}
doc['inference_v5']={'bootstrap':[gev_bootstrap(annual_swep,'rx1day_mm'),gev_bootstrap(annual_swep,'rx5day_mm'),gev_bootstrap(annual_swep,'txx_c')],'trends':[trend_test(annual_swep,k) for k in ('rx1day_mm','rx5day_mm','cdd_days','wet_days','summer_days','frost_days','txx_c','tnn_c')],'nonstationarity':[gev_nonstationary(annual_swep,k) for k in ('rx1day_mm','rx5day_mm','txx_c')],'interpretation_boundary':{'association':'Statistical co-variation, trend or change-point evidence. It does not establish a physical mechanism.','attribution':'Requires a declared causal design, counterfactual or process model and is not performed automatically.','geological_causality':'Requires local lithology, permeability, slope, pore pressure or displacement observations plus temporal alignment. Climate forcing alone is insufficient.','prohibited_claim':'No chart may state that rainfall or temperature caused a geological response unless the response itself was observed and the causal design was declared.'}}
def sensitivity_table():
 out=[]
 for spatial,rows in [('South West England & Wales',daily['swep_daily']),('England & Wales',daily['ewp_daily'])]:
  for baseline in ((1961,1990),(1981,2010),(1991,2020)):
   for wet_threshold in (.5,1,2):
    wet=[r['value'] for r in rows if baseline[0]<=int(r['date'][:4])<=baseline[1] and r['value']>=wet_threshold]
    p95=float(np.percentile(wet,95));annual,_=annual_extremes(rows,daily['cet_daily_mean'],daily['cet_daily_max'],daily['cet_daily_min']);valid=[r for r in annual if r.get('rx1day_mm') is not None]
    tr=trend_test(valid,'rx1day_mm');out.append({'spatial_support':spatial,'baseline':f'{baseline[0]}–{baseline[1]}','wet_day_threshold_mm':wet_threshold,'wet_day_p95_mm':round(p95,2),'rx1day_sen_per_decade':tr['theil_sen_per_decade'],'rx1day_mk_p':tr['mann_kendall_p']})
 return out
doc['sensitivity_v5']=sensitivity_table()
doc['evidence_gaps'].append({'item':'HadUK-Grid 1 km Bristol cells','status':'AUTHENTICATED CEDA ACCESS REQUIRED','reason':'Official v1.3.2.ceda, DOI 10.5285/789b3065d74a4c948ab05d33556c86d0, is 463 GB and requires CEDA registration/login. No values are represented as imported until authenticated retrieval is completed.'})
all_dates=[r.get('date') or r.get('dateTime','')[:10] for series in doc['series'].values() for r in series if isinstance(r,dict) and (r.get('date') or r.get('dateTime'))]
doc['summary']={'records':sum(len(x) for x in doc['series'].values())+sum(len(x) for x in daily.values()),'monthly_records':sum(len(x) for x in doc['series'].values()),'daily_records':sum(len(x) for x in daily.values()),'cet_start':doc['series']['cet_mean_c'][0]['date'],'swep_start':doc['series']['swep_rain_mm'][0]['date'],'ewp_start':doc['series']['ewp_rain_mm'][0]['date'],'reanalysis_start':doc['series']['bristol_reanalysis'][0]['date'],'latest':max(all_dates)}
from bristol_quality import finish_document
finish_document(doc,groundwater,storm_rain,flow_15m)
daily_doc['updated_at']=doc['updated_at']
daily_bytes=json.dumps(daily_doc,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()
daily_temp=OUT/'bristol-climate-daily.json.tmp'
daily_temp.write_bytes(daily_bytes)
daily_temp.replace(OUT/'bristol-climate-daily.json')
(OUT/'bristol-climate-daily.json.sha256').write_text(hashlib.sha256(daily_bytes).hexdigest()+'  bristol-climate-daily.json\n')
rawb=json.dumps(doc,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode(); temp=OUT/'bristol-climate.json.tmp';temp.write_bytes(rawb);temp.replace(OUT/'bristol-climate.json'); (OUT/'bristol-climate.json.sha256').write_text(hashlib.sha256(rawb).hexdigest()+'  bristol-climate.json\n'); print(doc['summary'])
