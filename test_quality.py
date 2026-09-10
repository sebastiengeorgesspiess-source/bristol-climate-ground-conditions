import ast
import calendar
from datetime import date, timedelta
from pathlib import Path

src = ast.parse(Path(__file__).with_name('update_climate.py').read_text())
safe = ast.Module(body=[n for n in src.body if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef))], type_ignores=[])
ns = {}
exec(compile(safe, 'functions-only', 'exec'), ns)
assert len(ns['matrix']('2026 1 2 3 4 5 6 7 8', 1950)) == 8
assert ns['matrix']('2025 -99 nan inf 1',1950)==[{'date':'2025-04','value':1.0}]
rows=[{'date':f'2025-01-{i:02d}','value':1} for i in range(1,31)]
assert ns['monthly_summary'](rows,'sum')[0]['value'] is None
rows.append({'date':'2025-01-31','value':1})
assert ns['monthly_summary'](rows,'sum')[0]['value']==31
days=[{'date':(date(2025,1,1)+timedelta(days=i)).isoformat(),'value':1} for i in range(365)]
partial=[{'date':'2026-01-01','value':999}]
annual,_=ns['annual_extremes'](days+partial,days,days,days)
assert len(annual)==1 and annual[0]['rx5day_mm']==5
bad,_=ns['annual_extremes'](days[:-1],days,days,days)
assert 'rx1day_mm' not in bad[0]
chunks=[{'time':[r['date'] for r in rows], 'temperature_2m_mean':[1]*31, 'precipitation_sum':[1]*30+[None]}]
monthly=ns['monthly_era'](chunks)[0]
assert monthly['temperature_2m_mean']==1 and monthly['precipitation_sum'] is None
print('PASS: partial NAO year, missing sentinels, complete months, per-variable coverage, full annual extremes, partial-year exclusion')
