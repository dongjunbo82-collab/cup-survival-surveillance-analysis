"""Verify the observed SEER*Stat aggregate export; no trend model is fitted.
Run: python verify_incidence.py [path/to/aggregate_export.txt]
Default input is the aggregate text export in this script's directory.
"""
import csv,json,sys
from pathlib import Path

base=Path(__file__).resolve().parent
src=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__import__('os').environ.get('CUP_INCIDENCE_EXPORT', str(base/'incidence_observed.txt')))
with src.open(encoding='utf-8-sig',newline='') as f:
    raw=list(csv.DictReader(f,delimiter='\t'))
annual=[]
for row in raw:
    year=row['Year of diagnosis']
    if not (year.isdigit() and 2004<=int(year)<=2022):continue
    z={k:(v if k=='Year of diagnosis' else float(v)) for k,v in row.items()}
    assert z['Population']>0 and z['Standard Error']>0
    assert z['Lower Confidence Interval']<=z['Age-Adjusted Rate']<=z['Upper Confidence Interval']
    annual.append(z)
assert [int(z['Year of diagnosis']) for z in annual]==list(range(2004,2023))
tot=next(z for z in raw if z['Year of diagnosis']=='2000-2023')
assert sum(z['Count'] for z in annual)==float(tot['Count'])
assert sum(z['Population'] for z in annual)==float(tot['Population'])
assert {z['Std Pop'] for z in annual}=={195850985}
base=Path(__import__('os').environ.get('CUP_INCIDENCE_OUTPUT', 'outputs/incidence'));base.mkdir(parents=True,exist_ok=True)
out=base/'annual_observed_rates.csv'
with out.open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(annual[0]));w.writeheader();w.writerows(annual)
summary={'years':len(annual),'cases':int(sum(z['Count'] for z in annual)),
         'population_person_years':int(sum(z['Population'] for z in annual)),
         'first_year':annual[0],'last_year':annual[-1],
         'endpoint_change_percent':100*(annual[-1]['Age-Adjusted Rate']/annual[0]['Age-Adjusted Rate']-1),
         'rate_type':'Observed age-adjusted; not delay-adjusted',
         'inference':'Exported annual SE and 95% CI only; no APC, AAPC or trend test'}
(base/'verification_results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
