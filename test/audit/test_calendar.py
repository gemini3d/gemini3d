"""Independent Gregorian calendar oracle for the compiled time utilities."""
import argparse,datetime as d,json,math,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--exe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();rows=[]
def args(t):return [t.year,t.month,t.day,(t-t.replace(hour=0,minute=0,second=0,microsecond=0)).total_seconds()]
def run(mode,values,expected=None):
 r=subprocess.run([str(a.exe),mode,*map(str,values)],capture_output=True,text=True,timeout=5)
 ok=r.returncode!=0 if expected is None else r.returncode==0
 if expected is not None and ok:
  got=[float(x) for x in r.stdout.split()];ok=len(got)==len(expected) and all(math.isclose(x,y,abs_tol=1e-7,rel_tol=0) for x,y in zip(got,expected))
 rows.append({'mode':mode,'input':values,'passed':ok,'returncode':r.returncode})
for start in [d.datetime(1999,12,31,23,59,45),d.datetime(2000,2,28,23,59),d.datetime(1900,2,28,23,59),d.datetime(2026,1,31,23,59,45),d.datetime(2026,9,16,0,0,30)]:
 for seconds in [0,0.5,1,15,60,180,86400,3*86400+35]:
  end=start+d.timedelta(seconds=seconds)
  for cadence in [0.5,15,60,3600]:
   offset=math.floor(seconds/cadence)*cadence
   run('last',args(start)+args(end)+[cadence],args(start+d.timedelta(seconds=offset)))
   run('elapsed',args(start)+args(end)+[cadence],[offset])
for offset in [-172801,-86400,-60,-0.5,0,0.5,60,86400,172801]:
 start=d.datetime(2000,3,1,0,0,30)
 run('shift',args(start)+[offset],args(start+d.timedelta(seconds=offset)))
for cadence in [0,-1,'NaN',86401]:run('last',[2026,1,1,0,2026,1,2,0,cadence])
run('last',[2026,2,30,0,2026,3,1,0,60])
run('elapsed',[2026,1,2,0,2026,1,1,0,60])
result={'count':len(rows),'passed':all(x['passed'] for x in rows),'cases':rows}
a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='cases'}))
raise SystemExit(0 if result['passed'] else 1)
