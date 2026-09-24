import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from run_research import input_signature
from check_transport import check,check_continuity
from research_matrix import compare,restart_advanced,set_duration
import h5py
import numpy as np

class Profile(unittest.TestCase):
    def test_duration_rewrite_handles_indented_assignment_and_preserves_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'inputs').mkdir();cfg=p/'inputs/config.nml'
            for assignment in ['tdur = 300','  tdur = 300','\tTDUR\t= 3e2,',
                               '  tdur = 3D2','\ttdur = 3d2,']:
                with self.subTest(assignment=assignment):
                    original='&base\n! tdur = 900\n'+assignment+' ! duration\ntcfl = .5\n/\n'
                    cfg.write_text(original)
                    signature=input_signature(p)['sha256']
                    set_duration(p,120)
                    self.assertIn('= 120 ! duration\n',cfg.read_text())
                    self.assertEqual(input_signature(p)['sha256'],signature)
                    set_duration(p,300,extend=True)
                    self.assertIn('= 300 ! duration\n',cfg.read_text())
                    self.assertEqual(input_signature(p)['sha256'],signature)
                    self.assertIn('! tdur = 900\n',cfg.read_text())
                    self.assertIn('tcfl = .5\n',cfg.read_text())

    def test_duration_rewrite_rejects_missing_duplicate_and_malformed_values(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'inputs').mkdir();cfg=p/'inputs/config.nml'
            invalid=['! tdur = 300','tdur: 300','tdur = 300\n  TDUR = 120',
                     'tdur =','tdur = nope','tdur = 300, tcfl = .5',
                     'tdur = 300,,','tdur = -1','tdur = 0','tdur = nan','tdur = inf']
            for assignment in invalid:
                with self.subTest(assignment=assignment):
                    text='&base\n'+assignment+'\n/\n';cfg.write_text(text)
                    with self.assertRaisesRegex(ValueError,'tdur'):set_duration(p,120)
                    self.assertEqual(cfg.read_text(),text)
            text='&base\n  tdur = 120\n/\n'
            for duration in [0,-1,float('nan'),float('inf'),'invalid',None]:
                with self.subTest(duration=duration):
                    cfg.write_text(text)
                    with self.assertRaisesRegex(ValueError,'tdur'):set_duration(p,duration)
                    self.assertEqual(cfg.read_text(),text)
            for duration in [60,120]:
                with self.subTest(continuation=duration):
                    cfg.write_text(text)
                    with self.assertRaisesRegex(ValueError,'extend'):set_duration(p,duration,extend=True)
                    self.assertEqual(cfg.read_text(),text)

    def test_restart_requires_preserved_seed_and_a_later_checkpoint(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);first='20130220_00120.000000.h5';last='20130220_00300.000000.h5'
            self.assertFalse(restart_advanced(p,set()))
            (p/first).touch()
            self.assertFalse(restart_advanced(p,set()))
            previous={first}
            self.assertFalse(restart_advanced(p,previous))
            (p/'20130220_00060.000000.h5').touch()
            self.assertFalse(restart_advanced(p,previous))
            (p/last).touch()
            self.assertTrue(restart_advanced(p,previous))
            (p/first).unlink()
            self.assertFalse(restart_advanced(p,previous))

    def test_trajectory_requires_all_frames_and_compatible_shapes(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);a=p/'continuous';b=p/'split';a.mkdir();b.mkdir()
            budget=dict(absolute_floor=dict(ns=1e-2,Ts=1e-7,vs1=1e-7,Phi=1e-8),trajectory_relative_l2=1e-10)
            def frame(path):
                with h5py.File(path,'w') as f:
                    for key in ['ns','Ts','vs1']:f['restart_core/'+key]=np.ones((7,2,2,2))
                    f['restart_core/Phi']=np.ones((2,2,2))
            for name in ['20130220_00060.000000.h5','20130220_00120.000000.h5']:
                frame(a/name);frame(b/name)
            self.assertTrue(compare(a,b,budget)['passed'])
            self.assertFalse(restart_advanced(b,{p.name for p in b.glob('*.h5')}))
            last=b/'20130220_00120.000000.h5';last.unlink()
            self.assertFalse(compare(a,b,budget)['passed'])
            frame(last)
            with h5py.File(last,'r+') as f:
                del f['restart_core/Phi'];f['restart_core/Phi']=np.ones((1,2,2))
            self.assertFalse(compare(a,b,budget)['passed'])

    def test_only_duration_can_change(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'inputs').mkdir();cfg=p/'inputs/config.nml'
            cfg.write_text('&base\ntdur = 120 ! duration\ntcfl = .5\n/\n')
            original=input_signature(p)['sha256']
            cfg.write_text('&base\ntdur = 300 ! duration\ntcfl = .5\n/\n')
            self.assertEqual(input_signature(p)['sha256'],original)
            cfg.write_text('&base\ntdur = 300 ! duration\ntcfl = .9\n/\n')
            self.assertNotEqual(input_signature(p)['sha256'],original)
            (p/'inputs/linked').symlink_to('/dev/null')
            with self.assertRaisesRegex(ValueError,'symlinks'):input_signature(p)
    def test_budget_checks_reject_bad_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);f=p/'transport-r00000000.csv'
            fields=['t','dt','quantity','species','axis','before','after','delta','outward','min_face','max_face','scale','residual','normalized_residual']
            rows=[]
            for q in [1,2,3]:
                for species in range(1,8 if q==3 else 7):
                    for axis in [1,2,3]:
                        rows.append(dict(zip(fields,[0,1,q,species,axis,1,1,0,0,0,0,1,0,0])))
            def write():
                with f.open('w') as h:
                    w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(rows)
            write();self.assertTrue(check(p,layout=[1,1])['passed'])
            # A complete first step must not hide a missing stage in a later step.
            complete=[dict(r) for r in rows]
            rows.extend(dict(r,t=1) for r in complete[:-1]);write()
            self.assertFalse(check(p)['passed'])
            rows[:]=complete
            rows[0]['delta']=.001;write();self.assertFalse(check(p)['passed'])
            rows.pop();write();self.assertFalse(check(p)['passed'])
            rows[0]['delta']=float('nan');write()
            with self.assertRaisesRegex(ValueError,'Nonfinite'):check(p)
    def test_continuity_does_not_hide_cleanup(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);f=p/'continuity-r00000000.csv'
            fields=['t','dt','species','before','after','source','outward','cleanup','residual','scale','normalized_residual','charge_fraction']
            rows=[dict(zip(fields,[0,1,s,100,102,1,0,1,0,102,0,0])) for s in range(1,7)]
            def write():
                with f.open('w') as h:
                    w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(rows)
            write();self.assertTrue(check_continuity(p)['passed'])
            complete=[dict(r) for r in rows]
            rows.extend(dict(r,t=1) for r in complete[:-1]);write()
            self.assertFalse(check_continuity(p)['passed'])
            rows[:]=complete
            rows[0]['cleanup']=0;write();self.assertFalse(check_continuity(p)['passed'])
if __name__=='__main__':unittest.main()
