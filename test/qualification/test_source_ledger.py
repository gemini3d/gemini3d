import csv
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from check_sources import check,EXPECTED


class Sources(unittest.TestCase):
    def write(self,p):
        cols='t,dt,quantity,species,before,after,delta,production,external_production,integrated_loss,residual,scale'.split(',')
        data=[dict(zip(cols,[t,1,q,s,100,105,5,10,2,5,0,105])) for t in (0,1) for q,s in sorted(EXPECTED)]
        with (p/'sources-r00000000.csv').open('w') as f:
            w=csv.DictWriter(f,cols);w.writeheader();w.writerows(data)
        (p/'continuity-r00000000.csv').write_text('t,dt\n0,1\n1,1\n')
        (p/'temperature-floor-r00000000.csv').write_text('t,dt,stage,species,energy_added_J\n'+
            ''.join(f'{t},1,{stage},{s},0\n' for t in (0,1) for stage in (1,2) for s in range(1,8)))
        return p/'sources-r00000000.csv'

    def test_known_balance_corruption_and_per_step_coverage(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);path=self.write(p)
            self.assertTrue(check(p,1)['passed'])
            text=path.read_text();path.write_text(text.replace(',105,5,10,',',105,7,10,',1))
            self.assertFalse(check(p,1)['passed'])
            path.write_text('\n'.join(text.splitlines()[:-1])+'\n')
            with self.assertRaisesRegex(ValueError,'Missing'):check(p,1)
            path.write_text(text+text.splitlines()[1]+'\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):check(p,1)

    def test_missing_rank_nan_floor_and_clock(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);path=self.write(p)
            with self.assertRaisesRegex(ValueError,'inventory'):check(p,2)
            path.write_text(path.read_text().replace(',105,5,',',nan,5,',1))
            with self.assertRaisesRegex(ValueError,'Nonfinite'):check(p,1)
            self.write(p);floor=p/'temperature-floor-r00000000.csv';floor.write_text('t,dt,stage,species,energy_added_J\n')
            with self.assertRaisesRegex(ValueError,'Missing'):check(p,1)

    def test_four_rank_regression_preserves_all_ledger_checks(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.write(p)
            for rank in range(1,4):
                for stem in ('sources','continuity','temperature-floor'):
                    (p/f'{stem}-r{rank:08d}.csv').write_bytes((p/f'{stem}-r00000000.csv').read_bytes())
            result=check(p,4)
            self.assertTrue(result['passed'])
            self.assertEqual(result['rows'],2*4*len(EXPECTED))
            self.assertEqual(result['temperature_floor_rows'],2*4*2*7)
            last=p/'sources-r00000003.csv';original=last.read_text()
            last.write_text(original.replace(',105,5,10,',',105,7,10,',1))
            self.assertFalse(check(p,4)['passed'])
            last.write_text(original)
            continuity=p/'continuity-r00000003.csv'
            continuity.write_text('t,dt\n0,1\n')
            with self.assertRaisesRegex(ValueError,'Source/continuity'):check(p,4)
            continuity.write_bytes((p/'continuity-r00000000.csv').read_bytes())
            floor=p/'temperature-floor-r00000003.csv';floor_text=floor.read_text()
            floor.write_text('\n'.join(floor_text.splitlines()[:-1])+'\n')
            with self.assertRaisesRegex(ValueError,'Missing'):check(p,4)
            floor.write_text(floor_text);last.unlink()
            with self.assertRaisesRegex(ValueError,'inventory'):check(p,4)
            last.write_text(original)
            (p/'sources-r00000004.csv').write_text(original)
            with self.assertRaisesRegex(ValueError,'inventory'):check(p,4)


if __name__=='__main__':unittest.main()
