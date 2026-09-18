import csv
from pathlib import Path
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from check_energy_operators import check

class EnergyLedger(unittest.TestCase):
    def fixture(self,p):
        keys=['t','dt','stage','species','before_J','after_J','delta_J',*[f'term{i}_J' for i in range(1,6)],'scale_J','residual_J']
        records=[]
        for t in (0,1):
            for stage in (1,2):
                for species in range(1,8):
                    values=[t,1,stage,species,10,11,1,1,0,0,0,0,11,0]
                    records.append(dict(zip(keys,values)))
        self.save(p,records)
        (p/'continuity-r00000000.csv').write_text('t,dt,species\n0,1,1\n1,1,1\n')
        return records
    def save(self,p,records):
        with (p/'energy-operators-r00000000.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    def test_complete_and_missing_later_step(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);r=self.fixture(p);self.assertTrue(check(p,1)['passed'])
            self.save(p,r[:-1])
            with self.assertRaisesRegex(ValueError,'Missing'):check(p,1)
    def test_faults_cannot_hide_in_logged_residual(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);r=self.fixture(p);r[-1]['term3_J']=.01;self.save(p,r)
            self.assertFalse(check(p,1)['passed'])
            r=self.fixture(p);r[0]['delta_J']=float('nan');self.save(p,r)
            with self.assertRaisesRegex(ValueError,'Nonfinite'):check(p,1)
    def test_rank_clock_scale_and_duplicate(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);r=self.fixture(p)
            with self.assertRaisesRegex(ValueError,'inventory'):check(p,2)
            r[0]['scale_J']=1;self.save(p,r)
            with self.assertRaisesRegex(ValueError,'scale'):check(p,1)
            r=self.fixture(p);r.append(r[-1]);self.save(p,r)
            with self.assertRaisesRegex(ValueError,'duplicate'):check(p,1)
if __name__=='__main__':unittest.main()
