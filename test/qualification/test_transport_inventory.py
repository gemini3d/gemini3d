"""A complete surviving step cannot conceal missing numerical evidence."""
import csv
from pathlib import Path
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from check_transport import check

class TransportInventory(unittest.TestCase):
    def fixture(self,directory):
        fields='t,dt,quantity,species,axis,before,after,delta,outward,min_face,max_face,scale,residual,normalized_residual'.split(',')
        records=[dict(zip(fields,[t,1,q,s,a,100,99,-1,1,0,1,100,0,0]))
                 for t in (0,1) for q in (1,2,3)
                 for s in range(1,8 if q==3 else 7) for a in (1,2,3)]
        path=directory/'transport-r00000000.csv'
        def save(rows):
            with path.open('w') as stream:
                writer=csv.DictWriter(stream,fields);writer.writeheader();writer.writerows(rows)
        save(records)
        (directory/'continuity-r00000000.csv').write_text('t,dt\n0,1\n1,1\n')
        return records,save

    def test_entire_later_step_is_required(self):
        with tempfile.TemporaryDirectory() as d:
            directory=Path(d);records,save=self.fixture(directory)
            self.assertTrue(check(directory,layout=[1,1])['passed'])
            save([r for r in records if r['t']==0])
            with self.assertRaisesRegex(ValueError,'Transport/continuity'):check(directory,layout=[1,1])

    def test_recorded_state_change_must_match_increment(self):
        with tempfile.TemporaryDirectory() as d:
            directory=Path(d);records,save=self.fixture(directory)
            # The declared scale remains valid, and delta + outward still closes.
            records[0]['after']=50;save(records)
            self.assertFalse(check(directory,layout=[1,1])['passed'])

    def test_invalid_budgets_and_scales_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            directory=Path(d);records,save=self.fixture(directory)
            for budget in (0,-1,float('inf'),float('nan')):
                with self.subTest(budget=budget),self.assertRaisesRegex(ValueError,'budget'):
                    check(directory,budget=budget,layout=[1,1])
            records[0]['scale']=1;save(records)
            with self.assertRaisesRegex(ValueError,'scale'):check(directory,layout=[1,1])

if __name__=='__main__':unittest.main()
