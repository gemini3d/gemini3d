"""Joint evidence deletion must not turn unexecuted steps into passes."""
from pathlib import Path
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from check_sources import check_clock_inventory


class ClockInventory(unittest.TestCase):
    def test_complete_run_and_float64_clock_roundoff(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'inputs').mkdir()
            (p/'inputs/config.nml').write_text('&base\n tdur = 3D0, ! fixed duration\n/\n')
            check_clock_inventory(p,{(0.,1.),(1.,1.),(2.,1.)})
            check_clock_inventory(p,{(0.,1.),(1.,1.),(2.+1e-15,1.)})

    def test_missing_first_middle_and_final_steps(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'inputs').mkdir()
            (p/'inputs/config.nml').write_text('tdur = 3\n')
            for clocks,message in [({(1.,1.),(2.,1.)},'begin'),
                                    ({(0.,1.),(2.,1.)},'Noncontiguous'),
                                    ({(0.,1.),(1.,1.)},'duration')]:
                with self.subTest(clocks=clocks),self.assertRaisesRegex(ValueError,message):
                    check_clock_inventory(p,clocks)

    def test_overlaps_nonadvancing_clocks_and_invalid_duration(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for clocks in [{(0.,2.),(1.,1.)},{(0.,0.)},{(-1.,1.)},{(0.,float('nan'))}]:
                with self.subTest(clocks=clocks),self.assertRaises(ValueError):check_clock_inventory(p,clocks)
            (p/'inputs').mkdir()
            for value in ['nan','inf','0','-1','1,2','']:
                (p/'inputs/config.nml').write_text('tdur = '+value+'\n')
                with self.subTest(value=value),self.assertRaises(ValueError):check_clock_inventory(p,{(0.,1.)})

if __name__=='__main__':unittest.main()
