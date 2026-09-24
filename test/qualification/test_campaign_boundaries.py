from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from campaign import conformal,coverage,qualify_coverage,rollout_metrics
from observation_operator import radar_los,linear_analysis
class CampaignBoundaries(unittest.TestCase):
    def test_masks_cannot_convert_missingness_into_validity(self):
        for mask in ([np.nan],[2],['valid'],[.5]):
            with self.assertRaisesRegex(ValueError,'mask'):radar_los([[1,0,0]],[[1,0,0]],[1],mask)
    def test_empty_predictions_and_observations_reject(self):
        with self.assertRaises(ValueError):rollout_metrics(np.empty((0,3)),np.empty((0,3)),np.ones(3))
        with self.assertRaises(ValueError):linear_analysis([1],[[1]],np.empty((0,1)),[],np.empty((0,0)))
    def test_finite_sample_confidence_cannot_be_replaced_with_empirical_coverage(self):
        c=conformal([.1]*19,list(range(19)),.05)
        small=qualify_coverage([0]*10,list(range(100,110)),c,list(range(19)),coverage_target=.95,width_cap=1,confidence=.95)
        self.assertEqual(small['empirical_coverage'],1)
        self.assertFalse(small['passed'])
        self.assertAlmostEqual(small['population_coverage_lower_bound'],.05**.1)
        large=qualify_coverage([0]*59,list(range(100,159)),c,list(range(19)),coverage_target=.95,width_cap=1,confidence=.95)
        self.assertTrue(large['passed'])
        wide=qualify_coverage([0]*59,list(range(100,159)),c,list(range(19)),coverage_target=.95,width_cap=.1,confidence=.95)
        self.assertFalse(wide['passed'])
        none=qualify_coverage([1]*59,list(range(100,159)),c,list(range(19)),coverage_target=.95,width_cap=1,confidence=.95)
        self.assertEqual(none['population_coverage_lower_bound'],0)
    def test_bad_calibration_inventory_and_radius(self):
        c=conformal([.1]*19,list(range(19)),.05)
        with self.assertRaises(ValueError):coverage([0],[100],c,[0]*19)
        c['radius']=-1
        with self.assertRaises(ValueError):coverage([0],[100],c,list(range(19)))
if __name__=='__main__':unittest.main()
