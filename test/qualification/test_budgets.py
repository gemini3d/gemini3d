"""Manufactured curvilinear control volume, explicit flux and clipping checks."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
p=Path(__file__).resolve().parents[2]/'scripts/qualification/budgets.py'
s=importlib.util.spec_from_file_location('budgets',p);b=importlib.util.module_from_spec(s);s.loader.exec_module(b)
class Budgets(unittest.TestCase):
    def test_physical_volume_and_units(self):
        shape=(2,2,1);one=np.ones(shape);mask=np.ones(shape,dtype=bool);mask[1,1,0]=False
        vol=b.cell_volumes(2*one,3*one,4*one,[1,2],[3,4],[5],mask)
        self.assertEqual(vol.sum(),24*5*(1*3+1*4+2*3))
        n=np.ones(shape+(2,))*np.array([2,2]);t=np.full_like(n,10);v=np.ones(n.shape+(3,))
        n[~mask]=np.nan;t[~mask]=np.nan;v[~mask]=np.nan
        out=b.integrals(n,v,t,vol,[2,4],[1,-1],[5/3,5/3],1)
        np.testing.assert_allclose(out['number'],[3120,3120])
        self.assertEqual(out['charge'].sum(),0)
        np.testing.assert_allclose(out['thermal'],[46800,46800])
        np.testing.assert_allclose(out['momentum'][0],[6240]*3)
        floor=b.clipping_budget(n,vol,3)
        np.testing.assert_equal(floor['cells_clipped'],[3,3]);np.testing.assert_equal(floor['number_added'],[1560,1560])
    def test_open_control_volume(self):
        # 100 initially; source 20, outward transport 7, floor adds 3 => 116.
        self.assertTrue(b.balance([100],[116],[20],[7],[3],atol=1e-12,rtol=1e-12)['passed'])
        self.assertFalse(b.balance([100],[116],[20],[7],[0],atol=1e-12,rtol=1e-12)['passed'])
        self.assertFalse(b.balance([100],[116],[20],[-7],[3],atol=1e-12,rtol=1e-12)['passed'])
    def test_rejects_missing_terms_and_bad_metrics(self):
        with self.assertRaises(ValueError):b.balance([1],[1],[0],[np.nan],[0],atol=0,rtol=0)
        with self.assertRaises(ValueError):b.cell_volumes(np.ones((2,1,1)),np.ones((2,1,1)),np.ones((2,1,1)),[1,-1],[1],[1],np.ones((2,1,1),dtype=bool))
if __name__=='__main__':unittest.main()
