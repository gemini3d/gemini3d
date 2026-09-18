import datetime as dt
import importlib.util
from pathlib import Path
import tempfile
import unittest
import sys
import h5py
import numpy as np
p=Path(__file__).resolve().parents[2]/'scripts/qualification/validate_driver.py'
sys.path.insert(0,str(p.parent))
s=importlib.util.spec_from_file_location('driver',p);d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
class Driver(unittest.TestCase):
    def test_midnight_gaps_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as path:
            root=Path(path);names=['20261231_86340.000000.h5','20270101_00000.000000.h5','20270101_00060.000000.h5']
            for name in names:
                with h5py.File(root/name,'w') as f:f['flux']=np.zeros((2,3))
            start=d.utc('2026-12-31T23:59:00');stop=d.utc('2027-01-01T00:01:00')
            self.assertTrue(d.validate(root,start,stop,60)['passed'])
            with h5py.File(root/names[1],'r+') as f:f['flux'][0,0]=np.nan
            with self.assertRaisesRegex(ValueError,'Nonfinite'):d.validate(root,start,stop,60)
            (root/names[1]).unlink()
            with self.assertRaisesRegex(ValueError,'cadence'):d.validate(root,start,stop,60)
            with self.assertRaisesRegex(ValueError,'cover'):d.validate(root,start,stop+dt.timedelta(seconds=60),120)
    def test_external_link_and_shape_change(self):
        with tempfile.TemporaryDirectory() as path:
            root=Path(path)
            with h5py.File(root/'20260101_00000.000000.h5','w') as f:f['flux']=[1,2]
            with h5py.File(root/'20260101_00060.000000.h5','w') as f:f['flux']=[1,2,3]
            args=(root,d.utc('2026-01-01'),d.utc('2026-01-01T00:01:00'),60)
            with self.assertRaisesRegex(ValueError,'shape'):d.validate(*args)
            with h5py.File(root/'20260101_00060.000000.h5','w') as f:f['flux']=h5py.ExternalLink('missing.h5','/flux')
            with self.assertRaisesRegex(ValueError,'links'):d.validate(*args)
    def test_timestamp_mismatch(self):
        with tempfile.TemporaryDirectory() as path:
            root=Path(path)
            for sec in [0,60]:
                with h5py.File(root/f'20260101_{sec:05d}.000000.h5','w') as f:
                    f['flux']=[1.,2.];f['time/ymd']=[2026,1,1];f['time/UTsec']=sec
            args=(root,d.utc('2026-01-01'),d.utc('2026-01-01T00:01:00'),60)
            self.assertTrue(d.validate(*args)['passed'])
            with h5py.File(root/'20260101_00060.000000.h5','r+') as f:f['time/UTsec'][...]=0
            with self.assertRaisesRegex(ValueError,'metadata'):d.validate(*args)
class RequiredSchema(unittest.TestCase):
    def fixture(self,root,kind):
        with h5py.File(root/'simsize.h5','w') as f:f['Nlat']=2;f['Nlon']=3
        with h5py.File(root/'simgrid.h5','w') as f:f['mlat']=[65.,66.];f['mlon']=[200.,201.,202.]
        for sec in [0,60]:
            with h5py.File(root/f'20260101_{sec:05d}.000000.h5','w') as f:
                f['time/ymd']=[2026,1,1];f['time/UTsec']=sec
                if kind=='precip':f['Qp']=np.ones((2,3));f['E0p']=np.ones((2,3))*1000
                else:
                    for k in ['Exit','Eyit','Vminx1it','Vmaxx1it']:f[k]=np.zeros((2,3))
                    for k in ['Vminx2ist','Vmaxx2ist']:f[k]=np.zeros(2)
                    for k in ['Vminx3ist','Vmaxx3ist']:f[k]=np.zeros(3)
                    f['flagdirich']=1
        return (root,d.utc('2026-01-01'),d.utc('2026-01-01T00:01:00'),60,kind)
    def test_required_fields_axis_order_and_flags(self):
        for mutation,expected in [('missing','required field'),('transpose','shape'),('flag','categorical'),('time','metadata')]:
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as path:
                root=Path(path);args=self.fixture(root,'efield')
                self.assertTrue(d.validate(*args)['required_schema_checked'])
                with h5py.File(root/'20260101_00000.000000.h5','r+') as f:
                    if mutation=='missing':del f['Exit']
                    if mutation=='transpose':del f['Exit'];f['Exit']=np.zeros((3,2))
                    if mutation=='flag':f['flagdirich'][...]=9
                    if mutation=='time':del f['time/UTsec']
                with self.assertRaisesRegex(ValueError,expected):d.validate(*args)
    def test_precipitation_sign_and_coordinate_grid(self):
        for mutation,expected in [('flux','Negative'),('energy','positive'),('coordinate','increasing'),('link','link')]:
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as path:
                root=Path(path);args=self.fixture(root,'precip')
                self.assertTrue(d.validate(*args)['passed'])
                if mutation=='coordinate':
                    with h5py.File(root/'simgrid.h5','r+') as f:f['mlon'][...]=[202.,201.,200.]
                elif mutation=='link':
                    with h5py.File(root/'simgrid.h5','r+') as f:del f['mlat'];f['mlat']=h5py.SoftLink('/mlon')
                else:
                    with h5py.File(root/'20260101_00000.000000.h5','r+') as f:
                        f['Qp' if mutation=='flux' else 'E0p'][0,0]=-1
                with self.assertRaisesRegex(ValueError,expected):d.validate(*args)
if __name__=='__main__':unittest.main()
