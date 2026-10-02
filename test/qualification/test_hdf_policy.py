from pathlib import Path
import sys
import tempfile
import unittest
import h5py
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from hdf_policy import numeric,inspect_dataset


class Policy(unittest.TestCase):
    def test_external_storage_vds_and_link_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            with h5py.File(p/'test.h5','w') as f:
                ext=f.create_dataset('external',(10,),dtype='f8',external=[('never-read.bin',0,80)])
                with self.assertRaisesRegex(ValueError,'External'):inspect_dataset(ext)
                layout=h5py.VirtualLayout((10,),dtype='f8');layout[:]=h5py.VirtualSource('missing.h5','x',shape=(10,))
                vds=f.create_virtual_dataset('vds',layout)
                with self.assertRaisesRegex(ValueError,'Virtual'):inspect_dataset(vds)
                f['linked']=h5py.ExternalLink('missing.h5','group')
                with self.assertRaisesRegex(ValueError,'indirect'):numeric(f,'linked/value')

    def test_limit_checked_without_allocating_and_builtin_filter(self):
        with tempfile.TemporaryDirectory() as d:
            with h5py.File(Path(d)/'test.h5','w') as f:
                sparse=f.create_dataset('large',(2**30,),dtype='f8',chunks=(1024,))
                with self.assertRaisesRegex(ValueError,'allocation'):inspect_dataset(sparse)
                f.create_dataset('compressed',data=np.arange(30),compression='gzip',shuffle=True,fletcher32=True)
                np.testing.assert_array_equal(numeric(f,'compressed'),np.arange(30))
                f.create_dataset('lossy',data=np.arange(10,dtype=float),scaleoffset=2)
                with self.assertRaisesRegex(ValueError,'filter'):numeric(f,'lossy')


if __name__=='__main__':unittest.main()
