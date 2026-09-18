import sys
from datetime import datetime,timezone
from pathlib import Path
import tempfile
import unittest
import h5py
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from state_exchange import export,read_exchange,restore_output,cvtwin_mhd6,FIELDS,native_time
from conservative_exchange import remap,volumes,reservoir_step


class Exchange(unittest.TestCase):
    def test_conservative_remap_nonuniform_and_reject_holes(self):
        source=[np.array([0,.2,.5,1]),np.array([0,.1,1]),np.array([0,.7,1])]
        target=[np.array([0,.1,.4,.8,1]),np.array([0,.4,.8,1]),np.array([0,.2,.6,1])]
        q=np.random.default_rng(912).normal(size=(3,2,2,4))
        mapped=remap(q,source,target)
        np.testing.assert_allclose(np.sum(q*volumes(source)[...,None],axis=(0,1,2)),
                                   np.sum(mapped*volumes(target)[...,None],axis=(0,1,2)),rtol=1e-13,atol=1e-14)
        np.testing.assert_allclose(remap(np.ones_like(q)*7,source,target),7,rtol=1e-14)
        np.testing.assert_allclose(remap(q,source,source),q,rtol=2e-15,atol=1e-15)
        with self.assertRaisesRegex(ValueError,'Masked'):remap(q,source,target,valid=np.zeros((3,2,2),bool))
        target[0][-1]=1.1
        with self.assertRaisesRegex(ValueError,'domains'):remap(q,source,target)

    def test_stiff_passive_exchange_conserves_and_dissipates(self):
        for dt in [0,1e-12,.01,1,1e8]:
            a,b,t=reservoir_step(300,100,2,3,7,dt)
            self.assertAlmostEqual(2*a+3*b,900,places=11)
            self.assertLessEqual(abs(a-b),200)
            self.assertAlmostEqual(a-b,200*np.exp(-7*(1/2+1/3)*dt),places=11)
            self.assertGreaterEqual(t,0)
        with self.assertRaises(ValueError):reservoir_step(1,2,-1,1,1,1)
        with self.assertRaises(ValueError):cvtwin_mhd6(None)

    def fixture(self,p):
        shape=(2,3,4);grid=p/'grid.h5';frame=p/'20130220_18060.000000.h5'
        with h5py.File(grid,'w') as f:
            f['nullpts']=np.zeros(shape,dtype='i4');f['nullpts'][0,0,0]=1
            for i in range(3):
                f[f'e{i+1}']=np.broadcast_to(np.eye(3,dtype='f4')[i][:,None,None,None],(3,*shape))
            for key in ['r','theta','phi','glat','glon','alt','Bmag']:f[key]=np.ones(shape,dtype='f4')
            for key in ['h1','h2','h3']:f[key]=np.ones(tuple(n+4 for n in shape),dtype='f4')
            for key,n in zip(['x1','x2','x3'],shape[::-1]):f[key]=np.arange(n+4,dtype='f4')
        with h5py.File(frame,'w') as f:
            f['flagoutput']=1;f['time/ymd']=[2013,2,20];f['time/UThour']=5+1/60
            for original in FIELDS:
                sizes=(7,*shape) if original in ('nsall','Tsall','vs1all') else (shape[:2] if original=='Phiall' else shape)
                f[original]=np.random.default_rng(513).normal(size=sizes).astype('f4')
        return frame,grid

    def test_roundtrip_units_mask_basis_and_timestamp(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);frame,grid=self.fixture(p);out=p/'exchange.h5'
            export(frame,grid,out);_,mask,_,_=read_exchange(out)
            self.assertFalse(mask[0,0,0]);self.assertEqual(mask.sum(),23)
            restored=restore_output(out,p/'roundtrip.h5')
            with h5py.File(frame) as a,h5py.File(restored) as b:
                for key in FIELDS:
                    self.assertEqual(a[key].dtype,b[key].dtype)
                    np.testing.assert_array_equal(a[key][...],b[key][...])
            with h5py.File(out,'r+') as f:f['state/number_density'].attrs['units']='cm^-3'
            with self.assertRaisesRegex(ValueError,'Unit'):read_exchange(out)
            with h5py.File(frame,'r+') as f:f['time/UThour'][...]=6
            with self.assertRaisesRegex(ValueError,'UTC'):export(frame,grid,p/'bad.h5')

    def test_rejects_corrupt_geometry(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);frame,grid=self.fixture(p)
            with h5py.File(grid,'r+') as f:f['e1'][...]=0
            with self.assertRaisesRegex(ValueError,'orthonormal'):export(frame,grid,p/'bad.h5')

    def test_native_filename_rounding_and_calendar_rollover(self):
        # The native initial frame is written one microsecond after its nominal
        # filename. Check that exact observed value and all calendar boundaries.
        cases=[([2013,2,20],18000.000001,'2013-02-20T05:00:00+00:00'),
               ([2013,2,20],.055,'2013-02-20T00:00:00.060000+00:00'),
               ([2013,2,20],18060.004,'2013-02-20T05:01:00+00:00'),
               ([2013,2,20],18060.006,'2013-02-20T05:01:00.010000+00:00'),
               ([2013,2,28],86399.999,'2013-03-01T00:00:00+00:00'),
               ([2012,2,28],86399.999,'2012-02-29T00:00:00+00:00'),
               ([2013,12,31],86399.999,'2014-01-01T00:00:00+00:00')]
        for date,seconds,stamp in cases:
            hour=np.asarray(seconds/3600,dtype='f8')
            metadata=native_time(np.asarray(date),hour,datetime.fromisoformat(stamp))
            self.assertEqual(float.fromhex(metadata['native_ut_hour_hex']),float(hour))
            self.assertEqual(metadata['filename_utc'],stamp)
        instant=datetime(2013,2,20,5,1,tzinfo=timezone.utc)
        with self.assertRaisesRegex(ValueError,'mismatch'):
            native_time(np.array([2013,2,20]),np.asarray(18060.006/3600),instant)
        with self.assertRaisesRegex(ValueError,'mismatch'):
            native_time(np.array([2013,2,20]),np.asarray(18060.001/3600),
                        instant.replace(microsecond=1000))
        for date,hour in [(np.array([2013,2,29]),5.),(np.array([2013.,2.,20.]),5.),
                          (np.array([2013,2,20]),24.),(np.array([2013,2,20]),float('nan'))]:
            with self.assertRaisesRegex(ValueError,'UTC'):native_time(date,np.asarray(hour),instant)


if __name__=='__main__':unittest.main()
