import copy
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from campaign import validate_manifest,admit_corpus,matched_benchmark,conformal,uq_decision,coverage
from observation_operator import radar_los,linear_analysis


class Campaign(unittest.TestCase):
    def manifest(self,p):
        events=[]
        for i,split in enumerate(['train','validation','calibration','test']):
            record=dict(id=str(i),split=split,solver_commit='a'*40,physical_event=str(i),lineage=str(i))
            for key in ['artifact','inputs','grid','forcing','initial_state']:
                path=p/f'{i}-{key}';path.write_text(str(i)+key)
                digest=hashlib.sha256(path.read_bytes()).hexdigest()
                record[key]=dict(path=path.name,sha256=digest)
            record.update(data_sha256=record['artifact']['sha256'],forcing_sha256=record['forcing']['sha256'],
                          initial_state_sha256=record['initial_state']['sha256'])
            events.append(record)
        return dict(schema='gemini.corpus.1',candidate_commit='a'*40,target_population='synthetic fixture only',
                    events=events,normalization_fit_events=['0'])

    def test_lineage_content_leakage_normalization_and_hashes(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);m=self.manifest(p);self.assertTrue(validate_manifest(m,p)['passed'])
            for key in ['physical_event','lineage','forcing_sha256','initial_state_sha256','data_sha256']:
                bad=copy.deepcopy(m);bad['events'][1][key]=bad['events'][0][key]
                with self.assertRaisesRegex(ValueError,'leakage'):validate_manifest(bad,p)
            bad=copy.deepcopy(m);bad['normalization_fit_events'].append('3')
            with self.assertRaisesRegex(ValueError,'Normalization'):validate_manifest(bad,p)
            (p/m['events'][0]['artifact']['path']).write_text('corruption')
            with self.assertRaisesRegex(ValueError,'hash'):validate_manifest(m,p)

    def test_real_gate_prerequisites_block_training(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);m=self.manifest(p)
            register=dict(candidate_commit='a'*40,gates=[dict(id=f'R{i:02d}',status='open') for i in range(1,19)])
            with self.assertRaisesRegex(ValueError,'R02,R04,R07'):admit_corpus(m,p,register,p)

    def test_matched_rollouts_baseline_and_missing_frames(self):
        truth=np.arange(12).reshape(4,3);prediction=truth+.1
        candidate=dict(prediction=prediction,model_sha256='a'*64,data_sha256='b'*64,end_to_end_seconds=.01)
        result=matched_benchmark(truth,[0,1,2],{'fixture':candidate},np.ones(3),
                                 dict(normalized_rms_max=.2,wall_seconds_max=1))
        self.assertTrue(result['passed']);self.assertLess(result['candidates']['fixture']['persistence_rms_ratio'],1)
        candidate['prediction']=prediction[:-1]
        with self.assertRaisesRegex(ValueError,'shapes'):matched_benchmark(truth,[0,1,2],{'x':candidate},np.ones(3),{})

    def test_event_conformal_insufficiency_width_shift_and_test_leakage(self):
        c=conformal([.1]*9,list(range(9)),.1)
        self.assertAlmostEqual(c['radius'],.1)
        self.assertEqual(conformal([.1]*8,list(range(8)),.1)['radius'],float('inf'))
        with self.assertRaises(ValueError):conformal([1,2],['same','same'],.1)
        self.assertEqual(uq_decision(c,[.5],[0],[1],1,'auroral','auroral')['status'],'eligible_for_research')
        r=uq_decision(c,[2],[0],[1],.1,'auroral','quiet')
        self.assertEqual(len(r['reasons']),3);self.assertFalse(r['control_enabled'])
        self.assertEqual(coverage([.05,.2],[20,21],c,list(range(9)))['empirical_coverage'],.5)
        with self.assertRaisesRegex(ValueError,'leakage'):coverage([.1],[0],c,[0])

    def test_los_sign_mask_and_linear_gaussian_analysis(self):
        v=np.array([[10,20,30],[-5,0,0]])
        self.assertAlmostEqual(radar_los(v,[[1,0,0],[1,0,0]],[.5,.5],[True,True]),2.5)
        self.assertAlmostEqual(radar_los(v,[[-1,0,0],[-1,0,0]],[.5,.5],[True,True]),-2.5)
        with self.assertRaisesRegex(ValueError,'invalid'):radar_los(v,[[1,0,0],[1,0,0]],[.5,.5],[True,False])
        a=linear_analysis([0],[[4]],[[1]],[3],[[1]])
        np.testing.assert_allclose(a['posterior'],[2.4]);np.testing.assert_allclose(a['covariance'],[[.8]])
        self.assertAlmostEqual(a['normalized_innovation_squared'],1.8)
        with self.assertRaisesRegex(ValueError,'Positive'):linear_analysis([0],[[4]],[[1]],[3],[[-1]])


if __name__=='__main__':unittest.main()
