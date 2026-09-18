import copy
from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from hosted_evidence import verify,REQUIRED,WORKFLOW,REQUIRED_STEPS
class HostedEvidence(unittest.TestCase):
    def fixture(self):
        run=dict(repository={'full_name':'example/research'},head_sha='a'*40,path=WORKFLOW,status='completed',conclusion='success',run_attempt=2,id=99)
        jobs=[dict(name=name,status='completed',conclusion='success',run_id=99,run_attempt=2,head_sha='a'*40,
                   steps=[dict(name=n,conclusion='success') for n in REQUIRED_STEPS[name]]) for name in REQUIRED]
        return run,jobs
    def test_exact_success_and_stale_identity(self):
        run,jobs=self.fixture();self.assertTrue(verify(run,jobs,'example/research','a'*40,2)['passed'])
        for key,value in [('head_sha','b'*40),('run_attempt',1),('run_id',98)]:
            bad=copy.deepcopy(jobs);bad[0][key]=value
            self.assertFalse(verify(run,bad,'example/research','a'*40,2)['passed'])
    def test_missing_skipped_duplicate_and_masked_failure(self):
        run,jobs=self.fixture()
        for bad in [jobs[:-1],jobs+[jobs[0]]]:self.assertFalse(verify(run,bad,'example/research','a'*40,2)['passed'])
        for verdict in ['skipped','failure','cancelled']:
            bad=copy.deepcopy(jobs);bad[0]['conclusion']=verdict
            self.assertFalse(verify(run,bad,'example/research','a'*40,2)['passed'])
        jobs[0]['steps'][0]['conclusion']='failure'
        self.assertFalse(verify(run,jobs,'example/research','a'*40,2)['passed'])
    def test_skipped_required_step_fails_successful_job(self):
        run,jobs=self.fixture();jobs[0]['steps'][0]['conclusion']='skipped'
        self.assertFalse(verify(run,jobs,'example/research','a'*40,2)['passed'])
if __name__=='__main__':unittest.main()
