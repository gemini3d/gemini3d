import copy
from pathlib import Path
import json,re,sys,unittest
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
    def test_workflow_preserves_required_executable_contract(self):
        root=Path(__file__).resolve().parents[2]
        workflow=(root/WORKFLOW).read_text()
        self.assertRegex(workflow,r'(?m)^  pull_request:\s*$')
        self.assertRegex(workflow,r'(?m)^  push:\s*$')
        self.assertNotIn('paths-ignore:',workflow)
        self.assertNotIn('continue-on-error:',workflow)
        jobs=dict(re.findall(r'(?ms)^  (native|current-hdf5|sanitizers):\n(.*?)(?=^  \w[\w-]*:\n|\Z)',workflow))
        self.assertEqual(set(jobs),{'native','current-hdf5','sanitizers'})
        for name,steps in REQUIRED_STEPS.items():
            block=jobs['native' if name.startswith('native ') else name]
            for step in steps:
                self.assertEqual(block.count('- name: '+step+'\n'),1,(name,step))
            self.assertIn('uses: ./.github/workflows/qualification-evidence',block)
            self.assertIn('ref: ${{ env.CANDIDATE_COMMIT }}',block)
        self.assertIn('build_type: Debug',jobs['native'])
        self.assertIn('build_type: Release',jobs['native'])
        self.assertIn('research_matrix.py',jobs['native'])
        for case in ['mini2dns_fang','mini2dew_fang','mini3d_fang','mini2dns_fang_cpp',
                     'mini2dns_glow','mini2dew_glow','mini3d_glow','mini2dns_glow_cpp','mini2dns_msis2_fang']:
            self.assertIn('"'+case+'"',jobs['native'])
        self.assertNotIn('--quick',workflow)
        self.assertIn('host_capabilities.py --only leak',jobs['sanitizers'])
        self.assertIn('detect_leaks=1:halt_on_error=1',jobs['sanitizers'])
        self.assertIn('kernel_memory.py',jobs['current-hdf5'])
        self.assertIn('for repetition in 1 2 3',jobs['current-hdf5'])
        presets=json.loads((root/'CMakePresets.json').read_text())['configurePresets']
        qualification=next(p for p in presets if p['name']=='qualification')['cacheVariables']
        for key in ['BUILD_TESTING','gemini3d_BUILD_TESTING','gemini3d_require_qualification',
                    'gemini3d_test_simulations','gemini3d_glow','gemini3d_msis2']:
            self.assertIs(qualification[key],True)
        self.assertIs(qualification['gemini3d_hwm14'],False)
        self.assertIn("mini2dns_fang_cpp",jobs['sanitizers'])
if __name__=='__main__':unittest.main()
