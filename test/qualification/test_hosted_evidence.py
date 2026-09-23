import copy
from pathlib import Path
import json,re,shlex,shutil,subprocess,sys,unittest
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
    def test_debug_requires_native_numerical_budgets(self):
        run,jobs=self.fixture()
        job=next(j for j in jobs if j['name']=='native (Debug)')
        step='Qualify exact restart and native numerical budgets'
        self.assertIn(step,REQUIRED_STEPS[job['name']])
        for conclusion in ['skipped',None]:
            with self.subTest(conclusion=conclusion):
                bad=copy.deepcopy(jobs)
                debug=next(j for j in bad if j['name']=='native (Debug)')
                if conclusion is None:
                    debug['steps']=[s for s in debug['steps'] if s['name']!=step]
                else:
                    next(s for s in debug['steps'] if s['name']==step)['conclusion']=conclusion
                self.assertFalse(verify(run,bad,'example/research','a'*40,2)['passed'])
    def test_workflow_pipelines_fail_closed(self):
        workflow=(Path(__file__).resolve().parents[2]/WORKFLOW).read_text()
        self.assertRegex(workflow,r'(?m)^defaults:\n  run:\n    shell: bash\s*$')
    @unittest.skipUnless(shutil.which('bash'),'Bash is required to exercise GitHub shell semantics')
    def test_explicit_github_bash_rejects_pipeline_failure(self):
        # GitHub adds pipefail only when the Bash shell is explicitly selected.
        result=subprocess.run(['bash','--noprofile','--norc','-e','-o','pipefail','-c',
                               'false | tee'],capture_output=True,text=True,timeout=10)
        self.assertNotEqual(result.returncode,0)
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
        research=re.search(r'(?ms)^      - name: Qualify exact restart and native numerical budgets\n'
                           r'        run: \|\n(.*?)(?=^      - |\Z)',jobs['native'])
        self.assertIsNotNone(research)
        command=shlex.split(research[1].replace('\\\n',' '),comments=True)
        command=command[:command.index('|')]
        self.assertEqual(command[:2],['python','test/qualification/research_matrix.py'])
        self.assertIn('--extended-layouts',command)
        self.assertNotIn('--quick',command)
        self.assertEqual(command[command.index('--work')+1],'build/research')
        for case in ['mini2dns_fang','mini2dew_fang','mini3d_fang','mini2dns_fang_cpp',
                     'mini2dns_glow','mini2dew_glow','mini3d_glow','mini2dns_glow_cpp','mini2dns_msis2_fang']:
            self.assertIn('"'+case+'"',jobs['native'])
        self.assertNotIn('--quick',workflow)
        self.assertIn('host_capabilities.py --only leak',jobs['sanitizers'])
        self.assertIn('detect_leaks=1:halt_on_error=1',jobs['sanitizers'])
        self.assertIn('LSAN_OPTIONS: exitcode=23:suppressions=${{ github.workspace }}/.github/lsan-openmpi.supp',
                      jobs['sanitizers'])
        self.assertEqual((root/'.github/lsan-openmpi.supp').read_text().strip().splitlines(),
                         ['# OpenMPI 4.x on Ubuntu 24.04 leaks hwloc allocations during MPI_Init().',
                          'leak:ompi_mpi_init'])
        self.assertIn('kernel_memory.py',jobs['current-hdf5'])
        self.assertIn('for repetition in 1 2 3',jobs['current-hdf5'])
        presets=json.loads((root/'CMakePresets.json').read_text())['configurePresets']
        qualification=next(p for p in presets if p['name']=='qualification')['cacheVariables']
        for key in ['BUILD_TESTING','gemini3d_BUILD_TESTING','gemini3d_require_qualification',
                    'gemini3d_test_simulations','gemini3d_glow','gemini3d_msis2']:
            self.assertIs(qualification[key],True)
        self.assertIs(qualification['gemini3d_hwm14'],False)
        self.assertIn("mini2dns_fang_cpp",jobs['sanitizers'])
    def test_install_smoke_exercises_real_entrypoints(self):
        root=Path(__file__).resolve().parents[2]
        workflow=(root/'.github/workflows/install-smoke.yml').read_text()
        self.assertRegex(workflow,r'(?m)^  pull_request:\s*$')
        self.assertRegex(workflow,r'(?m)^  push:\s*$')
        self.assertNotIn('paths-ignore:',workflow)
        self.assertNotIn('continue-on-error:',workflow)
        for platform in ['ubuntu-24.04','macos-14','Ubuntu-24.04']:
            self.assertIn(platform,workflow)
        jobs=dict(re.findall(r'(?ms)^  (unix|wsl):\n(.*?)(?=^  \w[\w-]*:\n|\Z)',workflow))
        self.assertEqual(set(jobs),{'unix','wsl'})
        self.assertRegex(jobs['wsl'],r'(?m)^    runs-on: windows-(?:latest|\d{4})\s*$')
        bash='bash scripts/install-local.sh --system-deps --root "$PWD/build/local-install" --jobs 2'
        self.assertEqual(workflow.count(bash),1)
        self.assertEqual(jobs['unix'].count(bash),1)
        self.assertIn('& ./scripts/install-local.ps1 -Distribution Ubuntu-24.04 -SystemDeps',jobs['wsl'])
        self.assertIn('-Root "$($linuxRoot.Trim())/build/local-install" -Jobs 2',jobs['wsl'])
        self.assertRegex(jobs['wsl'],r'(?m)^      - name: Install through the WSL local entrypoint\n        shell: pwsh\s*$')
        self.assertEqual(workflow.count('scripts/local_environment.py check'),2)
        self.assertEqual(workflow.count('build/local-install/environment.json'),2)
        self.assertEqual(workflow.count('if-no-files-found: error'),2)
if __name__=='__main__':unittest.main()
