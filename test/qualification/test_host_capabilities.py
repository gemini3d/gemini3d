from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
from types import SimpleNamespace
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from host_capabilities import leak_probe

class LeakControls(unittest.TestCase):
    def run_probe(self,leak_output,clean_code=0):
        with tempfile.TemporaryDirectory() as d:
            work=Path(d)/'probe'
            def fake(command,**kwargs):
                if '-fsanitize=address' in command:
                    (work/'leak_control').write_bytes(b'fixture')
                    return SimpleNamespace(returncode=0,stdout='',stderr='')
                is_leak=command[-1]=='leak'
                return SimpleNamespace(returncode=23 if is_leak else clean_code,stdout='',stderr=leak_output if is_leak else '')
            with patch('host_capabilities.subprocess.run',side_effect=fake):return leak_probe('cc',work)
    def test_detected_positive_control(self):
        self.assertTrue(self.run_probe('LeakSanitizer: detected memory leaks\n1237 byte(s) leaked')['passed'])
    def test_runtime_crash_is_not_leak_detection(self):
        self.assertFalse(self.run_probe('LeakSanitizer has encountered a fatal error')['passed'])
        self.assertFalse(self.run_probe('LeakSanitizer: detected memory leaks\n1237 byte(s) leaked',23)['passed'])
if __name__=='__main__':unittest.main()
