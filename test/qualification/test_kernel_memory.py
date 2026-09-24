from pathlib import Path
import os
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from kernel_memory import measure


class KernelMemory(unittest.TestCase):
    @unittest.skipUnless(os.name=='posix','cgroup preexec hook requires POSIX')
    def test_child_pid_is_registered_before_exec(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);group=p/'gemini-qualification-test';group.mkdir()
            for name,value in {'cgroup.procs':'','cgroup.events':'populated 0\n','memory.peak':'0',
                               'memory.events':'oom_kill 0\n','memory.max':'max'}.items():
                (group/name).write_text(value)
            child_pid=p/'child.pid'
            command=[sys.executable,'-c',
                     'import os;from pathlib import Path;'
                     'Path('+repr(str(child_pid))+').write_text(str(os.getpid()));'
                     'Path('+repr(str(group/'memory.peak'))+').write_text("1024")']
            with patch('kernel_memory.uuid.uuid4') as identity, \
                 patch.object(Path,'mkdir'),patch.object(Path,'rmdir'):
                identity.return_value.hex='test'
                result=measure(command,p,p/'result.json',600,2048,660)
            self.assertTrue(result['passed'],result)
            registered_pid=int((group/'cgroup.procs').read_text())
            self.assertEqual(registered_pid,int(child_pid.read_text()))
            self.assertNotEqual(registered_pid,os.getpid())

    def test_missing_controller_blocks_before_launch(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);marker=p/'must-not-exist'
            command=[sys.executable,'-c','from pathlib import Path;Path('+repr(str(marker))+').touch()']
            result=measure(command,p,p/'result.json',600,2048,660)
            self.assertEqual(result['status'],'blocked');self.assertFalse(result['passed'])
            self.assertFalse(marker.exists());self.assertEqual(list(p.glob('gemini-qualification-*')),[])
            with self.assertRaises(ValueError):measure(command,p,p/'bad.json',float('nan'),2048,660)


if __name__=='__main__':unittest.main()
