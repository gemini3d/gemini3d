from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/qualification'))
from kernel_memory import measure


class KernelMemory(unittest.TestCase):
    def test_missing_controller_blocks_before_launch(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);marker=p/'must-not-exist'
            command=[sys.executable,'-c','from pathlib import Path;Path('+repr(str(marker))+').touch()']
            result=measure(command,p,p/'result.json',600,2048,660)
            self.assertEqual(result['status'],'blocked');self.assertFalse(result['passed'])
            self.assertFalse(marker.exists());self.assertEqual(list(p.glob('gemini-qualification-*')),[])
            with self.assertRaises(ValueError):measure(command,p,p/'bad.json',float('nan'),2048,660)


if __name__=='__main__':unittest.main()
