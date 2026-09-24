import copy
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
p = Path(__file__).resolve().parents[2] / 'scripts/qualification/evaluate_gates.py'
s = importlib.util.spec_from_file_location('gates', p)
g = importlib.util.module_from_spec(s); s.loader.exec_module(g)


class Gates(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root/'proof.txt').write_text('Synthetic test evidence, not project qualification.\n')
        item = dict(path='proof.txt', sha256=hashlib.sha256((self.root/'proof.txt').read_bytes()).hexdigest())
        self.register = dict(candidate_commit='a'*40, gates=[dict(id=f'R{i:02d}', status='closed',
            tested_commit='a'*40, evidence=[copy.deepcopy(item)]) for i in range(1,19)])

    def test_hash_and_commit(self):
        self.assertTrue(g.evaluate(self.register,self.root)['qualified'])
        self.register['gates'][0]['tested_commit']='b'*40
        self.assertFalse(g.evaluate(self.register,self.root)['qualified'])
        self.register['gates'][0]['tested_commit']='a'*40
        (self.root/'proof.txt').write_text('changed')
        self.assertFalse(g.evaluate(self.register,self.root)['qualified'])

    def test_omission_and_escape(self):
        self.register['gates'].pop()
        with self.assertRaisesRegex(ValueError,'18 gates'):g.evaluate(self.register,self.root)
        self.register['gates'][0]['evidence'][0]['path']='../outside'
        with self.assertRaisesRegex(ValueError,'inside'):g.evaluate(self.register,self.root)
        self.register['candidate_commit']=''
        with self.assertRaisesRegex(ValueError,'commit'):g.evaluate(self.register,self.root)

    def test_approval_and_exclusion(self):
        gate=self.register['gates'][0];gate['requires_independent_approval']=True
        gate['approval_record']=True
        self.assertFalse(g.evaluate(self.register,self.root)['qualified'])
        gate['applicable']=False;gate['scope_exclusion_approved']=True
        self.assertFalse(g.evaluate(self.register,self.root)['qualified'])
        gate['scope_exclusion_record']=gate['evidence'][0]
        self.assertTrue(g.evaluate(self.register,self.root)['qualified'])


if __name__=='__main__':unittest.main()
