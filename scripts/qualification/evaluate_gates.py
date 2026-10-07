"""Check qualification evidence integrity, failing closed on open or stale gates.

This verifies records and hashes, not the identity or authority of a reviewer.
The release owner must independently authenticate every approval record.
Exit 2 means the candidate remains unqualified; no status is changed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re


def evaluate(register, root):
    root = root.resolve()
    commit = register.get('candidate_commit', '')
    if not re.fullmatch(r'[0-9a-f]{40}', commit):
        raise ValueError('An exact candidate commit is required')
    seen, issues = set(), []
    expected = {f'R{i:02d}' for i in range(1, 19)}

    def record(gid, item, kind):
        if not isinstance(item, dict) or not {'path', 'sha256'} <= item.keys():
            issues.append(dict(id=gid, reason=kind + ' must name a file and SHA256'))
            return
        rel = Path(item['path'])
        path = (root / rel).resolve()
        if rel.is_absolute() or not path.is_relative_to(root):
            raise ValueError('Evidence must remain inside package')
        if not path.is_file():
            issues.append(dict(id=gid, reason=kind + ' missing: ' + str(rel)))
        elif hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            issues.append(dict(id=gid, reason=kind + ' hash changed: ' + str(rel)))

    for gate in register['gates']:
        gid = gate['id']
        if gid in seen or gid not in expected:
            raise ValueError('Duplicate or unknown gate ' + gid)
        seen.add(gid)
        if gate.get('applicable', True) is False:
            if gate.get('scope_exclusion_approved') is not True:
                issues.append(dict(id=gid, reason='Scope exclusion needs approval'))
            record(gid, gate.get('scope_exclusion_record'), 'Scope exclusion')
            continue
        if gate.get('status') != 'closed':
            issues.append(dict(id=gid, reason=gate.get('required_change', 'Gate remains open')))
            continue
        if not gate.get('evidence'):
            issues.append(dict(id=gid, reason='No evidence attached'))
        for evidence in gate.get('evidence', []):
            record(gid, evidence, 'Evidence')
        if gate.get('requires_independent_approval', False):
            record(gid, gate.get('approval_record'), 'Independent approval')
        if gate.get('tested_commit') != commit:
            issues.append(dict(id=gid, reason='Evidence does not name the candidate commit'))
    if seen != expected:
        raise ValueError('Register must retain all 18 gates')
    return dict(qualified=not issues, open_gates=issues, candidate_commit=commit,
                scope='Record integrity only; approval authority is verified by the release owner')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('register', type=Path)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    try:
        result = evaluate(json.loads(a.register.read_text()), a.root)
    except (ValueError, KeyError, TypeError, OSError) as e:
        result = dict(qualified=False, error=str(e))
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return 0 if result['qualified'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
