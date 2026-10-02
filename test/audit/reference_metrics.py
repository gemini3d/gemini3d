"""Report every saved species, without changing pass tolerances or references.

Full-grid unweighted norms are diagnostics, not physical validation. Species
axis is the first HDF5 dimension. No null-cell/volume weighting is inferred.
"""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np

SPECIES=['O+','NO+','N2+','O2+','N+','H+','e-']
FIELDS=['nsall','Tsall','vs1all','Phiall','J1all','J2all','J3all','v2avgall','v3avgall']

def compare(new,reference):
    frames=sorted(reference.glob('????????_?????.??????.h5'))
    records=[]
    if not frames:raise ValueError('No reference frames: '+str(reference))
    for ref in frames:
        candidate=new/ref.name
        if not candidate.exists():raise ValueError('Missing candidate frame: '+str(candidate))
        with h5py.File(candidate) as n,h5py.File(ref) as r:
            for field in FIELDS:
                if field not in r:continue
                if field not in n:raise ValueError('Missing candidate field: '+field)
                a=n[field][...].astype(float);b=r[field][...].astype(float)
                if a.shape!=b.shape:raise ValueError('Shape mismatch: '+field)
                blocks=[(s,a[i],b[i]) for i,s in enumerate(SPECIES)] if field in FIELDS[:3] else [(None,a,b)]
                for sp,aa,bb in blocks:
                    delta=aa-bb
                    records.append(dict(frame=ref.name,field=field,species=sp,
                        candidate_finite=bool(np.isfinite(aa).all()),reference_finite=bool(np.isfinite(bb).all()),
                        candidate_min=float(aa.min()),reference_min=float(bb.min()),
                        rel_l2=float(np.linalg.norm(delta)/max(np.linalg.norm(bb),1e-30)),
                        max_abs=float(np.max(np.abs(delta))),bitwise_equal=bool(np.array_equal(aa,bb))))
    return records

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--new',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();records=compare(a.new,a.reference)
    payload=dict(schema='gemini.audit.species_metrics.1',new=str(a.new),reference=str(a.reference),
        scope='full-grid unweighted saved fields; no acceptance thresholds; no reference updates',records=records)
    a.output.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print('Wrote',len(records),'field/species/frame metrics')

if __name__=='__main__':main()
