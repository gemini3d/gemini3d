"""Verified Cartesian cell-average remapping and two-reservoir exchange kernel.

This is a manufactured coupling component, not a GEMINI/CVTWIN physical law.
Curved metrics, open boundaries and masked cut cells require a different operator.
"""
import numpy as np


def edges(value):
    x=np.asarray(value,dtype=float)
    if x.ndim!=1 or len(x)<2 or not np.isfinite(x).all() or not np.all(np.diff(x)>0):
        raise ValueError('Finite increasing cell edges required')
    return x


def overlap(source,target):
    a,b=edges(source),edges(target)
    # Exact matching endpoints prevent silent uncovered-volume extrapolation.
    if a[0]!=b[0] or a[-1]!=b[-1]:raise ValueError('Source and target domains must match')
    return np.maximum(0.,np.minimum(b[1:,None],a[None,1:])-np.maximum(b[:-1,None],a[None,:-1]))


def volumes(grid):
    if len(grid)!=3:raise ValueError('Three Cartesian axes required')
    return np.einsum('i,j,k->ijk',*[np.diff(edges(x)) for x in grid])


def remap(cell_average,source,target,valid=None):
    """Array axes follow the provided Cartesian edges; trailing axis is channels."""
    value=np.asarray(cell_average,dtype=float);vs=volumes(source);vt=volumes(target)
    if value.ndim!=4 or value.shape[:3]!=vs.shape or not np.isfinite(value).all():
        raise ValueError('Expected finite cell averages with explicit channel axis')
    if valid is not None and (np.shape(valid)!=vs.shape or not np.asarray(valid,dtype=bool).all()):
        raise ValueError('Masked or cut cells are outside this exchange contract')
    matrices=[overlap(a,b) for a,b in zip(source,target)]
    total=np.einsum('ai,bj,ck,ijkq->abcq',*matrices,value,optimize=True)
    return total/vt[...,None]


def reservoir_step(a,b,capacity_a,capacity_b,conductance,dt):
    """Exact passive exchange C_a da/dt=-k(a-b), C_b db/dt=k(a-b)."""
    values=np.asarray([a,b,capacity_a,capacity_b,conductance,dt],float)
    if not np.isfinite(values).all() or min(capacity_a,capacity_b)<=0 or min(conductance,dt)<0:
        raise ValueError('Invalid passive exchange parameters')
    transfer=(a-b)*(-np.expm1(-conductance*(1/capacity_a+1/capacity_b)*dt))/(1/capacity_a+1/capacity_b)
    return a-transfer/capacity_a,b+transfer/capacity_b,transfer
