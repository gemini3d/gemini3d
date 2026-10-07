"""Volume-weighted physical-cell budgets with explicit sources and boundaries.

Inputs use solver order (x1,x2,x3,species), SI units, and a caller-supplied
physical mask. Singleton dimensions need an explicit physical thickness.
Output differences alone do not establish conservation: the independent
source/boundary/clipping integrals must be supplied to balance().
"""
import numpy as np


def cell_volumes(h1, h2, h3, dx1, dx2, dx3, physical):
    physical=np.asarray(physical)
    if physical.dtype!=np.bool_ or physical.ndim!=3:
        raise ValueError('physical must be an explicit three-dimensional boolean mask')
    shape=physical.shape
    metrics=[np.asarray(x,dtype=float) for x in (h1,h2,h3)]
    widths=[np.asarray(x,dtype=float) for x in (dx1,dx2,dx3)]
    if any(m.shape!=shape for m in metrics) or any(w.shape!=(shape[i],) for i,w in enumerate(widths)):
        raise ValueError('Metric, mask or cell-width shape mismatch')
    if any(not np.isfinite(x).all() or np.any(x<=0) for x in metrics+widths):
        raise ValueError('Metrics and physical widths must be finite and positive')
    v=metrics[0]*metrics[1]*metrics[2]*widths[0][:,None,None]*widths[1][None,:,None]*widths[2][None,None,:]
    if not np.isfinite(v).all():raise ValueError('Volume overflow')
    return np.where(physical,v,0.0)


def integrals(n, velocity, temperature, volume, masses, charges, gammas, kb):
    n=np.asarray(n,dtype=float);v=np.asarray(velocity,dtype=float);t=np.asarray(temperature,dtype=float)
    vol=np.asarray(volume,dtype=float);m=np.asarray(masses);q=np.asarray(charges);g=np.asarray(gammas)
    if n.ndim!=4 or n.shape[:3]!=vol.shape or t.shape!=n.shape or v.shape!=n.shape+(3,):
        raise ValueError('State/mask shape mismatch; ghost cells must be removed')
    if any(a.shape!=(n.shape[-1],) for a in (m,q,g)):raise ValueError('Species metadata mismatch')
    if not all(np.isfinite(a).all() for a in (vol,m,q,g)) or np.any(vol<0) or np.any(m<=0) or np.any(g<=1):
        raise ValueError('Invalid metrics or species constants')
    if not np.isfinite(kb) or kb<=0:raise ValueError('Invalid Boltzmann constant')
    physical=vol>0
    if not all(np.isfinite(a[physical]).all() for a in (n,v,t)):
        raise ValueError('Nonfinite state in physical cells')
    if np.any(n[physical]<0) or np.any(t[physical]<=0):raise ValueError('Invalid physical state')
    # Index first: multiplying a masked NaN by zero would still propagate NaN.
    weights=vol[physical,None];ns=n[physical];vs=v[physical];ts=t[physical]
    number=np.sum(weights*ns,axis=0)
    momentum=np.sum(weights[:,:,None]*ns[:,:,None]*m[None,:,None]*vs,axis=0)
    thermal=np.sum(weights*ns*kb*ts/(g-1),axis=0)
    kinetic=np.sum(weights*ns*m*np.sum(vs**2,axis=-1)/2,axis=0)
    result=dict(number=number,mass=number*m,charge=number*q,momentum=momentum,thermal=thermal,kinetic=kinetic)
    if not all(np.isfinite(a).all() for a in result.values()):raise ValueError('Integral overflow')
    return result


def clipping_budget(n, volume, floor):
    n=np.asarray(n,dtype=float);volume=np.asarray(volume,dtype=float)
    if n.ndim!=4 or n.shape[:3]!=volume.shape or not np.isfinite(floor) or floor<=0:
        raise ValueError('Invalid clipping budget inputs')
    if not np.isfinite(volume).all() or np.any(volume<0):raise ValueError('Invalid volume')
    physical=volume>0
    if not np.isfinite(n[physical]).all():raise ValueError('Nonfinite physical density')
    delta=np.maximum(floor-n[physical],0)
    return dict(number_added=np.sum(delta*volume[physical,None],axis=0),cells_clipped=np.sum(delta>0,axis=0))


def balance(before, after, source_integral, outward_boundary_integral, clipping_integral, *, atol, rtol):
    # All terms are integrals over the same interval and physical control volume.
    arrays=[np.asarray(a,dtype=float) for a in (before,after,source_integral,outward_boundary_integral,clipping_integral)]
    if any(a.shape!=arrays[0].shape or not np.isfinite(a).all() for a in arrays):
        raise ValueError('All five budget terms must have equal shape and finite values')
    if not np.isfinite([atol,rtol]).all() or min(atol,rtol)<0:raise ValueError('Invalid acceptance bounds')
    b,a,s,f,c=arrays
    residual=a-b-s+f-c
    scale=np.maximum.reduce([np.abs(x) for x in arrays])
    bound=atol+rtol*scale
    return dict(residual=residual,bound=bound,passed=bool(np.all(np.abs(residual)<=bound)))
