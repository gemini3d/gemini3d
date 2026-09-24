"""Calibrated linear observation building blocks; no raw radar fitting or assimilation claim."""
import numpy as np


def radar_los(velocity_ecef,look_ecef,weights,valid):
    v,l,w=[np.asarray(x,float) for x in (velocity_ecef,look_ecef,weights)]
    valid=np.asarray(valid)
    if valid.dtype.kind not in 'biu' or not np.isin(valid,[0,1]).all():
        raise ValueError('Explicit Boolean or zero/one validity mask required')
    valid=valid.astype(bool)
    if v.ndim!=2 or v.shape[1]!=3 or l.shape!=v.shape or w.shape!=v.shape[:1] or valid.shape!=w.shape:
        raise ValueError('Matched ECEF cell vectors, quadrature weights and mask required')
    if not all(np.isfinite(x).all() for x in (v,l,w)) or np.any(w<0) or not np.isclose(w.sum(),1,atol=1e-12,rtol=0):
        raise ValueError('Finite normalized nonnegative instrument weights required')
    if np.any((~valid)&(w!=0)):raise ValueError('Instrument support intersects invalid cells')
    if not np.allclose(np.linalg.norm(l,axis=1),1,atol=1e-12,rtol=0):raise ValueError('Unit line-of-sight vectors required')
    return float(np.sum(w*np.einsum('ij,ij->i',v,l)))


def linear_analysis(prior,covariance,operator,observed,noise_covariance):
    x,P,H,y,R=[np.asarray(a,float) for a in (prior,covariance,operator,observed,noise_covariance)]
    if x.ndim!=1 or y.ndim!=1 or not x.size or not y.size or P.shape!=(x.size,x.size) or H.shape!=(y.size,x.size) or R.shape!=(y.size,y.size):
        raise ValueError('Observation/analysis shape mismatch')
    if not all(np.isfinite(a).all() for a in (x,P,H,y,R)):raise ValueError('Nonfinite analysis data')
    for cov in (P,R):
        if not np.allclose(cov,cov.T,rtol=1e-12,atol=0):raise ValueError('Covariance must be symmetric')
        try:np.linalg.cholesky(cov)
        except np.linalg.LinAlgError as e:raise ValueError('Positive definite covariance required') from e
    innovation=y-H@x;S=H@P@H.T+R
    K=np.linalg.solve(S,H@P).T;A=np.eye(x.size)-K@H
    posterior=x+K@innovation;post_cov=A@P@A.T+K@R@K.T
    return dict(posterior=posterior,covariance=(post_cov+post_cov.T)/2,
                innovation=innovation,innovation_covariance=S,
                normalized_innovation_squared=float(innovation@np.linalg.solve(S,innovation)),
                observation_dimension=y.size)
