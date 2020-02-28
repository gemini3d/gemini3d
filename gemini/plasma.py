"""
plasma functions
"""
import typing as T


def equilibrium_resample(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    raise NotImplementedError


def equilibrium_state(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    raise NotImplementedError


def Efield_BCs2d(p: T.Dict[str, T.Any]):
    raise NotImplementedError


def Efield_BCs3d(p: T.Dict[str, T.Any]):
    raise NotImplementedError


def particles_BCs(p: T.Dict[str, T.Any]):
    raise NotImplementedError
