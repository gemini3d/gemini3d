import functools
import typing as T
from pathlib import Path
import numpy as np
from datetime import datetime, timedelta

NaN = float("NaN")
R = Path(__file__).resolve().parents[1]

__all__ = ["read_config"]


@functools.lru_cache()
def read_config(path: Path) -> T.Dict[str, T.Any]:
    """
    read simulation input configuration

    .nml is strongly preferred, .ini is legacy.

    Parameters
    ----------
    path: pathlib.Path
        config file path

    Returns
    -------
    params: dict
        simulation parameters from config file
    """
    P: T.Dict[str, T.Any] = {}

    path = Path(path).expanduser().resolve()
    if path.is_file():
        file = path
        path = path.parent
        if file.suffix == ".nml":
            P = read_nml(file)
        elif file.suffix == ".ini":
            P = read_ini(file)
        else:
            raise ValueError(f"unsure how to read {file}")
    elif path.is_dir():
        for suff in (".nml", ".ini"):
            file = path / ("config" + suff)
            if file.is_file():
                P = read_config(file)
                break
    else:
        raise FileNotFoundError(f"could not find config file in {path}")

    if not P:
        raise FileNotFoundError(f"could not find config file in {path}")

    P["nml"] = file

    return P


def read_nml(fn: Path) -> T.Dict[str, T.Any]:
    """ parse .nml file
    for now we don't use the f90nml package, though maybe we will in the future.
    Just trying to keep Python prereqs reduced for this simple parsing.
    """

    if fn.is_dir():
        fn = fn / "config.nml"

    groups = ["base", "setup"]
    params = read_nml_group(fn, groups)

    groups = []
    if params["flagdneu"]:
        groups.append("neutral_perturb")
    if params["flagprecfile"]:
        groups.append("precip")
    if params["flagE0file"]:
        groups.append("efield")
    if params["flagglow"]:
        groups.append("glow")

    params.update(read_nml_group(fn, groups))

    return params


def read_nml_group(fn: Path, group: T.Sequence[str]) -> T.Dict[str, T.Any]:
    """ read a group from an .nml file """

    groups = tuple(group)
    if not groups:
        return {}

    raw: T.Dict[str, T.Sequence[str]] = {}

    with fn.open("r") as f:
        for line in f:
            ls = line.strip()
            if ls.startswith('&'):
                # a group name
                if not ls[1:].startswith(groups):
                    continue
            for line in f:
                if ls.startswith("/"):
                    # end of group
                    break
                vals = line.split("=")
                if len(vals) != 2:
                    # not a valid line
                    continue

                key = vals[0].strip()
                values = [v.strip().replace("'", "").replace('"', "") for v in vals[1].split("!")[0].split(",")]
                raw[key] = values[0] if len(values) == 1 else values

    if not raw:
        raise KeyError(f"did not find Namelist group(s) {groups} in {fn}")

    return parse_group(raw, groups)


def parse_group(raw: T.Dict[str, T.Any], group: T.Sequence[str]) -> T.Dict[str, T.Any]:
    """
    this is Gemini-specific
    don't resolve absolute paths here because that assumes same machine
    """

    if isinstance(group, str):
        group = [group]

    P: T.Dict[str, T.Any] = {}

    if "base" in group:
        P["t0"] = datetime(int(raw["ymd"][0]), int(raw["ymd"][1]), int(raw["ymd"][2])) + timedelta(seconds=float(raw["UTsec0"]))
        P["tdur"] = timedelta(seconds=float(raw["tdur"]))
        P["dtout"] = timedelta(seconds=float(raw["dtout"]))
        P["f107a"] = float(raw["activ"][0])
        P["f107"] = float(raw["activ"][1])
        P["Ap"] = float(raw["activ"][2])
        P["tcfl"] = float(raw["tcfl"])
        P["Teinf"] = float(raw["Teinf"])
        for k in ("potsolve", "flagperiodic", "flagoutput", "flagcap", "flagglow", "flagE0file", "flagdneu", "flagprecfile"):
            if k in raw:
                P[k] = int(raw[k])
            else:
                P[k] = 0
        for k in ("indat_file", "indat_grid", "indat_size"):
            P[k] = Path(raw[k])

    if "setup" in group:
        P["format"] = raw["format"]
        P["alt_scale"] = np.array(list(map(float, raw["alt_scale"])))
        for k in ("realbits", "lxp", "lyp"):
            P[k] = int(raw[k])
        for k in ("glat", "glon", "xdist", "ydist", "alt_min", "alt_max", "Bincl", "nmf", "nme"):
            if k in raw:
                P[k] = float(raw[k])

    if "neutral_perturb" in group:
        P["interptype"] = int(raw["interptype"])
        P["sourcedir"] = Path(raw["source_dir"])

        for k in ("sourcemlat", "sourcemlon", "dtneu", "dxn", "drhon", "dzn"):
            try:
                P[k] = float(raw[k])
            except KeyError:
                P[k] = NaN

    if "precip" in group:
        P["dtprec"] = float(raw["dtprec"])
        P["precdir"] = Path(raw["prec_dir"])

    if "efield" in group:
        P["dtE0"] = float(raw["dtE0"])
        P["E0dir"] = Path(raw["E0_dir"])

    if "glow" in group:
        P["dtglow"] = float(raw["dtglow"])
        P["dtglowout"] = float(raw["dtglowout"])

    if not P:
        raise ValueError(f"Not sure how to parse NML group {group}")

    return P


def read_ini(fn: Path) -> T.Dict[str, T.Any]:
    """ parse .ini file (legacy) """

    P: T.Dict[str, T.Any] = {}

    with fn.open("r") as f:
        date = list(map(int, f.readline().split()[0].split(",")))[::-1]
        sec = float(f.readline().split()[0])

        P["t0"] = datetime(*date) + timedelta(seconds=sec)  # type: ignore  # mypy bug

        P["tdur"] = timedelta(seconds=float(f.readline().split()[0]))

        P["dtout"] = timedelta(seconds=float(f.readline().split()[0]))

        P["f107a"], P["f107"], P["Ap"] = map(float, f.readline().split()[0].split(","))

        P["tcfl"] = float(f.readline().split()[0])
        P["Teinf"] = float(f.readline().split()[0])

        P["potsolve"] = int(f.readline().split()[0])
        P["flagperiodic"] = int(f.readline().split()[0])
        P["flagoutput"] = int(f.readline().split()[0])
        P["flagcap"] = int(f.readline().split()[0])

        for k in ("indat_size", "indat_grid", "indat_file"):
            P[k] = Path(f.readline().strip().replace("'", "").replace('"', "")).expanduser()

    return P
