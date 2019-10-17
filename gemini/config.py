import functools
import typing
from pathlib import Path
from datetime import datetime, timedelta

from . import raw

NaN = float("NaN")
R = Path(__file__).resolve().parents[1]


@functools.lru_cache()
def read_config(path: Path) -> typing.Dict[str, typing.Any]:
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

    # NOT P["indat_size"] because that assumes the reading computer has the same directory layout as HPC
    simsize_file = path / "simsize.dat"
    if simsize_file.is_file():
        P["lxs"] = raw.get_simsize(simsize_file)

    return P


def read_nml(fn: Path) -> typing.Dict[str, typing.Any]:
    """ parse .nml file
    for now we don't use the f90nml package, though maybe we will in the future.
    Just trying to keep Python prereqs reduced for this simple parsing.
    """

    params = {}
    params.update(read_nml_group(fn, "base"))
    if params["flagdneu"]:
        params.update(read_nml_group(fn, "neutral_perturb"))
    if params["flagprecfile"]:
        params.update(read_nml_group(fn, "precip"))
    if params["flagE0file"]:
        params.update(read_nml_group(fn, "efield"))
    if params["flagglow"]:
        params.update(read_nml_group(fn, "glow"))

    return params


def read_nml_group(fn: Path, group: str) -> typing.Dict[str, typing.Any]:
    """ read a group from an .nml file """

    raw: typing.Dict[str, typing.Any] = {}

    with fn.open("r") as f:
        for line in f:
            if not line.lstrip().startswith(f"&{group}"):
                continue
            for line in f:
                if line.strip().startswith("/"):
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
        raise KeyError(f"did not find group {group} in {fn}")

    return parse_group(raw, group)


def parse_group(raw: typing.Dict[str, str], group: str) -> typing.Dict[str, typing.Any]:
    """ this is Gemini-specific """

    if group == "base":
        P = {
            "t0": datetime(int(raw["ymd"][0]), int(raw["ymd"][1]), int(raw["ymd"][2])) + timedelta(seconds=float(raw["UTsec0"])),
            "tdur": timedelta(seconds=float(raw["tdur"])),
            "dtout": timedelta(seconds=float(raw["dtout"])),
            "f107a": float(raw["activ"][0]),
            "f107": float(raw["activ"][1]),
            "Ap": float(raw["activ"][2]),
            "tcfl": float(raw["tcfl"]),
            "Teinf": float(raw["Teinf"]),
        }
        for k in ("potsolve", "flagperiodic", "flagoutput", "flagcap", "flagglow", "flagE0file", "flagdneu", "flagprecfile"):
            try:
                P[k] = int(raw[k])
            except KeyError:
                P[k] = 0
        for k in ("indat_file", "indat_grid", "indat_size"):
            P[k] = make_abspath(raw[k])
    elif group == "neutral_perturb":
        P = {"interptype": int(raw["interptype"]), "sourcedir": make_abspath(raw["source_dir"])}

        for k in ("sourcemlat", "sourcemlon", "dtneu", "dxn", "drhon", "dzn"):
            try:
                P[k] = float(raw[k])
            except KeyError:
                P[k] = NaN
    elif group == "precip":
        P = {"dtprec": float(raw["dtprec"]), "precdir": make_abspath(raw["prec_dir"])}
    elif group == "efield":
        P = {"dtE0": float(raw["dtE0"]), "E0dir": make_abspath(raw["E0_dir"])}
    elif group == "glow":
        P = {"dtglow": float(raw["dtglow"]), "dtglowout": float(raw["dtglowout"])}
    else:
        raise ValueError(f"not sure how to handle group {group}")

    return P


def make_abspath(fn: str) -> Path:
    """ if relative path, make absolute based on this file """
    fn = Path(fn).expanduser()
    if not fn.is_absolute():
        fn = R / fn
    return fn


def read_ini(fn: Path) -> typing.Dict[str, typing.Any]:
    """ parse .ini file (legacy) """

    P: typing.Dict[str, typing.Any] = {}

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
