import functools
import typing as T
import re
from pathlib import Path
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

    file = get_config_filename(path)

    if file.suffix == ".nml":
        P = read_nml(file)
    elif file.suffix == ".ini":
        P = read_ini(file)
    else:
        raise ValueError(f"unsure how to read {file}")

    P["nml"] = file

    return P


def get_config_filename(path: Path) -> Path:
    """ given a path or config filename, return the full path to config file """

    path = Path(path).expanduser().resolve()

    if path.is_file():
        return path

    if path.is_dir():
        for p in (path, path / "inputs"):
            for suff in (".nml", ".ini"):
                for file in p.glob("config*" + suff):
                    if file.is_file():
                        return file

    raise FileNotFoundError(f"could not find config file in {path}")


def read_nml(fn: Path) -> T.Dict[str, T.Any]:
    """ parse .nml file
    for now we don't use the f90nml package, though maybe we will in the future.
    Just trying to keep Python prereqs reduced for this simple parsing.
    """

    fn = Path(fn).expanduser().resolve()

    if fn.is_dir():
        fn = fn / "config.nml"

    params = {}
    for n in ("base", "files", "flags", "setup"):
        params.update(read_namelist(fn, n))

    if namelist_exists(fn, "neutral_perturb"):
        params.update(read_namelist(fn, "neutral_perturb"))
    if namelist_exists(fn, "precip"):
        params.update(read_namelist(fn, "precip"))
    if namelist_exists(fn, "efield"):
        params.update(read_namelist(fn, "efield"))
    if namelist_exists(fn, "glow"):
        params.update(read_namelist(fn, "glow"))

    return params


def namelist_exists(fn: Path, namelist: str) -> bool:
    """ determines if a namelist exists in a file """

    pat = re.compile(r"^\s*&(" + namelist + ")$")

    with fn.open("r") as f:
        for line in f:
            if pat.match(line) is not None:
                return True

    return False


def read_namelist(fn: Path, namelist: str) -> T.Dict[str, T.Any]:
    """ read a namelist from an .nml file """

    raw: T.Dict[str, T.Sequence[str]] = {}
    nml_pat = re.compile(r"^\s*&(" + namelist + r")")
    end_pat = re.compile(r"^\s*/\s*$")
    val_pat = re.compile(r"^\s*(\w+)\s*=\s*['\"]?([^!'\"]*)['\"]?")

    with fn.open("r") as f:
        for line in f:
            if not nml_pat.match(line):
                continue

            for line in f:
                if end_pat.match(line):
                    # end of namelist
                    return parse_namelist(raw, namelist)
                val_mat = val_pat.match(line)
                if not val_mat:
                    continue

                key, vals = val_mat.group(1), val_mat.group(2).split(",")
                raw[key] = vals[0] if len(vals) == 1 else vals

    raise KeyError(f"did not find Namelist {namelist} in {fn}")


def parse_namelist(raw: T.Dict[str, T.Any], namelist: str) -> T.Dict[str, T.Any]:
    """
    this is Gemini-specific
    don't resolve absolute paths here because that assumes same machine
    """

    P: T.Dict[str, T.Any] = {}

    if namelist == "base":
        P["t0"] = datetime(int(raw["ymd"][0]), int(raw["ymd"][1]), int(raw["ymd"][2])) + timedelta(
            seconds=float(raw["UTsec0"])
        )
        P["tdur"] = timedelta(seconds=float(raw["tdur"]))
        P["dtout"] = timedelta(seconds=float(raw["dtout"]))
        P["f107a"] = float(raw["activ"][0])
        P["f107"] = float(raw["activ"][1])
        P["Ap"] = float(raw["activ"][2])
        P["tcfl"] = float(raw["tcfl"])
        P["Teinf"] = float(raw["Teinf"])
    elif namelist == "flags":
        for k in raw:
            P[k] = int(raw[k])
    elif namelist == "files":
        for k in ("indat_file", "indat_grid", "indat_size"):
            P[k] = Path(raw[k])

        if "file_format" in raw:
            P["format"] = raw["file_format"]
        else:
            # defaults to type of input
            P["format"] = P["indat_size"].suffix[1:]

        if "realbits" in raw:
            P["realbits"] = int(raw["realbits"])
        else:
            if P["format"] in ("raw", "dat"):
                P["realbits"] = 64
            else:
                P["realbits"] = 32
    elif namelist == "setup":
        P["alt_scale"] = list(map(float, raw["alt_scale"]))

        for k in ("lxp", "lyp"):
            P[k] = int(raw[k])
        for k in (
            "glat",
            "glon",
            "xdist",
            "ydist",
            "alt_min",
            "alt_max",
            "Bincl",
            "nmf",
            "nme",
            "precip_latwidth",
            "precip_lonwidth",
            "Qprecip",
            "Qprecip_background",
            "E0precip",
            "Etarg",
            "Jtarg",
            "Efield_latwidth",
            "Efield_lonwidth",
            # "Eflagdirich",  # future
        ):
            if k in raw:
                P[k] = float(raw[k])
        for k in ("eqdir",):
            if k in raw:
                P[k] = Path(raw[k])
    elif namelist == "neutral_perturb":
        P["interptype"] = int(raw["interptype"])
        P["sourcedir"] = Path(raw["source_dir"])

        for k in ("sourcemlat", "sourcemlon", "dtneu", "dxn", "drhon", "dzn"):
            try:
                P[k] = float(raw[k])
            except KeyError:
                P[k] = NaN
    elif namelist == "precip":
        P["dtprec"] = timedelta(seconds=float(raw["dtprec"]))
        P["precdir"] = Path(raw["prec_dir"])
    elif namelist == "efield":
        P["dtE0"] = timedelta(seconds=float(raw["dtE0"]))
        P["E0dir"] = Path(raw["E0_dir"])
    elif namelist == "glow":
        P["dtglow"] = timedelta(seconds=float(raw["dtglow"]))
        P["dtglowout"] = float(raw["dtglowout"])

    if not P:
        raise ValueError(f"Not sure how to parse NML namelist {namelist}")

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
