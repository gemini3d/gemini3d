#!/usr/bin/env python3
import pytest
from pathlib import Path
from datetime import datetime, timedelta

import gemini3d.config as config

Rc = Path(__file__).resolve().parents[2] / "unit_tests/config"


def test_nml_bad(tmp_path):
    blank = tmp_path / "foo"
    blank.touch()
    with pytest.raises(KeyError):
        config.read_namelist(blank, "base")

    blank.write_text(
        """
&base
 t0 =
/
"""
    )
    with pytest.raises(KeyError):
        config.read_namelist(blank, "base")


@pytest.mark.parametrize("group", ["base", ("base", "flags", "files", "precip", "efield")])
def test_namelist(group):

    assert config.namelist_exists(Rc / "config_example.nml", "base")


@pytest.mark.parametrize("namelist", ["base", "flags", "files", "precip", "efield"])
def test_nml_namelist(namelist):

    params = config.read_namelist(Rc / "config_example.nml", namelist)
    if "base" in namelist:
        assert params["t0"] == datetime(2013, 2, 20, 5)

    if "files" in namelist:
        assert params["format"] == "h5"

    if "precip" in namelist:
        assert params["dtprec"] == timedelta(seconds=5)

    if "efield" in namelist:
        assert params["dtE0"] == timedelta(seconds=1)


@pytest.mark.parametrize(
    "filename",
    [Rc / "config_example.nml", Rc / "config_example.ini"],
    ids=["nml", "ini"],
)
def test_read_config(filename):
    params = config.read_config(filename)
    assert params["t0"] == datetime(2013, 2, 20, 5)


if __name__ == "__main__":
    pytest.main([__file__])
