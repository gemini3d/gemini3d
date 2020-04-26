#!/usr/bin/env python3
import pytest
from pathlib import Path
from datetime import datetime, timedelta

import gemini3d.config as config

R = Path(__file__).resolve().parent
Rc = R / "config/test2d_fang"


def test_nml_bad(tmp_path):
    blank = tmp_path / "foo"
    blank.touch()
    with pytest.raises(KeyError):
        config.read_nml_group(blank, "base")

    blank.write_text(
        """
&base
 t0 =
/
"""
    )
    with pytest.raises(KeyError):
        config.read_nml_group(blank, "base")


@pytest.mark.parametrize("group", ["base", ("base", "flags", "files", "precip", "efield")])
def test_nml_group(group):

    params = config.read_nml_group(Rc / "config.nml", group)
    if "base" in group:
        assert params["t0"] == datetime(2013, 2, 20, 5)

    if "files" in group:
        assert params["format"] == "h5"

    if "precip" in group:
        assert params["dtprec"] == timedelta(seconds=5)

    if "efield" in group:
        assert params["dtE0"] == timedelta(seconds=1)


@pytest.mark.parametrize("filename", [Rc, Rc / "config.nml", R / "config/config_example.ini"], ids=["path", "nml", "ini"])
def test_read_config(filename):
    params = config.read_config(filename)
    assert params["t0"] == datetime(2013, 2, 20, 5)


if __name__ == "__main__":
    pytest.main([__file__])
