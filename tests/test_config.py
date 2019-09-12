#!/usr/bin/env python3
import pytest
from pytest import approx
from pathlib import Path
from datetime import datetime

import gemini.config as config

R = Path(__file__).resolve().parents[1]
nml = R / "initialize/test2d_fang/config.nml"


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


@pytest.mark.parametrize("group", ["base", "precip", "efield"])
def test_nml_group(group):
    params = config.read_nml_group(nml, group)
    if group == "base":
        assert params["t0"] == datetime(2013, 2, 20, 5)
    elif group == "precip":
        assert params["dtprec"] == approx(5.0)
    elif group == "efield":
        assert params["dtE0"] == approx(1.0)


def test_nml_read():
    params = config.read_config(nml)
    assert params["t0"] == datetime(2013, 2, 20, 5)


if __name__ == "__main__":
    pytest.main([__file__])
