"""Tests for InDat and INVariable classes"""

import pytest

from process.io.in_dat import INVariable


@pytest.mark.parametrize("value", ["1", "1.0", "1,2,3", "1.0, 2.0", "string"])
def test_invariable_equality(value):
    """A test to check equality between INVariable's"""
    name = "test"
    v_type = "Parameter"
    parameter_group = "test_group"

    v1 = INVariable(name, value, v_type, parameter_group, "")
    v2 = INVariable(name, value, v_type, parameter_group, "different comment")

    assert v1 == v2
