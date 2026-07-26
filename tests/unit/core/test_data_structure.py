from dataclasses import fields

import pytest

from process.core.model import DataStructure

DATA_STRUCTURES = {f.name: f.type() for f in fields(DataStructure)}


@pytest.mark.parametrize(
    "data_structure_class",
    DATA_STRUCTURES.values(),
    ids=DATA_STRUCTURES.keys(),
)
def test_arbitrary_assigns_disallowed(data_structure_class):
    """A test that checks arbitrary data cannot be assigned to a data structure.

    I.e. that the data class uses slots.
    """
    with pytest.raises(AttributeError):
        data_structure_class.my_property_that_doesnt_exist = 42
