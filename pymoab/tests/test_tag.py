
from pymoab import core
from pymoab import types

def test_tag_properties():
    mb = core.Core()

    tag_size = 16
    test_tag = mb.tag_get_handle("Test",tag_size,types.MB_TYPE_INTEGER,True)

    assert test_tag.get_length() == tag_size
