
from pymoab import core
from pymoab import types
from driver import test_driver

def test_tag_properties():
    mb = core.Core()

    tag_size = 16
    test_tag = mb.tag_get_handle("Test",tag_size,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)

    assert test_tag.get_length() == tag_size

if __name__ == "__main__":
    tests = [test_tag_properties,]
    test_driver(tests)
