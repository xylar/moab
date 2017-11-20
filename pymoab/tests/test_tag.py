
from pymoab import core
from pymoab import types
from driver import test_driver, CHECK_EQ

def test_tag_properties():
    mb = core.Core()

    tag_size = 16
    test_tag = mb.tag_get_handle("Test",tag_size,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)

    CHECK_EQ(test_tag.get_length(), tag_size)

def test_tag_conventions():

    CHECK_EQ(type(types.MATERIAL_SET_TAG_NAME)   , str)
    CHECK_EQ(type(types.DIRICHLET_SET_TAG_NAME)  , str)
    CHECK_EQ(type(types.NEUMANN_SET_TAG_NAME)    , str)
    CHECK_EQ(type(types.HAS_MID_NODES_TAG_NAME)  , str)
    CHECK_EQ(type(types.GEOM_DIMENSION_TAG_NAME) , str)
    CHECK_EQ(type(types.MESH_TRANSFORM_TAG_NAME) , str)
    CHECK_EQ(type(types.GLOBAL_ID_TAG_NAME)      , str)
    CHECK_EQ(type(types.CATEGORY_TAG_NAME)       , str)
    CHECK_EQ(type(types.NAME_TAG_NAME)           , str)
    
    CHECK_EQ(type(types.CATEGORY_TAG_SIZE)       , int)    
    CHECK_EQ(type(types.NAME_TAG_SIZE)           , int)

if __name__ == "__main__":
    tests = [test_tag_properties, test_tag_conventions]
    test_driver(tests)


    
