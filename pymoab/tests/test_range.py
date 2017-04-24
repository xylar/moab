from pymoab import core
from pymoab import types
from pymoab.rng import Range
from driver import test_driver
import numpy as np

def test_range():

    mb = core.Core()
    coord = np.array((1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),dtype='float64')
    verts = mb.create_vertices(coord)
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_INTEGER,True)
    data = np.array((1,2,3,4,5,6))
    mb.tag_set_data(test_tag,verts,data)

    dum = 0
    for v in verts:
        dum += 1
        if dum > 100: break
    assert len(verts) == dum

    #test creation
    verts_range_as_list = list(verts)
    new_range = Range(verts_range_as_list)
    assert len(new_range) == len(verts)
    for eh_orig,eh_new in zip(verts,new_range):
        assert eh_orig == eh_new
    
    #test slicing/indexing
    assert verts[-1] == verts[5]
    assert verts[-2] == verts[4]

    sliced_range = verts[0:2]
    assert len(sliced_range) == 2
    assert sliced_range[0] == verts[0]
    assert sliced_range[1] == verts[1]

    sliced_range = verts[0:5:2]
    assert len(sliced_range) == 3
    assert sliced_range[0] == verts[0]
    assert sliced_range[1] == verts[2]
    assert sliced_range[2] == verts[4]

    sliced_range = verts[0:]
    assert len(sliced_range) == len(verts)
    for eh_orig,eh_new in zip(verts,sliced_range):
        assert eh_orig == eh_new
        
    sliced_range = verts[::]
    assert len(sliced_range) == len(verts)
    for eh_orig,eh_new in zip(verts,sliced_range):
        assert eh_orig == eh_new

    sliced_range = verts[-2:]
    assert len(sliced_range) == 2
    assert sliced_range[0] == verts[4]
    assert sliced_range[1] == verts[5]

    sliced_range = verts[-4:-2]
    assert len(sliced_range) == 2
    assert sliced_range[0] == verts[2]
    assert sliced_range[1] == verts[3]

    #test copy
    verts_copy = Range(verts)
    assert len(verts_copy) == len(verts)
    for eh_orig,eh_new in zip(verts,verts_copy):
        assert eh_orig == eh_new
    
    first_handle = verts_copy.pop_front()
    assert len(verts_copy) == len(verts)-1
    assert first_handle == verts[0]
    last_handle = verts_copy.pop_back()
    assert len(verts_copy) == len(verts)-2
    assert last_handle == verts[-1]

    verts_copy = Range(verts)
    handle = verts_copy[0]
    verts_copy.erase(handle)
    assert len(verts_copy) == len(verts)-1
    assert handle not in verts_copy

    verts_copy.clear()
    assert len(verts_copy) == 0
    try:
        verts_copy[0]
    except StopIteration:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

if __name__ == "__main__":
    tests = [test_range,]
    test_driver(tests)
