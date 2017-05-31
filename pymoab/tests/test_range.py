from pymoab import core
from pymoab import types
from pymoab.rng import Range
from pymoab.rng import intersect, subtract, unite
from pymoab import types
from driver import test_driver, CHECK_EQ, CHECK_NOT_EQ
import numpy as np

def test_range():

    mb = core.Core()
    coord = np.array((1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),dtype='float64')
    verts = mb.create_vertices(coord)
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    data = np.array((1,2,3,4,5,6))
    mb.tag_set_data(test_tag,verts,data)

    dum = 0
    for v in verts:
        dum += 1
        if dum > 100: break
    CHECK_EQ(dum,len(verts))

    #test creation
    verts_range_as_list = list(verts)
    new_range = Range(verts_range_as_list)
    CHECK_EQ(len(new_range),len(verts))
    for eh_orig,eh_new in zip(verts,new_range):
        CHECK_EQ(eh_new,eh_orig)
    
    #test slicing/indexing
    CHECK_EQ(verts[-1],verts[5])
    CHECK_EQ(verts[-2],verts[4])

    sliced_range = verts[0:2]
    CHECK_EQ(len(sliced_range),2)
    CHECK_EQ(sliced_range[0],verts[0])
    CHECK_EQ(sliced_range[1],verts[1])

    sliced_range = verts[0:5:2]
    CHECK_EQ(len(sliced_range),3)
    CHECK_EQ(sliced_range[0],verts[0])
    CHECK_EQ(sliced_range[1],verts[2])
    CHECK_EQ(sliced_range[2],verts[4])

    sliced_range = verts[0:]
    CHECK_EQ(len(sliced_range),len(verts))
    for eh_orig,eh_new in zip(verts,sliced_range):
        CHECK_EQ(eh_new,eh_orig)
        
    sliced_range = verts[::]
    CHECK_EQ(len(sliced_range),len(verts))
    for eh_orig,eh_new in zip(verts,sliced_range):
        CHECK_EQ(eh_new,eh_orig)

    sliced_range = verts[-2:]
    CHECK_EQ(len(sliced_range),2)
    CHECK_EQ(sliced_range[0],verts[4])
    CHECK_EQ(sliced_range[1],verts[5])

    sliced_range = verts[-4:-2]
    CHECK_EQ(len(sliced_range),2)
    CHECK_EQ(sliced_range[0],verts[2])
    CHECK_EQ(sliced_range[1],verts[3])

    #test copy
    verts_copy = Range(verts)
    CHECK_EQ(len(verts_copy),len(verts))
    for eh_orig,eh_new in zip(verts,verts_copy):
        CHECK_EQ(eh_orig,eh_new)
    
    first_handle = verts_copy.pop_front()
    CHECK_EQ(len(verts_copy),len(verts)-1)
    CHECK_EQ(first_handle,verts[0])
             
    last_handle = verts_copy.pop_back()
    CHECK_EQ(len(verts_copy),len(verts)-2)
    CHECK_EQ(last_handle,verts[-1])

    verts_copy = Range(verts)
    handle = verts_copy[0]
    verts_copy.erase(handle)
    CHECK_EQ(len(verts_copy),len(verts)-1)
    assert(handle not in verts_copy)

    verts_copy.clear()
    CHECK_EQ(len(verts_copy),0)
    try:
        verts_copy[0]
    except StopIteration:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

def test_range_methods():
    mb = core.Core()
    coord = np.array((1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),dtype='float64')
    range_a = mb.create_vertices(coord)

    coord = np.array((2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7),dtype='float64')
    range_b = mb.create_vertices(coord)

    CHECK_EQ(range_a.all_of_dimension(0),True)    
    CHECK_EQ(range_b.all_of_dimension(0),True)

    CHECK_EQ(range_a.all_of_dimension(1),False)    
    CHECK_EQ(range_b.all_of_dimension(1),False)

    CHECK_EQ(range_a.num_of_dimension(0),range_a.size())
    CHECK_EQ(range_b.num_of_dimension(0),range_b.size())

    CHECK_EQ(range_a.num_of_dimension(1),0)
    CHECK_EQ(range_b.num_of_dimension(1),0)

    CHECK_EQ(range_a.num_of_type(types.MBVERTEX),range_a.size())
    CHECK_EQ(range_b.num_of_type(types.MBVERTEX),range_b.size())

    range_intersect = intersect(range_a, range_b)
    CHECK_EQ(range_intersect.size(),0)

    range_unite = unite(range_a, range_b)
    CHECK_EQ(range_unite.size(),12)

    range_subtract = subtract(range_a,range_b)
    CHECK_EQ(range_subtract.size(),range_a.size())
        
    range_a.erase(range_a[0])
    CHECK_EQ(range_a.size(),5)

    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    CHECK_EQ(all_verts.size(),12)

    range_intersect = intersect(all_verts, range_a)
    CHECK_EQ(range_intersect.size(),5)

    range_intersect = intersect(all_verts, range_b)
    CHECK_EQ(range_intersect.size(),6)

    range_unite = unite(all_verts, range_a)
    CHECK_EQ(range_unite.size(),12)

    range_unite = unite(all_verts, range_b)
    CHECK_EQ(range_unite.size(),12)

    range_subtract = subtract(all_verts, range_a)
    CHECK_EQ(range_subtract.size(),7)

    range_subtract = subtract(all_verts, range_b)
    CHECK_EQ(range_subtract.size(),6)

    range_a.merge(range_b)
    CHECK_EQ(range_a.size(),11)

if __name__ == "__main__":
    tests = [test_range,
             test_range_methods]
    test_driver(tests)
