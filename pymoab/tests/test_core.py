from pymoab import core
from pymoab import types
from pymoab.scd import ScdInterface
from pymoab.hcoord import HomCoord
from subprocess import call
import numpy as np
import os

def test_load_mesh():
    mb = core.Core()
    mb.load_file("cyl_grps.h5m")

def test_write_mesh():
    mb = core.Core()
    mb.create_vertices(np.ones(3))
    mb.write_file("outfile.h5m")
    assert os.path.isfile("outfile.h5m")

def test_get_tag():
    mb = core.Core()
    #Create new tag
    tag = mb.tag_get_handle("Test",1,types.MB_TYPE_DOUBLE,True)
    #Create tag with speicified storage type
    tag = mb.tag_get_handle("Test1",1,types.MB_TYPE_DOUBLE,True,types.MB_TAG_SPARSE)
    #Query for exisiting tags
    ret_tag = mb.tag_get_handle("Test")
    ret_tag = mb.tag_get_handle("Test1")
    #Query for unknown tag, without raising exception
    tagtag = mb.tag_get_handle("Fake Tag", exceptions=(types.MB_TAG_NOT_FOUND,))

    #A bunch of invalid tag queries which should fail
    try:
        tag = mb.tag_get_handle("Fake Tag")
    except:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

    try:
        tag = mb.tag_get_handle("Fake Tag",1)
    except:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

    try:
        tag = mb.tag_get_handle("Fake Tag",None,types.MB_TYPE_DOUBLE)
    except:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

def test_integer_tag():
    mb = core.Core()
    vh = vertex_handle(mb)
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_INTEGER,True)
    test_val = 4
    test_tag_data = np.array((test_val,))
    mb.tag_set_data(test_tag, vh, test_tag_data)
    data = mb.tag_get_data(test_tag, vh)

    assert len(data) == 1
    assert data[0] == test_val
    assert data.dtype == 'int32'

def test_double_tag():
    mb = core.Core()
    vh = vertex_handle(mb)
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_DOUBLE,True)
    test_val = 4.4
    test_tag_data = np.array((test_val))
    mb.tag_set_data(test_tag, vh, test_tag_data)
    data = mb.tag_get_data(test_tag, vh)

    assert len(data) == 1
    assert data[0] == test_val
    assert data.dtype == 'float64'

    #a couple of tests that should fail
    test_tag = mb.tag_get_handle("Test1",1,types.MB_TYPE_DOUBLE,True)
    test_val = 4.4
    test_tag_data = np.array((test_val),dtype='float32')
    try:
        mb.tag_set_data(test_tag, vh, test_tag_data)
    except AssertionError:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

    test_tag = mb.tag_get_handle("Test2",1,types.MB_TYPE_DOUBLE,True)
    test_val = 4.4
    test_tag_data = np.array((test_val),dtype='int32')
    try:
        mb.tag_set_data(test_tag, vh, test_tag_data)
    except AssertionError:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

def test_opaque_tag():
    mb = core.Core()
    vh = vertex_handle(mb)
    tag_length = 6
    test_tag = mb.tag_get_handle("Test",tag_length,types.MB_TYPE_OPAQUE,True)
    test_val = 'four'
    test_tag_data = np.array((test_val,))
    mb.tag_set_data(test_tag, vh, test_tag_data)
    data = mb.tag_get_data(test_tag, vh)

    assert len(data) == 1
    assert data.nbytes == tag_length
    assert data[0] == test_val
    assert data.dtype == '|S' + str(tag_length)

def test_create_meshset():
    mb = core.Core()
    msh = mb.create_meshset()
    vh = vertex_handle(mb)
    mb.add_entities(msh,vh)
    
def vertex_handle(core):
    """Convenience function for getting an arbitrary vertex element handle."""
    coord = np.array((1,1,1),dtype='float64')
    vert = core.create_vertices(coord)
    vert_copy = np.array((vert[0],),dtype='uint64')
    return vert_copy

def test_create_elements():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    assert 3 == verts.size()
    #create elements
    verts = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)
    assert 1 == len(tris)
    #check that the element is there via GLOBAL_ID tag
    global_id_tag = mb.tag_get_handle("GLOBAL_ID",1,types.MB_TYPE_INTEGER,True)
    tri_id = mb.tag_get_data(global_id_tag, tris)
    assert 1 == len(tri_id)
    assert 0 == tri_id[0]

def test_range():

    mb = core.Core()
    coord = np.array((1,1,1),dtype='float64')
    vert = mb.create_vertices(coord)
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_INTEGER,True)
    data = np.array((1,))
    mb.tag_set_data(test_tag,vert,data)

    dum = 0
    for v in vert:
        dum += 1
        if dum > 100: break
    assert vert.size() == dum
    assert len(vert) == vert.size()


def test_tag_failures():
    mb = core.Core()
    coord = np.array((1,1,1),dtype='float64')
    verts = mb.create_vertices(coord)
    verts_illicit_copy = np.array((verts[0],),dtype='uint32')
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_INTEGER,True)
    data = np.array((1,))

    #this operation should fail due to the entity handle data type
    mb.tag_set_data(test_tag,verts,data)
    try:
        mb.tag_set_data(test_tag, verts_illicit_copy, data)
    except AssertionError:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError


    global_id_tag = mb.tag_get_handle("GLOBAL_ID",1,types.MB_TYPE_INTEGER,True)
    #so should this one
    try:
        tri_id = mb.tag_get_data(global_id_tag, verts_illicit_copy)
    except AssertionError:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

def test_adj():

    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    #create elements
    verts = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)
    #get the adjacencies of the triangle of dim 1 (should return the vertices)
    adjs = mb.get_adjacencies(tris, 0, False)
    assert 3 is adjs.size()

    #check that the entities are of the correct type
    for adj in adjs:
        type = mb.type_from_handle(adj)
        assert type is types.MBVERTEX

    #now get the edges and ask MOAB to create them for us
    adjs = mb.get_adjacencies(tris, 1, True)
    assert 3 is adjs.size()

    for adj in adjs:
        type = mb.type_from_handle(adj)
        assert type is types.MBEDGE

def test_meshsets():
    mb = core.Core()
    parent_set = mb.create_meshset()
    for i in range(5):
        a = mb.create_meshset()
        mb.add_child_meshset(parent_set,a)
    children = mb.get_child_meshsets(parent_set)
    assert children.size() is 5

def test_rs():
    mb = core.Core()
    rs = mb.get_root_set()
    assert rs == 0

def test_get_coords():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    ret_coords = mb.get_coords(verts)
    print ret_coords
    for i in range(len(coords)):
            assert coords[i] ==  ret_coords[i]

def test_get_ents_by_type():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    ms = mb.create_meshset()
    mb.add_entities(ms,verts)
    ret_verts = mb.get_entities_by_type(ms,types.MBVERTEX, False)
    for i in range(verts.size()):
        assert verts[i] == ret_verts[i]

def test_get_ents_by_tnt():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    test_tag = mb.tag_get_handle("TestTag",1,types.MB_TYPE_INTEGER,True)
    mb.tag_set_data(test_tag,verts,np.array((0,1,2,)))

    test_tag1 = mb.tag_get_handle("TestTag1",1,types.MB_TYPE_DOUBLE,True)
    mb.tag_set_data(test_tag1,verts,np.array((3.0,4.0,5.0)))
    
    rs = mb.get_root_set()

    # any vertices tagged with test_tag with a value of 6 (none of them)
    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag,test_tag1)),np.array((0,3.0)))
    print entities.size()
    assert entities.size() == 1

    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag,)),np.array((1,)))
    print entities.size()
    assert entities.size() == 1
    
    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag,)),np.array((None,)))
    print entities.size()
    assert entities.size() == 3

    # intersection of any vertices with test_tag value 1 and any value for test_tag1
    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag,test_tag1,)),np.array((1,None,)))
    print entities.size()
    assert entities.size() == 1

    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag1,test_tag,)),np.array((None,1,)))
    print entities.size()
    assert entities.size() == 1

    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag1,test_tag,)),np.array((3.0,None,)))
    print entities.size()
    assert entities.size() == 1

    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag1,test_tag,)),np.array((None,None,)))
    print entities.size()
    assert entities.size() == 3
    
    # intersection of any vertices tagged with test_tag and test_tag1 for any value on either tag
    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag,test_tag1,)),np.array((None,None,)))
    print entities.size()
    assert entities.size() == 3

    # intersection of any vertices tagged with test_tag and test_tag1 with a value of 2
    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag,test_tag1,)),np.array((None,3.0,)))
    print entities.size()
    assert entities.size() == 1
    
    # intersection of any vertices with test_tag & test_tag1 and values of 6 & 10, respectively (none of them)
    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag,test_tag1,)),np.array((6,10.0,)))
    print entities.size()
    assert entities.size() == 0
    
    # intersection of any vertices with test_tag & test_tag1 tagged with values of 1 & 2 respectively (none of them)
    entities = mb.get_entities_by_type_and_tag(rs,types.MBVERTEX,np.array((test_tag,test_tag1,)),np.array((1,2,)))
    print entities.size()
    assert entities.size() == 0

    # any hex elements tagged with test_tag (no hex elements exist, there should be none)
    entities = mb.get_entities_by_type_and_tag(rs,types.MBHEX,np.array((test_tag,)),np.array((None,)))
    print entities.size()
    assert entities.size() == 0

def test_get_entities_by_handle():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    ms = mb.create_meshset()
    mb.add_entities(ms,verts)
    ret_verts = mb.get_entities_by_handle(ms, False)
    for i in range(verts.size()):
        assert verts[i] == ret_verts[i]

def test_get_entities_by_dimension():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    rs = mb.get_root_set()
    ret_verts = mb.get_entities_by_dimension(rs, 0)
    for i in range(verts.size()):
        assert verts[i] == ret_verts[i]
