from pymoab import core
from pymoab import types
from pymoab.rng import Range
from pymoab.scd import ScdInterface
from pymoab.hcoord import HomCoord
from subprocess import call
from driver import test_driver, CHECK, CHECK_EQ, CHECK_NOT_EQ, CHECK_ITER_EQ
import numpy as np
import os

bytes_per_char_ = np.array(["a"]).nbytes

def test_load_mesh():
    mb = core.Core()
    try:
        mb.load_file("cyl_grps.h5m")
    except:
        try:
            print( """
            WARNING: .h5m file load failed. If hdf5 support is enabled in this
            build there could be a problem.
            """)
            mb.load_file("cyl_grps.vtk")
        except:
            raise(IOError, "Failed to load MOAB file.")

    #load into file_set
    mb1 = core.Core()
    file_set = mb1.create_meshset()
    try:
        mb1.load_file("cyl_grps.h5m",file_set)
    except:
        try:
            print("""
            WARNING: .h5m file load failed. If hdf5 support is enabled in this
            build there could be a problem.
            """)
            mb1.load_file("cyl_grps.vtk",file_set)
        except:
            raise(IOError, "Failed to load MOAB file.")

    ents = mb1.get_entities_by_type(file_set,types.MBMAXTYPE)
    CHECK_NOT_EQ(len(ents),0)

def test_write_mesh():
    mb = core.Core()
    mb.create_vertices(np.ones(3))

    try:
        mb.write_file("outfile.h5m")
        assert os.path.isfile("outfile.h5m")
    except:
        try:
            print("""
            WARNING: .h5m file write failed. If hdf5 support is enabled in this
            build there could be a problem.
            """)
            mb.write_file("outfile.vtk")
            assert os.path.isfile("outfile.vtk")
        except:
            raise(IOError, "Failed to write MOAB file.")


def test_delete_mesh():
    mb = core.Core()
    mb.create_vertices(np.ones(9))
    rs = mb.get_root_set()
    ents = mb.get_entities_by_handle(rs)
    CHECK_EQ(len(ents),3)
    # now delete all mesh entities
    mb.delete_mesh()
    ents = mb.get_entities_by_handle(rs)
    CHECK_EQ(len(ents),0)

def test_get_tag():
    mb = core.Core()
    #Create new tag
    tag = mb.tag_get_handle("Test",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    #Create tag with speicified storage type
    tag = mb.tag_get_handle("Test1",1,types.MB_TYPE_DOUBLE,types.MB_TAG_SPARSE,True)
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
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

    try:
        tag = mb.tag_get_handle("Fake Tag",1)
    except:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

    try:
        tag = mb.tag_get_handle("Fake Tag",None,types.MB_TYPE_DOUBLE)
    except:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

    # attempt tag lookup with incomplete tag spceification)
    try:
        tag = mb.tag_get_handle("Bad Tag", 1, types.MB_TYPE_OPAQUE)
    except ValueError:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(ValueError)

    def_val_tag_name = "def_val_tag"
    def_value = [0.0,1.0]
    def_val_tag = mb.tag_get_handle(def_val_tag_name,2,types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True, default_value = def_value)

    #look up tag by name only
    tag = mb.tag_get_handle(def_val_tag_name)

    #look up tag with default value
    tag = mb.tag_get_handle(def_val_tag_name,2,types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, default_value = def_value)

    # shouldn't be able to change the default value
    try:
        tag = mb.tag_get_handle(def_val_tag_name,2,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE, True, default_value = [1.0,0.0])
    except:
        pass
    else:
        raise(AssertionError)

    # check the default value
    def_val_chck = mb.tag_get_default_value(def_val_tag)
    assert len(def_val_chck) == 2
    assert def_val_chck[0] == def_value[0]
    assert def_val_chck[1] == def_value[1]

    def_tag_length = mb.tag_get_length(def_val_tag)
    assert def_tag_length == 2

def test_integer_tag():
    mb = core.Core()
    vh = vertex_handle(mb)
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    test_val = 4
    test_tag_data = np.array((test_val,))
    mb.tag_set_data(test_tag, vh, test_tag_data)
    data = mb.tag_get_data(test_tag, vh)

    CHECK_EQ(len(data),1)
    CHECK_EQ(data[0],test_val)
    CHECK_EQ(data.dtype,'int32')

    #set tag data for single handle (non-iterable)
    new_test_value = 36
    mb.tag_set_data(test_tag, vh[0], new_test_value)
    data = mb.tag_get_data(test_tag, vh[0])

    CHECK_EQ(len(data),1)
    CHECK_EQ(new_test_value, data[0])
    CHECK_EQ(data.dtype, 'int32')

    tags = mb.tag_get_tags_on_entity(vh[0])
    assert len(tags) > 0

def test_double_tag():
    mb = core.Core()
    vh = vertex_handle(mb)
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    test_val = 4.4
    test_tag_data = np.array((test_val,))
    mb.tag_set_data(test_tag, vh, test_tag_data)
    data = mb.tag_get_data(test_tag, vh)

    CHECK_EQ(len(data),1)
    CHECK_EQ(data[0],test_val)
    CHECK_EQ(data.dtype,'float64')

    #a couple of tests that should fail
    test_tag = mb.tag_get_handle("Test1",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    test_val = 4.4
    test_tag_data = np.array((test_val),dtype='float32')
    try:
        mb.tag_set_data(test_tag, vh, test_tag_data)
    except AssertionError:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

    test_tag = mb.tag_get_handle("Test2",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    test_val = 4.4
    test_tag_data = np.array((test_val),dtype='int32')
    try:
        mb.tag_set_data(test_tag, vh, test_tag_data)
    except AssertionError:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

    #set tag data for single handle (non-iterable)
    new_test_value = 36.0
    mb.tag_set_data(test_tag, vh[0], new_test_value)
    data = mb.tag_get_data(test_tag, vh[0])

    CHECK_EQ(len(data),1)
    CHECK_EQ(new_test_value, data[0])
    CHECK_EQ(data.dtype, 'float64')

def test_opaque_tag():
    mb = core.Core()
    vh = vertex_handle(mb)
    tag_length = 6
    test_tag = mb.tag_get_handle("Test",tag_length,types.MB_TYPE_OPAQUE,types.MB_TAG_DENSE,True)
    test_val = 'four'
    test_tag_data = np.array((test_val,))
    mb.tag_set_data(test_tag, vh, test_tag_data)
    data = mb.tag_get_data(test_tag, vh)
    CHECK_EQ(len(data),1)
    CHECK_EQ(data.nbytes,tag_length*bytes_per_char_)
    CHECK_EQ(data[0],test_val)

def test_tag_list():
    mb = core.Core()
    vh = vertex_handle(mb)
    tag_length = 6
    test_tag = mb.tag_get_handle("Test",tag_length,types.MB_TYPE_OPAQUE,types.MB_TAG_DENSE,True)
    test_val = 'four'
    test_tag_data = np.array((test_val,))
    # convert vertex handle Range to a list
    vh = list(vh)
    mb.tag_set_data(test_tag, vh, test_tag_data)
    data = mb.tag_get_data(test_tag, vh)

    CHECK_EQ(len(data),1)
    CHECK_EQ(data.nbytes,tag_length*bytes_per_char_)
    CHECK_EQ(data[0],test_val)

def test_tag_delete():
    mb = core.Core()
    vh = vertex_handle(mb)
    test_tag = mb.tag_get_handle(
        "Test", 1, types.MB_TYPE_INTEGER, True, types.MB_TAG_SPARSE)
    test_val = 4
    test_tag_data = np.array((test_val,))
    mb.tag_set_data(test_tag, vh, test_tag_data)

    mb.tag_delete_data(test_tag, vh)
    try:
        mb.tag_get_data(test_tag, vh)
        raised = False
    except RuntimeError as e:
        er_val = e.args[0].error_value
        raised = True
    CHECK_EQ(raised, True)
    CHECK_EQ(er_val, types.MB_TAG_NOT_FOUND)

    mb.tag_delete(test_tag)
    try:
        mb.tag_get_data(test_tag, vh)
        raised = False
    except RuntimeError as e:
        er_val = e.args[0].error_value
        raised = True
    CHECK_EQ(raised, True)
    CHECK_EQ(er_val, types.MB_TAG_NOT_FOUND)

def test_create_meshset():
    mb = core.Core()
    msh = mb.create_meshset()
    vh = vertex_handle(mb)
    mb.add_entities(msh,vh)

def vertex_handle(core):
    """Convenience function for getting an arbitrary vertex element handle."""
    coord = np.array((1,1,1),dtype='float64')
    vert = core.create_vertices(coord)
    return vert

def test_create_elements():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    CHECK_EQ(len(verts),3)
    #create elements
    verts = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)
    CHECK_EQ(len(tris),1)
    #check that the element is there via GLOBAL_ID tag
    global_id_tag = mb.tag_get_handle("GLOBAL_ID",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    tri_id = mb.tag_get_data(global_id_tag, tris)
    CHECK_EQ(len(tri_id),1)
    CHECK_EQ(tri_id[0],0)



def test_tag_failures():
    mb = core.Core()

    coord = np.array((1,1,1),dtype='float64')
    msh = mb.create_meshset()    
    msh_illicit_copy = np.array((msh,),dtype='uint32')
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    data = np.array((1,))

    #this operation should fail due to the entity handle data type
    mb.tag_set_data(test_tag,np.array([msh,]),data)
    try:
        mb.tag_set_data(test_tag, msh_illicit_copy, data)
    except ValueError:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(ValueError)


    coord = np.array((1,1,1),dtype='float64')
    verts = mb.create_vertices(coord)
    verts_illicit_copy = np.array((verts[0],),dtype='uint32')
    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    data = np.array((1,))

    #this operation should fail due to the entity handle data type
    mb.tag_set_data(test_tag,verts,data)
    try:
        mb.tag_set_data(test_tag, verts_illicit_copy, data)
    except ValueError:
        pass
    else:
        pass
        print("Shouldn't be here. Test fails.")
        raise(ValueError)


    global_id_tag = mb.tag_get_handle("GLOBAL_ID",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    #so should this one
    try:
        tri_id = mb.tag_get_data(global_id_tag, verts_illicit_copy)
    except ValueError:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(ValueError)

def test_adj():

    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    #create elements
    conn = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,conn)
    #get the adjacencies of the triangle of dim 1 (should return the vertices)
    adjs = mb.get_adjacencies(tris, 0, False)
    CHECK_EQ(len(adjs),3)
    CHECK(adjs.all_of_type(types.MBVERTEX))

    
    #check that the entities are of the correct type
    for adj in adjs:
        type = mb.type_from_handle(adj)
        assert type is types.MBVERTEX

    #now get the edges and ask MOAB to create them for us
    adjs = mb.get_adjacencies(tris, 1, True)
    CHECK_EQ(len(adjs),3)

    CHECK(adjs.all_of_type(types.MBEDGE))
          
    adjs = mb.get_adjacencies(tris[0], 0, False)
    CHECK_EQ(len(adjs),3)

    # create another triangle (with reverse normal)
    new_coords = np.array((2,2,2,3,5,3), dtype = 'float64')
    new_verts = mb.create_vertices(new_coords)
    
    # create a new triangle that shares a vertex with the first
    new_tri_conn = np.array(((verts[2], new_verts[0], new_verts[1]),) , dtype = 'uint64')
    tris.merge(mb.create_elements(types.MBTRI, new_tri_conn))
    # confirm that we can get the adjacency intersection of the two triangles
    adjs = mb.get_adjacencies(tris, 0, False)
    CHECK_EQ(len(adjs), 1)
    CHECK(adjs.all_of_type(types.MBVERTEX))
    
    adjs = mb.get_adjacencies(tris, 1, False)
    CHECK_EQ(len(adjs), 0)
    CHECK(adjs.all_of_type(types.MBEDGE))

    # now check that we can get the union of the two triangle adjacencies
    adjs = mb.get_adjacencies(tris, 0, False, types.UNION)
    CHECK_EQ(len(adjs), 5)
    CHECK(adjs.all_of_type(types.MBVERTEX))

    # sanity check for number of edges
    adjs = mb.get_adjacencies(tris[1], 1, False, types.UNION)
    CHECK_EQ(len(adjs), 0)    
    
def test_get_conn():

    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    #create elements
    verts = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)
    #get the adjacencies of the triangle of dim 1 (should return the vertices)
    conn = mb.get_connectivity(tris, 0, False)
    CHECK_EQ(len(conn),3)

    #check that the entities are of the correct type
    for c in conn:
        type = mb.type_from_handle(c)
        assert type is types.MBVERTEX

    conn = mb.get_connectivity(tris[0])
    CHECK_EQ(len(conn),3)
    CHECK_EQ(conn, verts)

    conn = mb.get_connectivity(Range(tris))
    CHECK_EQ(len(conn),3)
    CHECK_EQ(conn, verts)

    msh = mb.create_meshset()
    
    try:
        mb.get_connectivity(msh)
    except IndexError:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(IndexErrorx)
    
def test_type_from_handle():
    mb = core.Core()
    # create vertices
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    # create a triangle from the vertices
    verts_as_array = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts_as_array)
    # check that triangle is the correct type
    entity_type = mb.type_from_handle(tris[0])
    assert entity_type == types.MBTRI
    # check that vertex is the correct type
    entity_type = mb.type_from_handle(verts[0])
    assert entity_type == types.MBVERTEX
    # check that an entityset is the correct type
    rs = mb.get_root_set()
    entity_type = mb.type_from_handle(rs)
    assert entity_type == types.MBENTITYSET

def test_meshsets():
    mb = core.Core()
    parent_set = mb.create_meshset()
    for i in range(5):
        a = mb.create_meshset()
        mb.add_child_meshset(parent_set,a)
    children = mb.get_child_meshsets(parent_set)
    CHECK_EQ(len(children),5)

def test_rs():
    mb = core.Core()
    rs = mb.get_root_set()
    assert rs == 0

def test_get_coords():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    ret_coords = mb.get_coords(verts)
    for i in range(len(coords)):
        CHECK_EQ(ret_coords[i],coords[i])

    # test for passing entity handle (non-iterable)
    ret_coords = mb.get_coords(verts[0])
    CHECK_EQ(ret_coords[0], coords[0])
    CHECK_EQ(ret_coords[1], coords[1])
    CHECK_EQ(ret_coords[2], coords[2])

def test_set_coords():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    ret_coords = mb.get_coords(verts)
    coords += 1.0 # Shift x/y/z coordinates by 1.0
    mb.set_coords(verts, coords)
    ret_coords2 = mb.get_coords(verts)
    for i in range(len(coords)):
        CHECK_EQ(ret_coords[i],ret_coords2[i]-1.0)

def test_get_ents_by_type():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    ms = mb.create_meshset()
    mb.add_entities(ms,verts)
    ret_verts = mb.get_entities_by_type(ms,types.MBVERTEX, False)
    for i in range(len(verts)):
        CHECK_EQ(ret_verts[i],verts[i])

def test_get_ents_by_tnt():

    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    int_test_tag = mb.tag_get_handle("IntegerTestTag",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    int_test_tag_values = np.array((0,1,2,))
    mb.tag_set_data(int_test_tag,verts,int_test_tag_values)

    dbl_test_tag = mb.tag_get_handle("DoubleTestTag",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    dbl_test_tag_values = np.array((3.0,4.0,5.0))
    mb.tag_set_data(dbl_test_tag,verts,dbl_test_tag_values)

    opaque_test_tag = mb.tag_get_handle("OpaqueTestTag",6,types.MB_TYPE_OPAQUE,types.MB_TAG_DENSE,True)
    opaque_test_tag_values = np.array(("Six","Seven","Eight",))
    mb.tag_set_data(opaque_test_tag,verts,opaque_test_tag_values)

    opaque_test_tag1 = mb.tag_get_handle("OpaqueTestTag1",6,types.MB_TYPE_OPAQUE,types.MB_TAG_DENSE,True)
    opaque_test_tag1_values = np.array(("Nine","Ten","Eleven",))
    mb.tag_set_data(opaque_test_tag1,verts,opaque_test_tag1_values)

    rs = mb.get_root_set()

    single_tag_test_cases = []
    ###INTEGER TAG TESTS###
    integer_tag_test_cases = [
        dict(tag_arr = [int_test_tag], value_arr = np.array([[1]]), expected_size = 1),        # existing value
        dict(tag_arr = [int_test_tag], value_arr = np.array([[16]]), expected_size = 0),       # nonexistant value
        dict(tag_arr = [int_test_tag], value_arr = np.array([[None]]), expected_size = 3) ]    # any value
    single_tag_test_cases += integer_tag_test_cases # add integer tag test cases to all tests
    ###DOUBLE TAG TESTS###
    double_tag_test_cases = [
        dict(tag_arr = [dbl_test_tag], value_arr = np.array([[4.0]]), expected_size = 1),      # existing value
        dict(tag_arr = [dbl_test_tag], value_arr = np.array([[16.0]]), expected_size = 0),     # nonexistant value
        dict(tag_arr = [dbl_test_tag], value_arr = np.array([[None]]), expected_size = 3) ]    # any value
    single_tag_test_cases += double_tag_test_cases # add double tag test cases to all tests
    ###OPAQUE TAG TESTS###
    opaque_tag_test_cases = [
        dict(tag_arr = [opaque_test_tag], value_arr = np.array([["Six"]]), expected_size = 1), # existing value
        dict(tag_arr = [opaque_test_tag], value_arr = np.array([["Ten"]]), expected_size = 0), # nonexistant value
        dict(tag_arr = [opaque_test_tag], value_arr = np.array([[None]]), expected_size = 3) ] # any value
    single_tag_test_cases += opaque_tag_test_cases # add opaque tag test cases to all tests

    for test_case in single_tag_test_cases:
        entities = mb.get_entities_by_type_and_tag(rs,
                                                   types.MBVERTEX,
                                                   test_case['tag_arr'],
                                                   test_case['value_arr'])
        CHECK_EQ(len(entities),test_case['expected_size'])

    try:
        entities = mb.get_entities_by_type_and_tag(rs,
                                                   types.MBVERTEX,
                                                   [int_test_tag,dbl_test_tag],
                                                   [[1],[4.0]])
    except:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

    # ###MIXED TAG TYPE TESTS###
    # mixed_tag_test_cases = [
    #     # existing values
    #     dict(tag_arr = [int_test_tag,opaque_test_tag],
    #          value_arr = np.array([[1],["Seven"]], dtype='O'), #dtype must be specified here
    #          expected_size = 1),
    #     # any and existing value
    #     dict(tag_arr = [dbl_test_tag,int_test_tag],
    #          value_arr = np.array([[None],[1]]),
    #          expected_size = 1),
    #     # any and existing value (None comes second)
    #     dict(tag_arr = [dbl_test_tag,int_test_tag],
    #          value_arr = np.array([[3.0],[None]]),
    #          expected_size = 1),
    #     # any values for both tags
    #     dict(tag_arr = [dbl_test_tag,int_test_tag],
    #          value_arr = np.array([[None],[None]]),
    #          expected_size = 3),
    #     # mixed types, both nonexisting
    #     dict(tag_arr = [int_test_tag,dbl_test_tag],
    #          value_arr = np.array([[6],[10.0]], dtype = 'O'),
    #          expected_size = 0),
    #     #mixed types, both existing but not on same entity
    #     dict(tag_arr = [int_test_tag,dbl_test_tag],
    #          value_arr = np.array([[1],[2.0]], dtype = 'O'),
    #          expected_size = 0),
    #     #mixed types, only one existing value
    #     dict(tag_arr = [int_test_tag,dbl_test_tag],
    #          value_arr = np.array([[1],[12.0]], dtype = 'O'),
    #          expected_size = 0),
    #     dict(tag_arr = [opaque_test_tag,opaque_test_tag1],
    #          value_arr = np.array([["Seven"],["Ten"]], dtype='O'),
    #          expected_size = 1)
    # ]

    # for test_case in mixed_tag_test_cases:
    #     entities = mb.get_entities_by_type_and_tag(rs,
    #                                                types.MBVERTEX,
    #                                                test_case['tag_arr'],
    #                                                test_case['value_arr'])
    #     CHECK_EQ(len(entities),test_case['expected_size'])

    ###VECTOR TAG TESTS###
    int_vec_test_tag = mb.tag_get_handle("IntegerVecTestTag",3,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    dbl_vec_test_tag = mb.tag_get_handle("DoubleVecTestTag",3,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)

    int_vec_test_tag_values = np.array([[0,1,2],[3,4,5],[6,7,8]])
    dbl_vec_test_tag_values = np.array([[9.0,10.0,11.0],[12.0,13.0,14.0],[15.0,16.0,17.0]])
    mb.tag_set_data(int_vec_test_tag,verts,int_vec_test_tag_values)
    mb.tag_set_data(dbl_vec_test_tag,verts,dbl_vec_test_tag_values)

    # # existing sets of values
    # entities = mb.get_entities_by_type_and_tag(rs,
    #                                            types.MBVERTEX,
    #                                            [int_vec_test_tag,dbl_vec_test_tag],
    #                                            np.array([[0,1,2],[9.0,10.0,11.0]],dtype='O'))
    # CHECK_EQ(len(entities),1)

    # #one non existant set of values, one existing
    # entities = mb.get_entities_by_type_and_tag(rs,
    #                                            types.MBVERTEX,
    #                                            [int_vec_test_tag,dbl_vec_test_tag],
    #                                            np.array([[0,1,2],[22.0,10.0,11.0]],dtype='O'))

    # CHECK_EQ(len(entities),0)

    # # any set of one tag
    # print "Running suspect test"
    # entities = mb.get_entities_by_type_and_tag(rs,
    #                                            types.MBVERTEX,
    #                                            [int_vec_test_tag,dbl_vec_test_tag],
    #                                            np.array([[None,None,None],[9.0,10.0,11.0]],dtype='O'))
    # CHECK_EQ(len(entities),1)

    # entities = mb.get_entities_by_type_and_tag(rs,
    #                                            types.MBVERTEX,
    #                                            [int_vec_test_tag,dbl_vec_test_tag],
    #                                            np.array([[None],[None]],dtype='O'))
    # CHECK_EQ(len(entities),3)


    # any set of one tag
    entities = mb.get_entities_by_type_and_tag(rs,
                                               types.MBVERTEX,
                                               dbl_vec_test_tag,
                                               np.array([9.0,10.0,11.0],dtype='O'))
    CHECK_EQ(len(entities),1)

    # any hex elements tagged with int_test_tag (no hex elements exist, there should be none)
    tag_test_vals = np.array([[None]])
    entities = mb.get_entities_by_type_and_tag(rs,
                                               types.MBHEX,
                                               np.array((int_test_tag,)),
                                               tag_test_vals)
    CHECK_EQ(len(entities) , 0)

def test_get_entities_by_handle():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    ms = mb.create_meshset()
    mb.add_entities(ms,verts)
    ret_verts = mb.get_entities_by_handle(ms, False)
    for i in range(len(verts)):
        CHECK_EQ(ret_verts[i],verts[i])

def test_get_entities_by_dimension():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    rs = mb.get_root_set()
    ret_verts = mb.get_entities_by_dimension(rs, 0)
    for i in range(len(verts)):
        CHECK_EQ(ret_verts[i],verts[i])

def test_parent_child():
    mb = core.Core()

    parent_set = mb.create_meshset()
    child_set = mb.create_meshset()

    mb.add_parent_meshset(child_set, parent_set)

    parent_sets = mb.get_parent_meshsets(child_set)
    CHECK_EQ(len(parent_sets),1)
    CHECK_EQ(parent_sets[0],parent_set)

    child_sets = mb.get_child_meshsets(parent_set)
    CHECK_EQ(len(child_sets),0)

    mb.add_child_meshset(parent_set,child_set)

    child_sets = mb.get_child_meshsets(parent_set)
    CHECK_EQ(len(child_sets),1)
    CHECK_EQ(child_sets[0],child_set)

    parent_set = mb.create_meshset()
    child_set = mb.create_meshset()

    parent_sets = mb.get_parent_meshsets(child_set)
    CHECK_EQ(len(parent_sets),0)
    child_sets = mb.get_child_meshsets(parent_set)
    CHECK_EQ(len(child_sets),0)

    mb.add_parent_child(parent_set,child_set)

    parent_sets = mb.get_parent_meshsets(child_set)
    CHECK_EQ(len(parent_sets),1)
    parent_sets[0] == parent_set
    child_sets = mb.get_child_meshsets(parent_set)
    CHECK_EQ(len(child_sets),1)
    child_sets[0] == child_set

def test_remove_ents():
    mb = core.Core()
    ms = mb.create_meshset()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    mb.add_entities(ms,verts)

    #make sure all vertices were added to the meshset
    entities = mb.get_entities_by_type(ms, types.MBVERTEX)
    CHECK_EQ(len(entities),3)
    CHECK_EQ(entities[0],verts[0])
    CHECK_EQ(entities[1],verts[1])
    CHECK_EQ(entities[2],verts[2])

    #remove first vertex from meshset
    vert_to_remove = np.array([verts[0],],dtype='uint64') #not proud of this...
    mb.remove_entities(ms,vert_to_remove)

    #make sure the right one got removed from the meshset
    entities = mb.get_entities_by_type(ms, types.MBVERTEX)
    CHECK_EQ(len(entities),2)
    CHECK_EQ(entities[0],verts[1])
    CHECK_EQ(entities[1],verts[2])

    #all verts should still be in the instance
    rs = mb.get_root_set()
    entities = mb.get_entities_by_type(rs, types.MBVERTEX)
    CHECK_EQ(len(entities),3)
    CHECK_EQ(entities[0],verts[0])
    CHECK_EQ(entities[1],verts[1])
    CHECK_EQ(entities[2],verts[2])

    #until it is deleted from the instance
    mb.delete_entities(vert_to_remove)
    entities = mb.get_entities_by_type(rs, types.MBVERTEX)
    CHECK_EQ(len(entities),2)
    CHECK_EQ(entities[0],verts[1])
    CHECK_EQ(entities[1],verts[2])

def test_iterables():

    mb = core.Core()
    # 1-D iterable
    coords = [0.,0.,0.,1.,0.,0.,1.,1.,1.]
    verts = mb.create_vertices(coords)
    CHECK_EQ(len(verts),3)
    # 2-D iterable w/ len 3 entries
    coords = [[0.,0.,0.],[1.,0.,0.],[1.,1.,1.],[0., 1., 0.]]
    verts = mb.create_vertices(coords)
    CHECK_EQ(len(verts),4)

    int_tag = mb.tag_get_handle("IntTag",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)

    #try to set data with bad array (contains int)
    int_data = [1, 2, 3, 4.0]

    try:
        mb.tag_set_data(int_tag,verts,int_data)
    except:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

    #now set with valid data and check
    int_data = [1, 2, 3, 4]

    mb.tag_set_data(int_tag,verts,int_data)

    return_data = mb.tag_get_data(int_tag,verts)

    CHECK_EQ(return_data[0],int_data[0])
    CHECK_EQ(return_data[1],int_data[1])
    CHECK_EQ(return_data[2],int_data[2])
    CHECK_EQ(return_data[3],int_data[3])

    #insert false vertex handle (not even correct type
    verts = [verts[0],23,verts[1]]
    try:
        mb.tag_set_data(int_tag,verts,int_data)
    except:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

    #insert correct type, but non-existant handle
    verts = [verts[0],int(23),verts[1]]
    try:
        mb.tag_set_data(int_tag,verts,int_data)
    except:
        pass
    else:
        print("Shouldn't be here. Test fails.")
        raise(AssertionError)

def test_unordered_tagging():

    mb = core.Core()
    # 1-D iterable
    coords = [0.,0.,0.,1.,0.,0.,1.,1.,1.]
    verts = mb.create_vertices(coords)
    CHECK_EQ(len(verts),3)
    # 2-D iterable w/ len 3 entries
    coords = [[0.,0.,0.],[1.,0.,0.],[1.,1.,1.]]
    verts = mb.create_vertices(coords)
    CHECK_EQ(len(verts),3)

    #create a tag
    int_tag = mb.tag_get_handle("IntTag",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)

    #try to set data with bad array (contains int)
    int_data = [1,2,3]

    # tag ordered vertices with value
    mb.tag_set_data(int_tag, verts, int_data)
        
    reordered_verts = [verts[1], verts[2], verts[0]]
    reordered_data =  [int_data[1], int_data[2], int_data[0]]

    # check that data array is correct for reordered vertex handles
    data = mb.tag_get_data(int_tag, reordered_verts)
    CHECK_EQ(data, reordered_data)
    
    
def test_vec_tags():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)

    int_vec_test_tag = mb.tag_get_handle("IntegerVecTestTag",3,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    dbl_vec_test_tag = mb.tag_get_handle("DoubleVecTestTag",3,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    opaque_vec_test_tag = mb.tag_get_handle("OPTag",10,types.MB_TYPE_OPAQUE,types.MB_TAG_DENSE,True)
    #should be able to successfully tag using a 2-D array
    int_vec_test_tag_values = np.array([[0,1,2],[3,4,5],[6,7,8]])
    dbl_vec_test_tag_values = np.array([[9.0,10.0,11.0],[12.0,13.0,14.0],[15.0,16.0,17.0]])
    opaque_vec_test_tag_values = np.array([["One"],["Two"],["Three"]])
    mb.tag_set_data(int_vec_test_tag,verts,int_vec_test_tag_values)
    mb.tag_set_data(dbl_vec_test_tag,verts,dbl_vec_test_tag_values)
    mb.tag_set_data(opaque_vec_test_tag,verts,opaque_vec_test_tag_values)

    #or a 1-D array
    int_vec_test_tag_values_flat = np.array([0,1,2,3,4,5,6,7,8])
    dbl_vec_test_tag_values_flat = np.array([9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0])
    opaque_vec_test_tag_values_flat = np.array(["One","Two","Three"])
    mb.tag_set_data(int_vec_test_tag,verts,int_vec_test_tag_values)
    mb.tag_set_data(dbl_vec_test_tag,verts,dbl_vec_test_tag_values)
    mb.tag_set_data(opaque_vec_test_tag,verts,opaque_vec_test_tag_values_flat)

    #these values should then be able to be retrieved as a 2-D array
    returned_int_test_tag_values = mb.tag_get_data(int_vec_test_tag,verts)
    CHECK_ITER_EQ(returned_int_test_tag_values,int_vec_test_tag_values)
    returned_dbl_test_tag_values = mb.tag_get_data(dbl_vec_test_tag,verts)
    CHECK_ITER_EQ(returned_dbl_test_tag_values,dbl_vec_test_tag_values)
    returned_opaque_test_tag_values = mb.tag_get_data(opaque_vec_test_tag,verts)
    CHECK_ITER_EQ(returned_opaque_test_tag_values,opaque_vec_test_tag_values)

    #or as a 1-D array
    returned_int_test_tag_values = mb.tag_get_data(int_vec_test_tag,verts,flat=True)
    CHECK_ITER_EQ(returned_int_test_tag_values,int_vec_test_tag_values_flat)
    returned_dbl_test_tag_values = mb.tag_get_data(dbl_vec_test_tag,verts,flat=True)
    CHECK_ITER_EQ(returned_dbl_test_tag_values,dbl_vec_test_tag_values_flat)
    returned_opaque_test_tag_values = mb.tag_get_data(opaque_vec_test_tag,verts,flat=True)
    CHECK_ITER_EQ(returned_opaque_test_tag_values,opaque_vec_test_tag_values_flat)

def test_create_element_iterable():
    mb = core.Core()
    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    CHECK_EQ(3,len(verts))
    #create tri using verts as a Range
    tris = mb.create_element(types.MBTRI,verts)
    #create another with the same vertices but in a list
    tri_verts = [verts[0], verts[1], verts[2]]
    tris = mb.create_element(types.MBTRI,verts)
    #make sure the right number of triangles is in the instance
    rs = mb.get_root_set()
    all_tris = mb.get_entities_by_type(rs,types.MBTRI)
    CHECK_EQ(len(all_tris),2)

def test_create_elements_iterable():
    mb = core.Core()
    #create some vertices for triangle 1
    coords1 = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts1 = mb.create_vertices(coords1)
    CHECK_EQ(len(verts1),3)
    #create some more vertices for triangle 2
    coords2 = np.array((1,2,3,4,5,6,1,1,1),dtype='float64')
    verts2 = mb.create_vertices(coords2)
    CHECK_EQ(len(verts2),3)
    tri_verts = [[verts1[0],verts1[1],verts1[2]],[verts2[0],verts2[1],verts2[2]]]
    tris = mb.create_elements(types.MBTRI,tri_verts)
    CHECK_EQ(len(tris),2)
    tri_verts = [verts1,verts2]
    tris = mb.create_elements(types.MBTRI,tri_verts)
    CHECK_EQ(len(tris),2)
    #make sure the right number of triangles is in the instance
    rs = mb.get_root_set()
    all_tris = mb.get_entities_by_type(rs,types.MBTRI)
    CHECK_EQ(len(all_tris),4)

def test_tag_root_set():
    # make sure that the root set can be tagged with data
    mb = core.Core()

    test_tag = mb.tag_get_handle("Test",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)

    root_set = mb.get_root_set()

    test_tag_data = np.array((4.,))
    mb.tag_set_data(test_tag, root_set, test_tag_data)

    data = mb.tag_get_data(test_tag, root_set, flat = True)
    CHECK_EQ(len(data), 1)
    CHECK_EQ(data[0],test_tag_data[0])

    mb.tag_delete_data(test_tag, root_set)

def test_entity_handle_tags():
    # make sure that the root set can be tagged with data
    mb = core.Core()

    eh_tag = mb.tag_get_handle("Test", 1, types.MB_TYPE_HANDLE, types.MB_TAG_SPARSE, True)
    dbl_tag = mb.tag_get_handle("Dbl", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
    
    meshset_a = mb.create_meshset()
    meshset_b = mb.create_meshset()

    # tag meshset a with meshset b
    mb.tag_set_data(eh_tag, meshset_a, meshset_b)
    # tag meshset b with a double value
    val = 16.0
    mb.tag_set_data(dbl_tag, meshset_b, val)

    eh = mb.tag_get_data(eh_tag, meshset_a)
    CHECK_EQ(eh, meshset_b)
    dbl_val = mb.tag_get_data(dbl_tag, eh)
    CHECK_EQ(dbl_val, val)

    eh = mb.tag_get_data(eh_tag, meshset_a, flat = True)[0]
    CHECK_EQ(eh, meshset_b)
    dbl_val = mb.tag_get_data(dbl_tag, eh, flat = True)[0]
    CHECK_EQ(dbl_val, val)
    
if __name__ == "__main__":
    tests = [test_load_mesh,
             test_write_mesh,
             test_delete_mesh,
             test_get_tag,
             test_integer_tag,
             test_double_tag,
             test_opaque_tag,
             test_tag_list,
             test_tag_delete,
             test_create_meshset,
             test_create_elements,
             test_tag_failures,
             test_adj,
             test_type_from_handle,
             test_meshsets,
             test_rs,
             test_get_coords,
             test_set_coords,
             test_get_ents_by_type,
             test_get_ents_by_tnt,
             test_get_entities_by_handle,
             test_get_entities_by_dimension,
             test_parent_child,
             test_remove_ents,
             test_iterables,
             test_vec_tags,
             test_create_element_iterable,
             test_create_elements_iterable,
             test_tag_root_set,
             test_entity_handle_tags]
    test_driver(tests)
