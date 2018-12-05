from ctypes import *

from pymoab import core, types, skinner
from driver import test_driver, CHECK_EQ, CHECK_NOT_EQ
import numpy as np

def test_get_geometric_skin():
    mb = core.Core()

    coords = np.array([[0,0,0], [1,0,0], [2,0,0],
                       [0,1,0], [1,1,0], [2,1,0],
                       [0,2,0], [1,2,0], [2,2,0]],dtype='float64')
    verts = mb.create_vertices(coords)
    verts = np.array((verts[0:8]),dtype='uint64')
    connect = np.array([[1, 2, 5, 4],
                        [2, 3, 6, 5],
                        [4, 5, 8, 7],
                        [5, 6, 9, 8]],dtype='uint64')
    quads = mb.create_elements(types.MBQUAD,connect)

    # create a geometric set for the quads
    surf_set = mb.create_meshset()

    geom_tag = mb.tag_get_handle(types.GEOM_DIMENSION_TAG_NAME,
                                 1,
                                 types.MB_TYPE_INTEGER,
                                 types.MB_TAG_SPARSE,
                                 create_if_missing = True)

    mb.tag_set_data(geom_tag, surf_set, 2)

    # place 2-D entities in this set
    mb.add_entities(surf_set, verts)
    mb.add_entities(surf_set, quads)

    # create a dummy volume set
    vol_set = mb.create_meshset()
    mb.tag_set_data(geom_tag, vol_set, 3)

    # set surface to volume parent-child relationship
    mb.add_parent_meshset(surf_set, vol_set)

    mskn = skinner.Skinner(mb)

    rs = mb.get_root_set()
    skin = mskn.find_geometric_skin(rs)

    CHECK_EQ(skin.num_of_type(types.MBVERTEX), 8)
    CHECK_EQ(skin.num_of_type(types.MBQUAD), 4)

def test_get_skin():
    mb = core.Core()
    rs = mb.get_root_set()

    coords = np.array([[0,0,0], [1,0,0], [2,0,0],
                       [0,1,0], [1,1,0], [2,1,0],
                       [0,2,0], [1,2,0], [2,2,0]],dtype='float64')
    verts = mb.create_vertices(coords)
    verts = np.array((verts[0:8]),dtype='uint64')
    connect = np.array([[1, 2, 5, 4],
                        [2, 3, 6, 5],
                        [4, 5, 8, 7],
                        [5, 6, 9, 8]],dtype='uint64')
    quads = mb.create_elements(types.MBQUAD,connect)

    nelems = mb.get_entities_by_dimension(rs, 2)
    CHECK_EQ(len(nelems), 4)

    mb.write_file('quads_test.vtk')

    mskn = skinner.Skinner(mb)

    skin_verts = mskn.find_skin(rs, quads, True, False)
    CHECK_EQ(len(skin_verts), 8)

    skin_edges = mskn.find_skin(rs, quads, False, False)
    CHECK_EQ(len(skin_edges), 8)

if __name__ == "__main__":
    tests = [
             test_get_geometric_skin,
             test_get_skin]
    test_driver(tests)
