from pymoab import core
from pymoab import types
from pymoab import topo_util
from driver import test_driver
from driver import CHECK_EQ, CHECK_NOT_EQ, CHECK_ITER_EQ
import numpy as np


def test_get_bridge_adjacencies():
    mb = core.Core()

    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    CHECK_EQ(len(verts), 3)
    conn = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI, conn)

    mtu = topo_util.MeshTopoUtil(mb)

    rs = mb.get_root_set()
    tris = mb.get_entities_by_dimension(rs, 2)
    CHECK_EQ(len(tris), 1)
    adj_verts = mtu.get_bridge_adjacencies(tris, 0, 0, 1)

    CHECK_ITER_EQ(adj_verts,verts)


def test_get_average_position():
    mb = core.Core()

    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    verts = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)

    mtu = topo_util.MeshTopoUtil(mb)

    rs = mb.get_root_set()
    tris = mb.get_entities_by_dimension(rs, 2)
    avg_pos = mtu.get_average_position(tris)

    np.testing.assert_almost_equal(avg_pos, coords.reshape(3,3).mean(axis=0))


def test_construct_aentities():
    mb = core.Core()

    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    verts_array = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    mb.create_elements(types.MBTRI,verts_array)

    mtu = topo_util.MeshTopoUtil(mb)

    edges_before = mb.get_entities_by_dimension(0, 1)

    mtu.construct_aentities(verts)

    edges_after = mb.get_entities_by_dimension(0, 1)

    CHECK_NOT_EQ(edges_before, edges_after)
    CHECK_EQ(len(edges_after), 3)

if __name__ == "__main__":
    tests = [test_get_bridge_adjacencies,
             test_get_average_position,
             test_construct_aentities]
    test_driver(tests)

