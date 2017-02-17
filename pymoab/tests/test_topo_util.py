from pymoab import core
from pymoab import types
from pymoab import topo_util
import numpy as np


def test_get_bridge_adjacencies():
    mb = core.Core()

    coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
    verts = mb.create_vertices(coords)
    verts = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)

    mtu = topo_util.MeshTopoUtil(mb)

    rs = mb.get_root_set()
    tris = mb.get_entities_by_dimension(rs, 2)
    adj_verts = mtu.get_bridge_adjacencies(tris, 0, 0, 1)

    assert (adj_verts == verts).all()


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
