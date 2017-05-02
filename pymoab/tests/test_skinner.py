
from ctypes import *
lib1 = cdll.LoadLibrary('/usr/software/moab/dev/pymoab/lib/libMOAB.so')

from pymoab import core
from pymoab import types
from pymoab import skinner
from driver import test_driver, CHECK_EQ, CHECK_NOT_EQ
import numpy as np


def test_get_geometric_skin():
    mb = core.Core()

    coords = np.array((0,0,0, 1,0,0, 2,0,0,
                       0,1,0, 1,1,0, 2,1,0,
                       0,2,0, 1,2,0, 2,2,0),dtype='float64')
    verts = mb.create_vertices(coords)
    verts = np.array((verts[0:8]),dtype='uint64')
    connect = np.array([[0, 1, 4, 3],
                        [1, 2, 5, 4],
                        [3, 4, 7, 6],
                        [4, 5, 8, 7]],dtype='uint64')
    quads = mb.create_elements(types.MBQUAD,connect)

    mskn = skinner.Skinner(mb)

    rs = mb.get_root_set()
    skin_verts = mskn.find_geometric_skin(rs)

    CHECK_EQ(len(skin_verts), 8)


def test_get_skin():
    mb = core.Core()
    rs = mb.get_root_set()

    coords = np.array((0,0,0, 1,0,0, 2,0,0,
                       0,1,0, 1,1,0, 2,1,0,
                       0,2,0, 1,2,0, 2,2,0),dtype='float64')
    verts = mb.create_vertices(coords)
    verts = np.array((verts[0:8]),dtype='uint64')
    connect = np.array([[0, 1, 4, 3],
                        [1, 2, 5, 4],
                        [3, 4, 7, 6],
                        [4, 5, 8, 7]],dtype='uint64')
    quads = mb.create_elements(types.MBQUAD,connect)

    nelems = mb.get_entities_by_dimension(rs, 2)
    CHECK_EQ(len(nelems), 4)

    mskn = skinner.Skinner(mb)

    skin_verts = mskn.find_skin(rs, quads, False, False)
    CHECK_EQ(len(skin_verts), 8)


if __name__ == "__main__":
    tests = [
#             test_get_geometric_skin,
             test_get_skin]
    test_driver(tests)

