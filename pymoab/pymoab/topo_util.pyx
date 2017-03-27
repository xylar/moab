""" MOAB MeshTopoUtil """
from cython.operator cimport dereference as deref

cimport numpy as np
import numpy as np

from .rng cimport Range
from .core cimport Core
from .types import check_error
from . import types


cdef class MeshTopoUtil(object):

    def __cinit__(self, Core c):
        self.interface  = <moab.Interface*> c.inst
        self.inst = new moab.MeshTopoUtil(self.interface)

    def __del__(self):
        del self.inst

    def get_bridge_adjacencies(self,
                               from_ent,
                               int bridge_dim,
                               int to_dim,
                               int num_layers = 1,
                               exceptions = ()):
        cdef moab.ErrorCode err
        cdef moab.EntityHandle ms_handle
        cdef Range r
        cdef Range adjs = Range()

        if isinstance(from_ent, Range):
            r = from_ent
            err = self.inst.get_bridge_adjacencies(deref(r.inst), bridge_dim, to_dim, deref(adjs.inst), num_layers)
        else:
            ms_handle = from_ent
            err = self.inst.get_bridge_adjacencies(ms_handle, bridge_dim, to_dim, deref(adjs.inst))

        check_error(err, exceptions)

        return adjs

    def get_average_position(self,
                             entity_handles,
                             exceptions = ()):
        cdef moab.ErrorCode err
        cdef moab.EntityHandle ms_handle
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef np.ndarray avg_position

        if isinstance(entity_handles, Range):
            r = entity_handles
            avg_position = np.empty((3,),dtype='float64')
            err = self.inst.get_average_position(deref(r.inst), <double*> avg_position.data)
            check_error(err, exceptions)
        elif isinstance(entity_handles,np.ndarray):
            assert entity_handles.dtype == 'uint64'
            arr = entity_handles
            avg_position = np.empty((3,),dtype='float64')
            err = self.inst.get_average_position(<unsigned long*> arr.data, len(entity_handles), <double*> avg_position.data)
            check_error(err, exceptions)
        else:
            check_error(types.MB_FAILURE)

        return avg_position

    def construct_aentities(self,
                            vertices,
                            exceptions = ()):
        """
        Generate all the AEntities bounding the vertices.

        Example
        -------
        root_set = mb.get_root_set()
        all_verts = mb.get_entities_by_dimension(root_set, 0)
        mtu.construct_aentities(all_verts)

        Parameters
        ----------
        vertices : Range or iterable of EntityHandles
            Vertices that will be used to generate the bounding AEntities
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)
        """
        cdef moab.ErrorCode err
        cdef Range r

        r = vertices
        err = self.inst.construct_aentities(deref(r.inst))

        check_error(err, exceptions)
