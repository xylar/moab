""" MOAB MeshTopoUtil """
from cython.operator cimport dereference as deref

cimport numpy as np
import numpy as np

from pymoab cimport moab
from .rng cimport Range
from .core cimport Core
from .types import check_error, np_tag_type, validate_type, _convert_array, _eh_array
from . import types

from libcpp.vector cimport vector
from libcpp.string cimport string as std_string
from libc.stdlib cimport malloc

cdef void* null = NULL

cdef class Skinner(object):

    def __cinit__(self, Core c):
        """ Constructor """
        self.interface  = <moab.Interface*> c.inst
        self.inst = new moab.Skinner(self.interface)

    def __del__(self):
        """ Destructor """
        del self.inst

    def find_geometric_skin(self, moab.EntityHandle ms_handle, exceptions = ()):
        """
        Add entities to the specified meshset. Entities can be provided either
        as a pymoab.rng.Range o bject or as an iterable of EntityHandles.

        If meshset has MESHSET_TRACK_OWNER option set, adjacencies are also
        added to entities in entities.

        Example
        -------
        vertices # a iterable of MOAB vertex handles
        new_meshset = mb.create_meshset()
        mb.add_entities(new_meshset, entities)

        Parameters
        ----------
        ms_handle : EntityHandle
            EntityHandle of the Meshset the entities will be added to
        entities : Range or iterable of EntityHandles
            Entities that to add to the Meshset
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        None

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        err = self.inst.find_geometric_skin(ms_handle, deref(ents.inst))
        check_error(err, exceptions)

        return ents


    def find_skin(self, const moab.EntityHandle ms_handle, entities, bint get_vertices=False, bint is_scd=False, exceptions = ()):
        """
        Add entities to the specified meshset. Entities can be provided either
        as a pymoab.rng.Range o bject or as an iterable of EntityHandles.

        If meshset has MESHSET_TRACK_OWNER option set, adjacencies are also
        added to entities in entities.

        Example
        -------
        vertices # a iterable of MOAB vertex handles
        new_meshset = mb.create_meshset()
        mb.add_entities(new_meshset, entities)

        Parameters
        ----------
        ms_handle : EntityHandle
            EntityHandle of the Meshset the entities will be added to
        entities : Range or iterable of EntityHandles
            Entities that to add to the Meshset
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        None

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """
        cdef moab.ErrorCode err
        cdef Range r
        cdef Range outhandles = Range()
        cdef bint create_vert_elem_adjs=True, create_skin_elements=True

        if isinstance(entities, Range):
            r = entities
        else:
            r = Range(entities)
        err = self.inst.find_skin(ms_handle, deref(r.inst), get_vertices, deref(outhandles.inst), NULL, create_vert_elem_adjs, create_skin_elements, is_scd)
        check_error(err, exceptions)

        return outhandles

