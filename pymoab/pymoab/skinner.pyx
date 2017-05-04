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
        Find the geometric skin of the domain
        This method requires that GEOM_DIMENSION tag is available on entities
        to specify the volumes that they belong to. Internally, the queries 
        will use CGM.

        Returns entities that form the geometric skin.

        Example
        -------
        ms_handle # the mesh set corresponding for which boundary is requested
        ms_handle = mb.get_rootset()
        skin_ents = find_geometric_skin(ms_handle)

        Parameters
        ----------
        ms_handle : EntityHandle
            EntityHandle of the Meshset that is being queried
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        entities : Range or iterable of EntityHandles
            Entities that are on the geometric boundary of the queried set

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if the meshset EntityHandle is not of the correct type
        """
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        err = self.inst.find_geometric_skin(ms_handle, deref(ents.inst))
        check_error(err, exceptions)

        return ents


    def find_skin(self, const moab.EntityHandle ms_handle, entities, bint get_vertices=False, bint is_scd=False, exceptions = ()):
        """
        Find the entities that are on the topological boundary of the specified meshset.
        User can query for either vertices or get the d-1 dimensional entities on the
        mesh boundary.

        If the meshset belongs to a SCD set, the information is forwarded to the 
        underlying call appropriately.

        Example
        -------
        ms_handle # the mesh set corresponding for which boundary is requested
        ms_handle = mb.get_rootset()
        quads = mb.get_entities_by_dimension(ms_handle, 2)
        skin_verts = mskn.find_skin(rs, quads, True, False)
        skin_edges = mskn.find_skin(rs, quads, False, False)

        Parameters
        ----------
        ms_handle : EntityHandle
            EntityHandle of the Meshset that is being queried
        get_vertices : Boolean
            This flag indicates whether the user wants the 
            boundary vertices or d-1 elements
        is_scd : Boolean
            This flag indicates whether the underlying mesh object
            is using the SCD interface
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        entities : Range or iterable of EntityHandles
            Entities that are on the mesh boundary of the queried set.
            The entities may contain either the skin vertices or the
            skin elements depending on the \p get_vertices flag value.

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if the meshset EntityHandle is not of the correct type
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

