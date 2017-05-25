""" Core Cython Header """

from pymoab cimport moab

cdef class Core:
    cdef moab.Core* inst
