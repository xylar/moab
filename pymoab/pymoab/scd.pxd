"""MOAB Structured Mesh Interface"""

from pymoab cimport moab

cdef class ScdParData:
    cdef moab.ScdParData *inst

cdef class ScdInterface:

    cdef moab.ScdInterface * inst
    cdef moab.Interface * interface

cdef class ScdBox:
    cdef moab.ScdBox* inst
