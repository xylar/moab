"""MOAB Tag Class"""

from pymoab cimport moab

cdef class Tag:
    cdef moab.TagInfo * inst
    cdef moab.Tag * ptr
 
cdef class _tagArray:
    cdef moab.TagInfo ** inst
    cdef moab.Tag * ptr
