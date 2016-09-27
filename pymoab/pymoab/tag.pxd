"""MOAB Tag Class"""

from pymoab cimport moab

cdef class Tag:
    cdef moab.TagInfo * inst