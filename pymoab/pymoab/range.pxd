"""Implements range functionality."""

from pymoab cimport moab

cdef class Range:

    cdef moab.Range * inst

