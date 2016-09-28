# Test that numpy works in Cython:
from numpy cimport ndarray
cimport numpy as np
import numpy as np
import glob
import os

# Test that libcpp module is present:
from libcpp.vector cimport vector
from libcpp.string cimport string as std_string

# Test that templates work:
cdef vector[int] array2vector_int(a):
    cdef vector[int] v
    for i in range(len(a)):
        v.push_back(a[i])
    return v

# Test that bool support works
from libcpp cimport bool

# Test that C++ support works
cdef extern from "TagInfo.hpp" namespace "moab":

    cdef cppclass TagInfo:
        TagInfo()

cdef extern from 'moab/Types.hpp' namespace "moab":

    ctypedef enum EntitySetProperty:
        MESHSET_TRACK_OWNER = 0x01
        MESHSET_SET = 0x02
        MESHSET_ORDERED = 0x03

    ctypedef TagInfo* Tag

cdef extern from "moab/Core.hpp" namespace "moab":

    cdef cppclass Core:
        # Constructors
        Core()


cdef class mbcore(object):

    cdef Core * inst

    def __cinit__(self):
        self.inst = new Core()
