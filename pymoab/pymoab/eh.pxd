"""Header file for MOAB EntityHandle type"""
cimport numpy as np

cdef extern from "moab/EntityHandle.hpp" namespace "moab":

    ctypedef np.uint64_t EntityID
    ctypedef np.uint64_t EntityHandle
