
from libcpp.string cimport string as std_string


cdef extern from "MBTagConventions.hpp":


    cdef std_string MATERIAL_SET_TAG_NAME
    cdef std_string DIRICHLET_SET_TAG_NAME
    cdef std_string NEUMANN_SET_TAG_NAME
    cdef std_string HAS_MID_NODES_TAG_NAME
    cdef std_string GEOM_DIMENSION_TAG_NAME
    cdef std_string MESH_TRANSFORM_TAG_NAME
    cdef std_string GLOBAL_ID_TAG_NAME
    cdef std_string CATEGORY_TAG_NAME
    cdef int        CATEGORY_TAG_SIZE
    cdef std_string NAME_TAG_NAME
    cdef int        NAME_TAG_SIZE

