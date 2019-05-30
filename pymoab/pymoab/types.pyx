"""Python wrappers for MOAB Types."""

from pymoab cimport moab
from pymoab cimport tag_conventions

cimport numpy as np
import numpy as np

_eh_py_type = np.uint64

cdef class MOABErrorCode:

    cdef readonly moab.ErrorCode error_value

    cdef readonly dict err_strings
    def __cinit__(self, value = 0):
        if isinstance(value, MOABErrorCode):
            self.error_value = value.error_value
        else:
            self.error_value = <moab.ErrorCode> value

        self.err_strings = { moab.MB_SUCCESS : "MB_SUCCESS",
                             moab.MB_INDEX_OUT_OF_RANGE : "MB_INDEX_OUT_OF_RANGE",
                             moab.MB_TYPE_OUT_OF_RANGE : "MB_TYPE_OUT_OF_RANGE",
                             moab.MB_MEMORY_ALLOCATION_FAILED : "MB_MEMORY_ALLOCATION_FAILED",
                             moab.MB_ENTITY_NOT_FOUND : "MB_ENTITY_NOT_FOUND",
                             moab.MB_MULTIPLE_ENTITIES_FOUND : "MB_MULTIPLE_ENTITIES_FOUND",
                             moab.MB_TAG_NOT_FOUND : "MB_TAG_NOT_FOUND",
                             moab.MB_FILE_DOES_NOT_EXIST : "MB_FILE_DOES_NOT_EXIST",
                             moab.MB_FILE_WRITE_ERROR : "MB_FILE_WRITE_ERROR",
                             moab.MB_NOT_IMPLEMENTED : "MB_NOT_IMPLEMENTED",
                             moab.MB_ALREADY_ALLOCATED : "MB_ALREADY_ALLOCATED",
                             moab.MB_VARIABLE_DATA_LENGTH : "MB_VARIABLE_DATA_LENGTH",
                             moab.MB_INVALID_SIZE : "MB_INVALID_SIZE",
                             moab.MB_UNSUPPORTED_OPERATION : "MB_UNSUPPORTED_OPERATION",
                             moab.MB_UNHANDLED_OPTION : "MB_UNHANDLED_OPTION",
                             moab.MB_STRUCTURED_MESH : "MB_STRUCTURED_MESH",
                             moab.MB_FAILURE : "MB_FAILURE" }
        
    def __richcmp__(self, other, op):
        if op == 2:
            if isinstance(other, MOABErrorCode):
                return self.error_value == other.error_value
            elif type(other) == int:
                return self.error_value == other
        else:
            return NotImplemented
        
    def __hash__(self):
        return self.error_value

    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        return "MOAB ErrorCode: "+self.err_strings[self.error_value]
        
    
# Error codes
MB_SUCCESS = MOABErrorCode(moab.MB_SUCCESS)
MB_INDEX_OUT_OF_RANGE = MOABErrorCode(moab.MB_INDEX_OUT_OF_RANGE)
MB_TYPE_OUT_OF_RANGE = MOABErrorCode(moab.MB_TYPE_OUT_OF_RANGE) 
MB_MEMORY_ALLOCATION_FAILED = MOABErrorCode(moab.MB_MEMORY_ALLOCATION_FAILED)     
MB_ENTITY_NOT_FOUND = MOABErrorCode(moab.MB_ENTITY_NOT_FOUND)
MB_MULTIPLE_ENTITIES_FOUND = MOABErrorCode(moab.MB_MULTIPLE_ENTITIES_FOUND)  
MB_TAG_NOT_FOUND = MOABErrorCode(moab.MB_TAG_NOT_FOUND)
MB_FILE_DOES_NOT_EXIST = MOABErrorCode(moab.MB_FILE_DOES_NOT_EXIST)  
MB_FILE_WRITE_ERROR = MOABErrorCode(moab.MB_FILE_WRITE_ERROR)
MB_NOT_IMPLEMENTED = MOABErrorCode(moab.MB_NOT_IMPLEMENTED)
MB_ALREADY_ALLOCATED = MOABErrorCode(moab.MB_ALREADY_ALLOCATED)    
MB_VARIABLE_DATA_LENGTH = MOABErrorCode(moab.MB_VARIABLE_DATA_LENGTH)     
MB_INVALID_SIZE = MOABErrorCode(moab.MB_INVALID_SIZE)
MB_UNSUPPORTED_OPERATION = MOABErrorCode(moab.MB_UNSUPPORTED_OPERATION)    
MB_UNHANDLED_OPTION = MOABErrorCode(moab.MB_UNHANDLED_OPTION)
MB_STRUCTURED_MESH = MOABErrorCode(moab.MB_STRUCTURED_MESH)
MB_FAILURE = MOABErrorCode(moab.MB_FAILURE)


cdef dict _ERROR_MSGS = {
    MB_INDEX_OUT_OF_RANGE: (IndexError, 'MOAB index out of range'),
    MB_TYPE_OUT_OF_RANGE: (TypeError, 'Incorrect MOAB type, out of range'),
    MB_MEMORY_ALLOCATION_FAILED: (MemoryError, 'MOAB memory allocation'),
    MB_ENTITY_NOT_FOUND: (RuntimeError, 'Entity not found'),
    MB_MULTIPLE_ENTITIES_FOUND: (RuntimeError, 'Multiple entities found'),
    MB_TAG_NOT_FOUND: (RuntimeError, 'Tag not found'),
    MB_FILE_DOES_NOT_EXIST: (IOError, 'File not found'),
    MB_FILE_WRITE_ERROR: (IOError, 'File write error'), 
    MB_NOT_IMPLEMENTED: (NotImplementedError, '[MOAB]'),
    MB_ALREADY_ALLOCATED: (MemoryError, 'already allocated'),
    MB_VARIABLE_DATA_LENGTH: (TypeError, 'variable length data'),
    MB_INVALID_SIZE: (ValueError, 'invalid size'),
    MB_UNSUPPORTED_OPERATION: (RuntimeError, 'unsupported operation'),
    MB_UNHANDLED_OPTION: (RuntimeError, 'unhandled option'),
    MB_STRUCTURED_MESH: (RuntimeError, 'structured mesh'),
    MB_FAILURE: (RuntimeError, '[MOAB] failure'),
    }

def check_error(err, tuple exceptions = (), **kwargs):
    """Checks error status code and raises error if needed."""
    for exception in exceptions:
        if exception == err:
            return
    if err == MB_SUCCESS:
        return
    errtype, msg = _ERROR_MSGS[err]
    if len(kwargs) > 0:
        msg += ': '
        msg += ', '.join(sorted(['{0}={1!r}'.format(k, v) for k, v in kwargs.items()]))
    raise errtype(MOABErrorCode(err))

# Data Types
MB_TYPE_OPAQUE = moab.MB_TYPE_OPAQUE
MB_TYPE_INTEGER = moab.MB_TYPE_INTEGER
MB_TYPE_DOUBLE = moab.MB_TYPE_DOUBLE
MB_TYPE_BIT = moab.MB_TYPE_BIT
MB_TYPE_HANDLE = moab.MB_TYPE_HANDLE
MB_MAX_DATA_TYPE = moab.MB_MAX_DATA_TYPE


_TAG_TYPE_STRS = {
    MB_TYPE_OPAQUE : "MB_TYPE_OPAQUE",
    MB_TYPE_INTEGER : "MB_TYPE_INTEGER",
    MB_TYPE_DOUBLE : "MB_TYPE_DOUBLE",
    MB_TYPE_BIT : "MB_TYPE_BIT",
    MB_TYPE_HANDLE : "MB_TYPE_HANDLE",
    MB_MAX_DATA_TYPE : "MB_MAX_DATA_TYPE"
}

_DTYPE_CONV = {
    MB_TYPE_OPAQUE: 'S',
    MB_TYPE_INTEGER: 'int32',
    MB_TYPE_DOUBLE: 'float64',
    MB_TYPE_BIT: 'bool',
    MB_TYPE_HANDLE: 'uint64',
    MB_MAX_DATA_TYPE: 'uint64'
}

_VALID_DTYPES= {
    MB_TYPE_OPAQUE: frozenset(['S','U','O','object']),
    MB_TYPE_INTEGER: frozenset(['i','int8','int16','int32','int64','O','object']),
    MB_TYPE_DOUBLE: frozenset(['float64','float','f8','f', 'O','object',]),
    MB_TYPE_BIT: frozenset(['f','int8','int16','int32','int64','S1','bool','O','object']),
    MB_TYPE_HANDLE: frozenset(['uint64','O','object']),
    MB_MAX_DATA_TYPE: frozenset(['uint64','O','object'])
}

_VALID_NATIVE_TYPES = {
    int: MB_TYPE_INTEGER,
    float: MB_TYPE_DOUBLE,
    str: MB_TYPE_OPAQUE,
    object: MB_TYPE_OPAQUE,
    np.uint64 : MB_TYPE_HANDLE
}

def pymoab_data_type(input_type):
    """
    Attempts to find a PyMOAB datatype given a Python native type or 
    NumPy dtype
    """

    # check native types
    try:
        t = _VALID_NATIVE_TYPES[input_type]
        return t
    except KeyError:
        pass

    # check valid dtypes
    for k in _VALID_DTYPES.keys():
        if input_type in _VALID_DTYPES[k]:
            return k

    raise ValueError("Could not determine the PyMOAB data type.")
    
def _convert_array(iterable, accepted_types, return_dtype):
    err_msg = "Incorrect datatype found in array: {}"
    #if this is already an array of the correct type, avoid the loop
    if isinstance(iterable, np.ndarray) and iterable.dtype == return_dtype:
        return  iterable
    #if not, each entry in the iterable should be verified
    for entry in iterable:
        assert (isinstance(entry, accepted_types)), err_msg.format(type(entry))
        #if this is true, then create an array from the iterable
    return np.fromiter(iterable, return_dtype)

def _eh_array(iterable):
    err_msg = """
               Invalid EntityHandle type is being used.  Please ensure all
               EntityHandles are either Python type long or NumPy dtype uint64.
              """
    # get dtype for EntityHandles
    EH_DTYPE = _DTYPE_CONV[MB_TYPE_HANDLE]
    # try to convert array
    try:
        arr = _convert_array(iterable, (_eh_py_type,), EH_DTYPE)
    except:
        raise ValueError(err_msg)
    # return array if successful
    return arr
    
def np_tag_type(type):
    return _DTYPE_CONV[type]

def validate_type(tag_type,tag_length,tag_data):

    #ensure this type is supported
    assert tag_type in _DTYPE_CONV.keys(), "Invalid Tag Type"
    #and that it is the correct shape
    assert tag_data.ndim == 0 or tag_data.ndim == 1, "Not a flat array"

    
    if MB_TYPE_OPAQUE == tag_type:
        #so long as the array is some kind of string type, we're happy
        is_valid = tag_data.dtype.char in _VALID_DTYPES[tag_type]
        final_type = _DTYPE_CONV[tag_type]+str(tag_length)
    else:
        is_valid = str(tag_data.dtype) in _VALID_DTYPES[tag_type]
        final_type = _DTYPE_CONV[tag_type]

    assert is_valid, "Data is invalid for Tag. Please verify data length and type."
    tag_data = np.asarray(tag_data, dtype=final_type)
    return tag_data

# Entity types
MBVERTEX = moab.MBVERTEX
MBEDGE = moab.MBEDGE
MBTRI = moab.MBTRI
MBQUAD = moab.MBQUAD
MBPOLYGON = moab.MBPOLYGON
MBTET = moab.MBTET
MBPYRAMID = moab.MBPYRAMID
MBPRISM = moab.MBPRISM
MBKNIFE = moab.MBKNIFE
MBHEX = moab.MBHEX
MBPOLYHEDRON = moab.MBPOLYHEDRON
MBENTITYSET = moab.MBENTITYSET
MBMAXTYPE = moab.MBMAXTYPE

# Tag Types
MB_TAG_BIT  = moab.MB_TAG_BIT   
MB_TAG_SPARSE = moab.MB_TAG_SPARSE
MB_TAG_DENSE = moab.MB_TAG_DENSE
MB_TAG_MESH = moab.MB_TAG_MESH  
MB_TAG_BYTES = moab.MB_TAG_BYTES 
MB_TAG_VARLEN = moab.MB_TAG_VARLEN
MB_TAG_CREAT = moab.MB_TAG_CREAT 
MB_TAG_EXCL = moab.MB_TAG_EXCL  
MB_TAG_STORE = moab.MB_TAG_STORE 
MB_TAG_ANY = moab.MB_TAG_ANY   
MB_TAG_NOOPQ = moab.MB_TAG_NOOPQ 
MB_TAG_DFTOK = moab.MB_TAG_DFTOK 

# Query selection types
INTERSECT = 0
UNION = 1

MESHSET_TRACK_OWNER = moab.MESHSET_TRACK_OWNER
MESHSET_SET = moab.MESHSET_SET
MESHSET_ORDERED = moab.MESHSET_ORDERED

MATERIAL_SET_TAG_NAME   = str(tag_conventions.MATERIAL_SET_TAG_NAME.decode())
DIRICHLET_SET_TAG_NAME  = str(tag_conventions.DIRICHLET_SET_TAG_NAME.decode())
NEUMANN_SET_TAG_NAME    = str(tag_conventions.NEUMANN_SET_TAG_NAME.decode())
HAS_MID_NODES_TAG_NAME  = str(tag_conventions.HAS_MID_NODES_TAG_NAME.decode())
GEOM_DIMENSION_TAG_NAME = str(tag_conventions.GEOM_DIMENSION_TAG_NAME.decode())
MESH_TRANSFORM_TAG_NAME = str(tag_conventions.MESH_TRANSFORM_TAG_NAME.decode())
GLOBAL_ID_TAG_NAME      = str(tag_conventions.GLOBAL_ID_TAG_NAME.decode())
CATEGORY_TAG_NAME       = str(tag_conventions.CATEGORY_TAG_NAME.decode())
CATEGORY_TAG_SIZE       = tag_conventions.CATEGORY_TAG_SIZE
NAME_TAG_NAME           = str(tag_conventions.NAME_TAG_NAME.decode())
NAME_TAG_SIZE           = tag_conventions.NAME_TAG_SIZE
