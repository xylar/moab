"""Python wrappers for MOAB Types."""

from pymoab cimport moab
cimport numpy as np
import numpy as np

cdef class MOABErrorCode:

    cdef readonly moab.ErrorCode error_value
    
    def __cinit__(self, value = 0):
        self.error_value = <moab.ErrorCode> value

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
    raise errtype(msg)

# Data Types
MB_TYPE_OPAQUE = moab.MB_TYPE_OPAQUE
MB_TYPE_INTEGER = moab.MB_TYPE_INTEGER
MB_TYPE_DOUBLE = moab.MB_TYPE_DOUBLE
MB_TYPE_BIT = moab.MB_TYPE_BIT
MB_TYPE_HANDLE = moab.MB_TYPE_HANDLE
MB_MAX_DATA_TYPE = moab.MB_MAX_DATA_TYPE

_DTYPE_CONV = {
    MB_TYPE_OPAQUE: 'S',
    MB_TYPE_INTEGER: 'int32',
    MB_TYPE_DOUBLE: 'float64',
    MB_TYPE_BIT: 'bool',
    MB_TYPE_HANDLE: 'uint64',
    MB_MAX_DATA_TYPE: 'uint64'
}

_VALID_DTYPES= {
    MB_TYPE_OPAQUE: frozenset(['S']),
    MB_TYPE_INTEGER: frozenset(['int8','int16','int32','int64']),
    MB_TYPE_DOUBLE: frozenset(['float64']),
    MB_TYPE_BIT: frozenset(['int8','int16','int32','int64','S1','bool']),
    MB_TYPE_HANDLE: frozenset(['uint64']),
    MB_MAX_DATA_TYPE: frozenset(['uint64'])
}

_VALID_DTYPES= {
    MB_TYPE_OPAQUE: frozenset(['S','O']),
    MB_TYPE_INTEGER: frozenset(['int8','int16','int32','int64','O','object']),
    MB_TYPE_DOUBLE: frozenset(['float64','float','O','object']),
    MB_TYPE_BIT: frozenset(['int8','int16','int32','int64','S1','bool','O','object']),
    MB_TYPE_HANDLE: frozenset(['uint64','O','object']),
    MB_MAX_DATA_TYPE: frozenset(['uint64','O','object'])
}

def _convert_array(iterable, accepted_dtypes, return_dtype):
    err_msg = "Incorrect type in EntityHandle Array"
    #if this is already an array of the correct type, avoid the loop
    if isinstance(iterable, np.ndarray) and iterable.dtype == return_dtype:
        return  iterable
    #if not, each entry in the iterable should be verified
    for entry in iterable:
        assert (type(entry) in accepted_dtypes), err_msg
    #if this is true, then create an array from the iterable
    return np.fromiter(iterable, return_dtype)

def _eh_array(iterable):
    EH_DTYPE = _DTYPE_CONV[MB_TYPE_HANDLE]
    return _convert_array(iterable, [long, EH_DTYPE], EH_DTYPE)

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
