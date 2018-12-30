"""MOAB Tag Class"""

from pymoab cimport moab
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc,free
from cython.operator cimport dereference as deref
from .types import _TAG_TYPE_STRS, _DTYPE_CONV

cdef class Tag(object):
    def __cinit__(self):
        self.inst = <moab.TagInfo*> malloc(sizeof(moab.TagInfo))
        self.ptr = &self.inst

    def __del__(self):
        free(self.inst)

    def get_length(self):
        """
        Returns the length (number of entries) for this Tag as an integer.
        """
        t = self.inst.get_data_type()
        type_byte_size = self.inst.size_from_data_type(t)
        total_byte_size = self.inst.get_size()
        return int(total_byte_size/type_byte_size)

    def get_data_type(self):
        """
        Returns the Tag's data type.
        """
        return self.get_type()

    def get_type(self):
        """
        Returns the Tag's data type.
        """
        return self.inst.get_data_type()

    def get_dtype(self):
        """
        Returns the Tag's numpy data type.
        """
        return _DTYPE_CONV[self.inst.get_data_type()]

    def get_name(self):
        """
        Returns the name of this Tag.
        """
        return str(self.inst.get_name().c_str().decode())

    def get_default_value(self):
        """
        Returns the default value of the tag.
        If the tag does not have a default value, None is returned.
        """
        tag_size = self.get_length()
        dtype = self.get_dtype()
        cdef np.ndarray arr = np.empty((tag_size,), dtype = dtype)

        cdef const void* data_ptr = self.inst.get_default_value()
        arr.data = <char *> data_ptr

        if tag_size == 1:
            return arr[0]
        else:
            return arr

    def __str__(self):
        outstr = "Name: " + self.get_name()
        outstr += ", Type: " + _TAG_TYPE_STRS[self.get_data_type()]
        outstr += ", Length: " + str(self.get_length())
        return outstr

    def __repr__(self):
        return self.__str__()

cdef class _tagArray(object):
    def __cinit__(self, tags = None):
        if tags == None:
            return

        cdef Tag t
        cdef int num_tags
        cdef int i = 0
        if isinstance(tags,Tag):
            self.inst =  <moab.TagInfo**> malloc(sizeof(moab.TagInfo*))
            t = <Tag> tags
            self.inst[0] = t.inst
        else:
            num_tags = len(tags)
            self.inst = <moab.TagInfo**> malloc(num_tags*sizeof(moab.TagInfo*))
            for i in range(0,num_tags):
                t = <Tag> tags[i]
                self.inst[i] = t.inst
        self.ptr = self.inst

    def __del__(self):
        free(self.inst)
