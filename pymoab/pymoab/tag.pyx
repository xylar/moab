"""MOAB Tag Class"""

from pymoab cimport moab 
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc,free
from cython.operator cimport dereference as deref

cdef class Tag(object): 
    def __cinit__(self):
        self.inst = <moab.TagInfo*> malloc(sizeof(moab.TagInfo))
        self.ptr = &self.inst
        
    def __del__(self):
        free(self.inst)

    def get_length(self):
        t = self.inst.get_data_type()
        type_byte_size = self.inst.size_from_data_type(t)
        total_byte_size = self.inst.get_size()
        return total_byte_size/type_byte_size

    def get_data_type(self): 
        return self.inst.get_data_type()

    def get_name(self):
        return self.inst.get_name()
    
cdef class TagArray(object): 
    def __cinit__(self, tags):
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
