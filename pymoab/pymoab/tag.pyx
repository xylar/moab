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

    def get_data_type(self): 
        return self.inst.get_data_type()

cdef class TagArray(object): 
     def __cinit__(self, tags):
         cdef int num_tags = len(tags)
         self.inst = <moab.TagInfo**> malloc(num_tags*sizeof(moab.TagInfo*))
         cdef int i = 0
         cdef Tag t
         for i in range(0,num_tags):
             t = <Tag> tags[i]
             self.inst[i] = t.inst            
         self.ptr = self.inst
         
     def __del__(self):
         free(self.inst)
