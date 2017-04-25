"""Implements range functionality."""
from cython.operator cimport dereference as deref

from pymoab cimport moab
from .types import _eh_array

def intersect(Range r1, Range r2):
    r = Range()
    cdef moab.Range i = moab.intersect(deref(r1.inst),deref(r2.inst))
    r.inst.merge(i)
    return r

def subtract(Range r1, Range r2):
    r = Range()
    cdef moab.Range i = moab.subtract(deref(r1.inst),deref(r2.inst))
    r.inst.merge(i)
    return r

def unite(Range r1, Range r2):
    r = Range()
    cdef moab.Range i = moab.unite(deref(r1.inst),deref(r2.inst))
    r.inst.merge(i)
    return r

cdef class Range(object):

    #def __cinit__(self, moab.EntityHandle val1=None, moab.EntityHandle val2=None):
    #    if val1 is None or val2 is None:
    #        self.inst = new moab.Range()
    #    else:
    #        self.inst = new moab.Range(val1, val2)

    def __cinit__(self, arg = None):
        self.inst = new moab.Range()
        if arg is None:
            return
        if type(arg) is long:
            self.inst.insert(arg)
        #hack to copy
        elif isinstance(arg, Range):
            for eh in arg:
                self.inst.insert(eh)
        #create from iterable
        elif arg is not None:
            entity_array = _eh_array(arg)
            for eh in entity_array:
                self.inst.insert(eh)
        else:
            raise ValueError, "Not a valid argument to Range constructor."
                

    def __del__(self):
        del self.inst

    def size(self):
        """The number of values this Ranges represents."""
        return len(self)

    def __len__(self):
        return self.inst.size()

    def psize(self):
        """The number of range pairs in the list."""
        return self.inst.psize()

    def empty(self):
        """Is the range empty?"""
        return self.inst.empty()

    def clear(self):
        """clears the contents of the list."""
        self.inst.clear()

    def erase(self, moab.EntityHandle val):
         self.inst.erase(val)

    def pop_front(self):
        return self.inst.pop_front()

    def pop_back(self):
        return self.inst.pop_back()
    
    def all_of_type(self, moab.EntityType t):
        return self.inst.all_of_type(t)

    def all_of_dimension(self, int dim):
        return self.inst.all_of_dimension(dim)

    def num_of_dimension(self, int dim):
        return self.inst.num_of_dimension(dim)

    def num_of_type(self, moab.EntityType t):
        return self.inst.num_of_type(t)

    def insert(self, moab.EntityHandle ent_handle):
        self.inst.insert(ent_handle)

    def merge(self, other):
        cdef Range r
        if isinstance(other, Range):
            r = other
            self.inst.merge(deref(r.inst))
        else:
            raise ValueError("Operation not valie for non-Range")
        
    def unite(self, other):
        cdef Range r
        if isinstance(other, Range):
            r = other

        
    def __iter__(self):
        cdef int i = 0
        for i in range(0, self.inst.size()):
            yield self[i]

    def __getitem__(self, key):
        cdef moab.EntityHandle rtn
        if isinstance(key, int):
            i = key if key >= 0 else len(self)+key
            rtn = deref(self.inst)[i]
            if i < self.size():
                return rtn
            else:
                raise StopIteration
        elif isinstance(key, slice):
            step = key.step if key.step is not None else 1
            start = key.start if key.start is not None else 0
            stop = key.stop if key.stop is not None else len(self)
            ents = list(self)[start:stop:step]
            return Range(ents)
        else:
            raise ValueError
        
    def __str__(self):
        prefix= "["
        delim = ", "
        suffix= "]"
        out_str=prefix
        for entity in self:
            out_str+=str(entity)
            out_str+=delim
        #replace final delimiter
        out_str=out_str[:-2]
        out_str+=suffix
        return out_str
            
