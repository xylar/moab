"""Implements range functionality."""
from cython.operator cimport dereference as deref
from pymoab cimport moab
from .types import _eh_array, _eh_py_type

cdef void *null = NULL

def intersect(Range r1, Range r2):
    """
    Returns a range that is the intersection of r1 and r2.
    """
    r = Range()
    cdef moab.Range i = moab.intersect(deref(r1.inst),deref(r2.inst))
    r.inst.merge(i)
    return r

def subtract(Range r1, Range r2):
    """
    Returns a range that is the subtraction of r2 from r1.
    """
    r = Range()
    cdef moab.Range i = moab.subtract(deref(r1.inst),deref(r2.inst))
    r.inst.merge(i)
    return r

def unite(Range r1, Range r2):
    """
    Returns a range that is the union of r1 and r2.
    """
    r = Range()
    cdef moab.Range i = moab.unite(deref(r1.inst),deref(r2.inst))
    r.inst.merge(i)
    return r

cdef class Range(object):

    def __cinit__(self, arg = None):
        """
        Constructor.

        Accepts either a range or an iterable of EntityHandles.

        If no argument is provided, an empty Range will be created and returned.
        """
        self.inst = new moab.Range()
        if arg is None:
            return
        if isinstance(arg, _eh_py_type):
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
            raise ValueError("Not a valid argument to Range constructor.")


    def __del__(self):
        """
        Destructor.
        """
        del self.inst

    def size(self):
        """The number of values this Range represents."""
        return len(self)

    def __len__(self):
        """The number of values this Range represents."""
        return self.inst.size()

    def psize(self):
        """The number of range pairs in the list."""
        return self.inst.psize()

    def empty(self):
        """Is the range empty?"""
        return self.inst.empty()

    def clear(self):
        """Clears the contents of the Range."""
        self.inst.clear()

    def erase(self, moab.EntityHandle val):
        """Removes the EntityHandle, val, from the Range if present."""
        self.inst.erase(val)

    def pop_front(self):
        """Removes the front-most EntityHandle in the Range and returns the EntityHandle."""
        return _eh_py_type(self.inst.pop_front())

    def pop_back(self):
        """Removes the back-most EntityHandle in the Range and returns the EntityHandle."""
        return _eh_py_type(self.inst.pop_back())

    def all_of_type(self, moab.EntityType t):
        """Returns True if all EntityHandles in the Range represent mesh entities of
        EntityType, t, and False otherwise."""
        return self.inst.all_of_type(t)

    def all_of_dimension(self, int dim):
        """Returns True if all EntityHandles in the Range represent mesh entities of
        of dimension, dim, and False otherwise."""
        return self.inst.all_of_dimension(dim)

    def num_of_dimension(self, int dim):
        """Returns the number of EntityHandles with dimension, dim, in the Range."""
        return self.inst.num_of_dimension(dim)

    def num_of_type(self, moab.EntityType t):
        """Returns the number of EntityHandles with EntityType, t, in the Range."""
        return self.inst.num_of_type(t)

    def insert(self, moab.EntityHandle eh):
        """Inserts the EntityHandle, eh, into the Range."""
        self.inst.insert(eh)

    def merge(self, other):
        """Merges this Range with another Range, other."""
        cdef Range r
        if isinstance(other, Range):
            r = other
            self.inst.merge(deref(r.inst))
        else:
            raise ValueError("Operation not valid for non-Range object")

    def contains(self, other):
        """Checks if this range contains all contents of another Range, other"""
        cdef Range r
        if isinstance(other, Range):
            r = other
            return self.inst.contains(deref(r.inst))
        else:
            raise ValueError("Operation not valid for non-Range object")

    def subset_by_type(self, t):
        """Returns the subset range containing only entities with the specified type"""
        cdef moab.Range mbr
        cdef moab.EntityType typ = t
        mbr = self.inst.subset_by_type(typ)
        cdef Range r = Range()
        # kind of a hack? allows one to return
        # a new range though
        r.inst.merge(mbr)
        return r

    def subset_by_dimension(self, int dim):
        """Returns the subset range containing only entities with the specified dimension"""
        cdef moab.Range mbr
        mbr = self.inst.subset_by_dimension(dim)
        cdef Range r = Range()
        # kind of a hack? allows one to return
        # a new range though
        r.inst.merge(mbr)
        return r

    def __iter__(self):
        """
        Iterator
        """
        cdef int i = 0
        for i in range(0, self.inst.size()):
            yield _eh_py_type(self[i])

    def __getitem__(self, key):
        """
        Index operator.
        """
        cdef moab.EntityHandle rtn
        if isinstance(key, int):
            i = key if key >= 0 else len(self)+key
            rtn = deref(self.inst)[i]
            if i < self.size():
                return _eh_py_type(rtn)
            else:
                raise StopIteration
        elif isinstance(key, slice):
            step = key.step if key.step is not None else 1
            start = key.start if key.start is not None else 0
            stop = key.stop if key.stop is not None else len(self)
            ents = list(self)[start:stop:step]
            return Range(ents)
        else:
            raise ValueError("Invalid key (type: {}) provided.".format(type(key)))

    def __richcmp__(self, other, op):
        cdef Range r
        if isinstance(other, Range):
            r = other
            result1 = subtract(self, r)
            result2 = subtract(r, self)
            result1.merge(result2)
            if op == 2: # ==
                return result1.empty()
            if op == 3: # !=
                return not result1.empty()
            else:
                NotImplementedError("This comparator isn't supported for Ranges at this time.")
        else:
            ValueError("Other is not a Range object. Cannot compare.")

    def __str__(self):
        """
        Range as a string
        """
        sout = self.inst.str_rep()
        return sout.decode()

    def __repr__(self):
        """
        Representation of class as a string
        """
        return self.__str__()
