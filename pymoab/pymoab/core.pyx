"""Implements core functionality."""
from cython.operator cimport dereference as deref

cimport numpy as np
import numpy as np
import ctypes

from pymoab cimport moab
from .tag cimport Tag, TagArray
from .rng cimport Range
from .types import check_error, np_tag_type, validate_type, _convert_array, _eh_array
from . import types
from libcpp.vector cimport vector
from libc.stdlib cimport malloc

cdef void* null = NULL

cdef class Core(object):

    def __cinit__(self):
        """ Constructor """
        self.inst = new moab.Core()

    def __del__(self):
        """ Destructor """
        del self.inst

    def impl_version(self):
        """MOAB implementation number as a float."""
        return self.inst.impl_version()

    def load_file(self, str fname, file_set = None, exceptions = ()):
        """
        Load or import a file.

        Load a MOAB-native file or import data from some other supported file format.

        Example
        -------
        mb = pymoab.core.Core()
        mb.load_file("/location/of/my_file.h5m")

        OR 

        fs = mb.create_meshset()
        mb = pymoab.core.Core()
        mb.load_file("/location/of/my_file.h5m", fs)

        Parameters
        ----------
        fname : string
            The location of the file to read.
        file_set : EntityHandle 
            If not None, this argument must be a valid
            entity set handle. All entities read from the file will be added to
            this set. File metadata will be added to tags on the set.
        exceptions : tuple 
            A tuple containing any error types that should
            be ignored. See pymoab.types module for more info.

        Returns
        -------
          None.

        Raises
        ------
          MOAB ErrorCode 
              if a MOAB error occurs
          ValueError
              if a parameter is not of the correct type
        """
        cfname = fname.decode()
        cdef moab.EntityHandle fset
        cdef moab.EntityHandle* ptr
        if file_set != None:
            fset = file_set
            ptr = &fset
        else:
            ptr = NULL
        cdef const char * file_name = cfname
        cdef moab.ErrorCode err = self.inst.load_file(fname,ptr)
        check_error(err, exceptions)

    def write_file(self, str fname, exceptions = ()):
        """Write sthe MOAB data to a file."""
        cfname = fname.decode()
        cdef const char * file_name = cfname
        cdef moab.ErrorCode err = self.inst.write_file(fname)
        check_error(err, exceptions)

    def create_meshset(self, unsigned int options = 0x02, exceptions = ()):
        cdef moab.EntityHandle ms_handle = 0
        cdef moab.EntitySetProperty es_property = <moab.EntitySetProperty> options
        cdef moab.ErrorCode err = self.inst.create_meshset(es_property, ms_handle)
        check_error(err, exceptions)
        return ms_handle

    def add_entities(self, moab.EntityHandle ms_handle, entities, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        if isinstance(entities, Range):
           r = entities
           err = self.inst.add_entities(ms_handle, deref(r.inst))
        else:
           arr = _eh_array(entities)
           err = self.inst.add_entities(ms_handle, <unsigned long*> arr.data, len(entities))
        check_error(err, exceptions)

    def remove_entities(self, moab.EntityHandle ms_handle, entities, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        if isinstance(entities, Range):
            r = entities
            err = self.inst.remove_entities(ms_handle, deref(r.inst))
        else:
            arr = _eh_array(entities)
            err = self.inst.remove_entities(ms_handle, <unsigned long*> arr.data, len(entities))
        check_error(err, exceptions)

    def delete_entities(self, entities, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        if isinstance(entities, Range):
            r = entities
            err = self.inst.delete_entities(deref(r.inst))
        else:
            arr = _eh_array(entities)
            err = self.inst.delete_entities(<unsigned long*> arr.data, len(entities))
        check_error(err, exceptions)
        

    def create_vertices(self, coordinates, exceptions = ()):
        cdef Range rng = Range()
        cdef np.ndarray[np.float64_t, ndim=1] coords_arr = _convert_array(coordinates, [np.float64], np.float64)
        assert len(coordinates)%3 == 0, "Incorrect number of coordinates provided."
        cdef moab.ErrorCode err = self.inst.create_vertices(<double *> coords_arr.data,
                                                            len(coordinates)//3,
                                                            deref(rng.inst))
        check_error(err, exceptions)
        return rng

    def create_element(self, int t, connectivity, exceptions = ()):
        cdef moab.EntityType typ = <moab.EntityType> t
        cdef moab.EntityHandle handle = 0
        cdef np.ndarray[np.uint64_t, ndim=1] conn_arr = _eh_array(connectivity)
        cdef int nnodes = len(connectivity)
        cdef moab.ErrorCode err = self.inst.create_element(typ,
            <unsigned long*> conn_arr.data, nnodes, handle)
        check_error(err, exceptions)
        return handle

    def create_elements(self, int t, np.ndarray[np.uint64_t, ndim=2] connectivity, exceptions = ()):
        cdef int i
        cdef moab.ErrorCode err
        cdef moab.EntityType typ = <moab.EntityType> t
        #cdef moab.EntityHandle handle = 0
        cdef int nelems = connectivity.shape[0]
        cdef int nnodes = connectivity.shape[1]
        cdef np.ndarray[np.uint64_t, ndim=1] connectivity_i
        cdef np.ndarray[np.uint64_t, ndim=1] handles = np.empty(nelems, 'uint64')
        for i in range(nelems):
            connectivity_i = connectivity[i]
            err = self.inst.create_element(typ, <unsigned long*> connectivity_i.data,
                                           nnodes, deref((<unsigned long*> handles.data)+i))
            check_error(err, exceptions)
        return handles

    def tag_get_handle(self,
                       const char* name,
                       size = None,
                       tag_type = None,
                       create_if_missing = False,
                       storage_type = types.MB_TAG_DENSE,
                       exceptions = ()):
        cdef Tag tag = Tag()
        cdef moab.ErrorCode err
        cdef moab.DataType tt
        cdef int s
        err = self.inst.tag_get_handle(name, tag.inst)
        if err == types.MB_TAG_NOT_FOUND and create_if_missing:
            if tag_type is None or size is None:
                print "ERROR: Not enough information provided to create tag."
                raise ValueError
            else:
                tt = tag_type
                s = size
                err = self.inst.tag_get_handle(name, s, tt, tag.inst, storage_type|types.MB_TAG_CREAT)
        check_error(err, exceptions)
        return tag

    def tag_set_data(self, Tag tag, entity_handles, data, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef moab.DataType tag_type = moab.MB_MAX_DATA_TYPE
        err = self.inst.tag_get_data_type(tag.inst, tag_type);
        check_error(err, ())
        cdef int length = 0
        err = self.inst.tag_get_length(tag.inst,length);
        check_error(err,())
        cdef np.ndarray data_arr = np.asarray(data)
        #if the data array is not flat it must be dimension 2 and have
        #as many entries as entity handles provided
        if data_arr.ndim > 1:
            assert data_arr.ndim == 2
            assert data_arr.shape[0] == len(entity_handles)
            #each entry must be equal to the tag length as well
            for entry in data_arr:
                len(entry) == length
            #if all of this is true, then flatten the array and continue
            data_arr = data_arr.flatten()
        error_str = "Incorrect data length"
        if types.MB_TYPE_OPAQUE == tag_type:
            assert data_arr.size == len(entity_handles), error_str
        else:
            assert data_arr.size == len(entity_handles)*length, error_str
        data_arr = validate_type(tag_type,length,data_arr)
        if isinstance(entity_handles,Range):
            r = entity_handles
            err = self.inst.tag_set_data(tag.inst, deref(r.inst), <const void*> data_arr.data)
            check_error(err, exceptions)
        else:
            arr = _eh_array(entity_handles)
            err = self.inst.tag_set_data(tag.inst, <unsigned long*> arr.data, len(entity_handles), <const void*> data_arr.data)
            check_error(err, exceptions)

    def tag_get_data(self, Tag tag, entity_handles, flat = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef moab.DataType tag_type = moab.MB_MAX_DATA_TYPE
        err = self.inst.tag_get_data_type(tag.inst, tag_type);
        check_error(err,())
        cdef int length = 0
        err = self.inst.tag_get_length(tag.inst,length);
        check_error(err,())
        cdef np.ndarray data
        if tag_type is types.MB_TYPE_OPAQUE:
            data = np.empty((len(entity_handles),),dtype='S'+str(length))
        else:
            data = np.empty((length*len(entity_handles),),dtype=np.dtype(np_tag_type(tag_type)))
        if isinstance(entity_handles,Range):
            r = entity_handles
            err = self.inst.tag_get_data(tag.inst, deref(r.inst), <void*> data.data)
            check_error(err, exceptions)
        else:
            arr = _eh_array(entity_handles)
            err = self.inst.tag_get_data(tag.inst, <unsigned long*> arr.data, len(entity_handles), <void*> data.data)
            check_error(err,exceptions)
        if flat:
            return data
        else:
            entry_len = 1 if tag_type == types.MB_TYPE_OPAQUE else length
            return data.reshape((len(entity_handles),entry_len))

    def get_adjacencies(self, entity_handles, int to_dim, bint create_if_missing = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef Range adj = Range()
        if isinstance(entity_handles, Range):
            r = entity_handles
            err = self.inst.get_adjacencies(deref(r.inst), to_dim, create_if_missing, deref(adj.inst))
        else:
            arr = _eh_array(entity_handles)
            err = self.inst.get_adjacencies(<unsigned long*> arr.data, len(entity_handles), to_dim, create_if_missing, deref(adj.inst))
        check_error(err, exceptions)
        return adj

    def type_from_handle(self, entity_handle):
        cdef moab.EntityType t
        t = self.inst.type_from_handle(<unsigned long> entity_handle)
        return t

    def get_child_meshsets(self, meshset_handle, num_hops = 1, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r = Range()
        err = self.inst.get_child_meshsets(<unsigned long> meshset_handle, deref(r.inst), num_hops)
        check_error(err, exceptions)
        return r

    def get_parent_meshsets(self, meshset_handle, num_hops = 1, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r = Range()
        err = self.inst.get_parent_meshsets(<unsigned long> meshset_handle, deref(r.inst), 0)
        check_error(err, exceptions)
        return r

    def add_parent_meshset(self, child_meshset, parent_meshset, exceptions = ()):
        cdef moab.ErrorCode err
        err = self.inst.add_parent_meshset(<unsigned long> child_meshset, <unsigned long> parent_meshset)
        check_error(err, exceptions)
    
    def add_child_meshset(self, parent_meshset, child_meshset, exceptions = ()):
        cdef moab.ErrorCode err
        err = self.inst.add_child_meshset(<unsigned long> parent_meshset, <unsigned long> child_meshset)
        check_error(err, exceptions)

    def add_parent_child(self, parent_meshset, child_meshset, exceptions = ()):
        cdef moab.ErrorCode err
        err = self.inst.add_parent_child(<unsigned long> parent_meshset, <unsigned long> child_meshset)
        check_error(err, exceptions)

    def get_root_set(self):
        return <unsigned long> 0

    def get_coords(self, entities, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef np.ndarray coords
        if isinstance(entities, Range):
            r = entities
            coords = np.empty((3*r.size(),),dtype='float64')
            err = self.inst.get_coords(deref(r.inst), <double*> coords.data)
        else:
            arr = _eh_array(entities)
            coords = np.empty((3*len(arr),),dtype='float64')
            err = self.inst.get_coords(<unsigned long*> arr.data, len(entities), <double*> coords.data)
        check_error(err, exceptions)
        return coords

    def get_entities_by_type(self, meshset, t, bint recur = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef vector[moab.EntityHandle] entities
        cdef moab.EntityType typ = t
        err = self.inst.get_entities_by_type(<unsigned long> meshset,
                                             typ,
                                             entities,
                                             recur)
        check_error(err, exceptions)
        return entities

    def get_entities_by_type_and_tag(self,
                                     meshset,
                                     t,
                                     tags,
                                     values,
                                     int condition = types.INTERSECT,
                                     bint recur = False,
                                     exceptions = ()):
        cdef np.ndarray vals = np.asarray(values, dtype='O')
        # overall dimension of the numpy array should be 2
        # one array for each tag passed to the function
        assert vals.ndim == 2
        #ensure that the correct number of arrays exist
        assert len(tags) == vals.shape[0]
        #setup some initial variables to pass to MOAB function
        cdef int num_tags = len(tags)
        cdef moab.EntityType typ = t
        #create tag array to pass to function
        cdef TagArray ta = TagArray(tags)
        #allocate memory for an appropriately sized void** array
        cdef void** arr = <void**> malloc(num_tags*sizeof(void*))
        
        #get the tag type
        cdef moab.ErrorCode err
        cdef moab.DataType this_tag_type = moab.MB_MAX_DATA_TYPE
        cdef int this_tag_length = 0
        cdef Tag this_tag
        cdef bytes val_str
        cdef np.ndarray this_data
        
        # assign values to void** array
        for i in range(num_tags):
            # make sure we're dealing with a 1-D array at this point
            assert vals[i].ndim == 1
            # if None is passed, set pointer to NULL and continue
            if vals[i][0] == None:
                arr[i] = NULL
            # otherwise get the tag type
            else:
                this_data = vals[i]
                this_tag = tags[i]
                err = self.inst.tag_get_data_type(this_tag.inst, this_tag_type)
                check_error(err)
                #and length
                err = self.inst.tag_get_length(this_tag.inst,this_tag_length);
                check_error(err)
                #check that the data for this tag is the correct length
                if this_tag_type == moab.MB_TYPE_OPAQUE:
                    #if this is an opaque tag, there should only be one
                    #string entry in the array
                    assert this_data.size == 1
                else:
                    assert this_data.size == this_tag_length
                #validate the array type and convert the dtype if necessary
                this_data = validate_type(this_tag_type, this_tag_length, this_data)
                #set the array value
                arr[i] = <void*> this_data.data
            
        #a range to hold returned entities
        cdef Range ents = Range()
        #here goes nothing
        err = self.inst.get_entities_by_type_and_tag(<unsigned long> meshset,
                                                     typ,
                                                     ta.ptr,
                                                     <const void**> arr,
                                                     len(tags),
                                                     deref(ents.inst),
                                                     condition,
                                                     recur)
        check_error(err, exceptions)
        # return entities found in the MOAB function call
        return ents


    def get_entities_by_handle(self, meshset, bint recur = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        err = self.inst.get_entities_by_handle(<unsigned long> meshset, deref(ents.inst), recur)
        check_error(err, exceptions)
        return ents

    def get_entities_by_dimension(self, meshset, int dimension, bint recur = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        err = self.inst.get_entities_by_dimension(<unsigned long> meshset, dimension, deref(ents.inst), recur)
        check_error(err, exceptions)
        return ents            
