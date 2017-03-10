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
        file_set : EntityHandle (default None)
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
        """
        Write or export a file.
        
        Write a MOAB-native file or export data to some other supported file format.
        
        Example
        -------
        mb.write_file("new_file.h5m")
        
        Parameters
        ----------
        fname : string
            thhe location of the file to write
        exceptions : tuple (default is empty tuple)
            tuple containing any error types that should
            be ignored (see pymoab.types module for more info)
        
        Returns
        -------
        None
        
        Raises
        ------
        MOAB ErrorCode 
            if a MOAB error occurs
        ValueError
            if a parameter is not of the correct type
        """
        cfname = fname.decode()
        cdef const char * file_name = cfname
        cdef moab.ErrorCode err = self.inst.write_file(fname)
        check_error(err, exceptions)

    def create_meshset(self, unsigned int options = 0x02, exceptions = ()):
        """
        Create a new mesh set.
        
        Create a new mesh set. Meshsets can store entities ordered or
        unordered. A set can include entities at most once (MESHSET_SET) or more
        than once. Meshsets can optionally track its members using adjacencies
        (MESHSET_TRACK_OWNER); if set, entities are deleted from tracking
        meshsets before being deleted. This adds data to mesh entities, which
        can be more memory intensive.
        
        Example
        -------
        # Create a standard meshset
        mb = pymoab.core.Core()
        new_meshset = mb.create_meshset()
        
        # Create a meshset which tracks entities
        mb = pymoab.core.Core()
        new_meshset = mb.create_meshset(pymoab.types.MESHSET_TRACK_OWNER)

        Parameters
        ----------
        options : MOAB EntitySet Property (default is MESHSET_SET)
            settings for the Meshset being created
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)
        
        Returns
        -------
        MOAB Meshset (EntityHandle)
        
        Raises
        ------
        MOAB ErrorCode 
            if a MOAB error occurs
        """        
        cdef moab.EntityHandle ms_handle = 0
        cdef moab.EntitySetProperty es_property = <moab.EntitySetProperty> options
        cdef moab.ErrorCode err = self.inst.create_meshset(es_property, ms_handle)
        check_error(err, exceptions)
        return ms_handle

    def add_entities(self, moab.EntityHandle ms_handle, entities, exceptions = ()):
        """
        Add entities to the specified meshset. Entities can be provided either
        as a pymoab.rng.Range o bject or as an iterable of EntityHandles.

        If meshset has MESHSET_TRACK_OWNER option set, adjacencies are also
        added to entities in entities.

        Example
        -------
        vertices # a iterable of MOAB vertex handles
        new_meshset = mb.create_meshset()
        mb.add_entities(new_meshset, entities)
        
        Parameters
        ----------
        ms_handle : EntityHandle
            EntityHandle of the Meshset the entities will be added to
        entities : Range or iterable of EntityHandles
            Entities that to add to the Meshset            
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)
        
        Returns
        -------
        None
        
        Raises
        ------
        MOAB ErrorCode 
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """        
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
        """
        Remove entities from the specified meshset. Entities can be provided either
        as a pymoab.rng.Range o bject or as an iterable of EntityHandles.

        If meshset has MESHSET_TRACK_OWNER option set, adjacencies in entities
        in entities are updated.

        Example
        -------
        vertices # iterable of MOAB vertex handles
        new_meshset = mb.create_meshset()
        mb.add_entities(new_meshset, entities)
        # keep only the first entity in this meshset
        mb.remove_entities(new_meshset, entities[1:])

        Parameters
        ----------
        ms_handle : EntityHandle
            EntityHandle of the Meshset the entities will be removed from
        entities : Range or iterable of EntityHandles
            Entities that to remove from the Meshset
        exceptions : tuple (default is empty tuple)
        A tuple containing any error types that should
        be ignored. (see pymoab.types module for more info)

        Returns
        -------
        None

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """
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
        """
        Delete entities from the database.

        If any of the entities are contained in any meshsets, it is removed from
        those meshsets which were created with MESHSET_TRACK_OWNER option

        Example
        -------
        mb = pymoab.core.Core()
        # generate mesh using other Core functions #
        mb.delete_entities(entities)

        Parameters
        ----------
        entities : Range or iterable of EntityHandles
            Entities that to delete from the database
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        None

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """
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
        """
        Create vertices using the specified x,y,z coordinates.

        Example
        -------
        mb = pymoab.core.Core()
        #create two vertices
        coords = np.array((0, 0, 0, 1, 1, 1),dtype='float64')
        verts = mb.create_vertices(coords)
        
        OR

        mb = pymoab.core.Core()
        #create two vertices
        coords = [[0.0, 0.0, 0.0],[1.0, 1.0, 1.0]]
        verts = mb.create_vertices(coords)       

        Parameters
        ----------
        coordinates : 1-D or 2-D iterable
            Coordinates can be provided as either a 1-D iterable of 3*n floats
            or as a 2-D array with dimensions (n,3) where n is the number of
            vertices being created.
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        None

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """
        cdef Range rng = Range()
        cdef np.ndarray coords_as_arr = np.asarray(coordinates)
        if coords_as_arr.ndim > 1:
            assert coords_as_arr.ndim == 2
            coords_as_arr = coords_as_arr.flatten()
            print len(coords_as_arr)
        cdef np.ndarray[np.float64_t, ndim=1] coords_arr = _convert_array(coords_as_arr, [float, np.float64], np.float64)
        assert len(coordinates)%3 == 0, "Incorrect number of coordinates provided."
        cdef moab.ErrorCode err = self.inst.create_vertices(<double *> coords_arr.data,
                                                            len(coords_arr)//3,
                                                            deref(rng.inst))
        check_error(err, exceptions)
        return rng

    def create_element(self, int entity_type, connectivity, exceptions = ()):
        """
        Create an elment of type, entity_type, using vertex EntityHandles in connectivity.

        Example
        -------
        mb = core.Core()
        # create some vertices
        coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
        verts = mb.create_vertices(coords)

        # create the triangle in the database
        tri_verts = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
        tris = mb.create_element(types.MBTRI,tri_verts)
        
        OR

        tri_verts = [verts[0],verts[1],verts[2]]
        tris = mb.create_element(types.MBTRI,tri_verts)


        Parameters
        ----------
        entity_type : MOAB EntityType (see pymoab.types module)
            type of entity to create (MBTRI, MBQUAD, etc.)
        coordinates : iterable of EntityHandles
            1-D array-like iterable of vertex EntityHandles the element is to be
            created from
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        EntityHandle of element created

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """
        cdef moab.EntityType typ = <moab.EntityType> entity_type
        cdef moab.EntityHandle handle = 0
        if isinstance(connectivity,Range):
            connectivity = list(connectivity)
        cdef np.ndarray[np.uint64_t, ndim=1] conn_arr = _eh_array(connectivity)
        cdef int nnodes = len(connectivity)
        cdef moab.ErrorCode err = self.inst.create_element(typ,
            <unsigned long*> conn_arr.data, nnodes, handle)
        check_error(err, exceptions)
        return handle

    def create_elements(self, int entity_type, connectivity, exceptions = ()):
        """
        Create an elments of type, entity_type, using vertex EntityHandles in connectivity.

        Example
        -------
        mb = core.Core()
        # create some vertices for triangle 1
        coords1 = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
        verts1 = mb.create_vertices(coords)
        # create some more vertices for triangle 2
        coords2 = np.array((1,2,3,4,5,6,1,1,1),dtype='float64')
        verts2 = mb.create_vertices(coords)

        # create the triangles in the database
        tri_verts = [[verts1[0],verts1[1],verts1[2]],[verts2[0],verts2[1],verts2[2]]
        tris = mb.create_elements(types.MBTRI,tri_verts)

        OR

        tri_verts = [verts1,verts2]
        tris = mb.create_elements(types.MBTRI,tri_verts)

        Parameters
        ----------
        entity_type : MOAB EntityType (see pymoab.types module)
            type of entity to create (MBTRI, MBQUAD, etc.)
        coordinates : iterable of EntityHandles
            2-D array-like iterable of vertex EntityHandles the elements are to
            be created from. Each entry in first dimension of the array should
            have an appropriate length for the element type being created (MBTRI
            means each entry has length 3).
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        EntityHandles of element created as a Range

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """
        cdef int i
        cdef moab.ErrorCode err
        cdef moab.EntityType typ = <moab.EntityType> entity_type
        cdef np.ndarray connectivity_as_arr = np.asarray(connectivity)
        assert connectivity_as_arr.ndim == 2 #required for now
        cdef int nelems = connectivity_as_arr.shape[0]
        cdef int nnodes = connectivity_as_arr.shape[1]
        cdef np.ndarray[np.uint64_t, ndim=1] connectivity_i
        cdef np.ndarray[np.uint64_t, ndim=1] handles = np.empty(nelems, 'uint64')
        for i in range(nelems):
            connectivity_i = _eh_array(connectivity[i])
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
        """
        Retrieve or create tag handles for storing data on mesh elements, vertices,
        meshsets, etc.

        Example - retrieving an existing tag handle
        -------
        # new MOAB core instance
        mb = core.Core()
        # 
        tag_type = pymoab.types.MB_TYPE_INTEGER # define the tag's data type
        tag_size = 1 # the tag size (1 integer value)
        tag_handle = mb.tag_get_handle("NewDataTag", 
                                       tag_size, 
                                       tag_type)

        Example - creating a new tag handle
        -------
        # new MOAB core instance
        mb = core.Core()
        # 
        tag_type = pymoab.types.MB_TYPE_INTEGER # define the tag's data type
        tag_size = 1 # the tag size (1 integer value)
        tag_handle = mb.tag_get_handle("NewDataTag", 
                                       tag_size, 
                                       tag_type, 
                                       create_if_missing = True)

        Parameters
        ----------
        name : string
            name of the tag
        size : int
            number of values the tag contains (or the number of characters in
            the case of an opaque tag)
        tag_type : MOAB DataType
            indicates the data type of the tag (MB_TYPE_DOUBLE, MB_TYPE_INTEGER,
            MB_TYPE_OPAQUE, etc.)
        create_if_missing : bool (default False)
            indicates to the database instance that the tag should be created
            if it cannot be found
        storage_type : MOAB tag storage type (default MB_TYPE_DENSE)
            in advanced use of the database, this flag controls how this tag's 
            data is stored in memory. The two most common storage types are

            MB_TYPE_SPARSE -  sparse tags are stored as a list of (entity
                handle, tag value) tuples, one list per sparse tag, sorted by
                entity handle

            MB_TYPE_DENSE - Dense tag values are stored in arrays which match
                arrays of contiguous entity handles. Dense tags are more
                efficient in both storage and memory if large numbers of
                entities are assigned the same tag. Storage for a given dense
                tag is not allocated until a tag value is set on an entity;
                memory for a given dense tag is allocated for all entities in a
                given sequence at the same time.  Sparse: Sparse tags are stored
                as a list of (entity handle, tag value) tuples, one list per
                sparse tag, sorted by entity handle.

            MB_TYPE_BIT - Bit tags are stored similarly to dense tags, but with
                special handling to allow allocation in bit-size amounts per
                entity.
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        MOAB TagHandle

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        """        
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
        """
        Create an elments of type, entity_type, using vertex EntityHandles in connectivity.

        This function will ensure that the data is of the correct type and
        length based on the tag it is to be applied to. Data can be passed as a
        1-D iterable of appropriate length or as 2-D iterable with entries for
        each entity handle equal to the length of the tag.

        Example
        -------
        # create moab core instance
        mb = core.Core()
        # create some vertices
        coords = np.array((0,0,0,1,0,0,1,1,1),dtype='float64')
        verts = mb.create_vertices(coords)
        # create a new tag for data
        data_tag = mb.tag_get_handle("Data",
                                     2,
                                     pymoab.types.MB_TYPE_DOUBLE,
                                     create_if_missing = True)
        # some sample data
        data = np.array([1.0,2.0,3.0,4.0,5.0,6.0], dtype = 'float')
        #use this function to tag vertices with the data
        mb.tag_set_data(data_tag, verts, data)

        OR

        # some sample data
        data = np.array([[1.0,2.0],[3.0,4.0],[5.0,6.0]], dtype = 'float')
        #use this function to tag vertices with the data
        mb.tag_set_data(data_tag, verts, data)

        OR

        ### pass data as a an iterable object ###
        # some sample data
        data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        #use this function to tag vertices with the data
        mb.tag_set_data(data_tag, verts, data)

        OR

        ### pass data as a nested iterable object ###
        ### (can be convenient for vector tags) ###
        # some sample data
        data = [[1.0, 2.0],[3.0, 4.0],[4.0, 5.0]]
        #use this function to tag vertices with the data
        mb.tag_set_data(data_tag, verts, data)

        OR

        ### tag a single vertex ###
        # get a single vertex
        vertex_to_tag = verts[0]
        data = [1.0, 2.0]
        #use this function to tag vertices with the data
        mb.tag_set_data(data_tag, vertex, data)

        Parameters
        ----------
        tag : MOAB TagHandle
            tag to which the data is applied
        entity_handles : iterable of MOAB EntityHandle's or single EntityHandle
            the EntityHandle(s) to tag the data on. This can be any iterable of
            EntityHandles or a single EntityHandle.
        entity_type : MOAB EntityType (see pymoab.types module)
            type of entity to create (MBTRI, MBQUAD, etc.)
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        None

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if a data entry is not of the correct type or cannot be converted to
            the correct type
        """
        cdef moab.ErrorCode err
        cdef Range r
        cdef moab.DataType tag_type = moab.MB_MAX_DATA_TYPE
        if isinstance(entity_handles,Range):
            r = entity_handles
        else:
            r = Range(entity_handles)        
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
            assert data_arr.shape[0] == len(r)
            #each entry must be equal to the tag length as well
            for entry in data_arr:
                len(entry) == length
            #if all of this is true, then flatten the array and continue
            data_arr = data_arr.flatten()
        error_str = "Incorrect data length"
        if types.MB_TYPE_OPAQUE == tag_type:
            assert data_arr.size == len(r), error_str
        else:
            assert data_arr.size == len(r)*length, error_str
        data_arr = validate_type(tag_type,length,data_arr)
        err = self.inst.tag_set_data(tag.inst, deref(r.inst), <const void*> data_arr.data)
        check_error(err, exceptions)

    def tag_get_data(self, Tag tag, entity_handles, flat = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef moab.DataType tag_type = moab.MB_MAX_DATA_TYPE
        #create a range
        if isinstance(entity_handles,Range):
            r = entity_handles
        else:
            r = Range(entity_handles)
        # get the tag type and length for validation
        err = self.inst.tag_get_data_type(tag.inst, tag_type);
        check_error(err,())
        cdef int length = 0
        err = self.inst.tag_get_length(tag.inst,length);
        check_error(err,())
        #create array to hold data
        cdef np.ndarray data
        if tag_type is types.MB_TYPE_OPAQUE:
            data = np.empty((len(r),),dtype='S'+str(length))
        else:
            data = np.empty((length*len(r),),dtype=np.dtype(np_tag_type(tag_type)))
        err = self.inst.tag_get_data(tag.inst, deref(r.inst), <void*> data.data)
        check_error(err, exceptions)
        # return data as user specifies
        if flat:
            return data
        else:
            entry_len = 1 if tag_type == types.MB_TYPE_OPAQUE else length
            return data.reshape((len(r),entry_len))

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
        assert isinstance(tags, Tag) or len(tags) == 1 , "Only single-tag queries are currently supported."
        # overall dimension of the numpy array should be 2
        # one array for each tag passed to the function
        # assert vals.ndim == 2 (was for support of multiple tags)
        #ensure that the correct number of arrays exist        
        assert len(tags) == vals.shape[0]
        cdef int num_tags = len(tags)        
        #allocate memory for an appropriately sized void** array
        cdef void** arr = <void**> malloc(num_tags*sizeof(void*))
        #some variables to help in setting up the void array
        cdef moab.ErrorCode err
        cdef moab.DataType this_tag_type = moab.MB_MAX_DATA_TYPE
        cdef int this_tag_length = 0
        cdef Tag this_tag
        cdef bytes val_str
        cdef np.ndarray this_data

        # if this is a single tag and the values array is 1-D
        # then validate and call function
        if (num_tags == 1) and (vals.ndim == 1):
            this_data = self._validate_data_for_tag(tags[0],vals)
            arr[0] = <void*> this_data.data if this_data is not None else NULL
        elif vals.ndim == 2:
            # assign values to void** array
            for i in range(num_tags):
                this_data = self._validate_data_for_tag(tags[i],vals[i])
                #set the array value
                arr[i] = <void*> this_data.data if this_data is not None else NULL
            
        #create tag array to pass to function
        cdef TagArray ta = TagArray(tags)
        #convert type to tag type
        cdef moab.EntityType typ = t        
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

    def _validate_data_for_tag(self, Tag tag, data_arr):
        cdef moab.ErrorCode err
        cdef moab.DataType tag_type = moab.MB_MAX_DATA_TYPE
        cdef int tag_length = 0
        # make sure we're dealing with a 1-D array at this point    
        assert data_arr.ndim == 1
        # if None is passed, set pointer to NULL and continue
        if data_arr[0] == None:
            return None
        else:
            # otherwise get the tag type            
            err = self.inst.tag_get_data_type(tag.inst, tag_type)
            check_error(err)
            #and length
            err = self.inst.tag_get_length(tag.inst,tag_length);
            check_error(err)
            #check that the data_arr for this tag is the correct length
            if tag_type == moab.MB_TYPE_OPAQUE:
                #if this is an opaque tag, there should only be one
                #string entry in the array
                assert data_arr.size == 1
            else:
                assert data_arr.size == tag_length
        #validate the array type and convert the dtype if necessary
        data_arr = validate_type(tag_type, tag_length, data_arr)
        return data_arr
