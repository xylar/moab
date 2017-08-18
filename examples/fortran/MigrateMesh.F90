
program MigrateMesh
  implicit none

#include "moab/MOABConfig.h"
#ifdef MOAB_HAVE_MPI
#  include "mpif.h"
#else
#  error "enable parallel build"
#endif


    ! init the parallel partition
    integer ierr, sz, rank, color, new_rank
    integer  newComm
    integer pid !  this is for physics id
    integer nghlay ! number of ghost layers for loading
    character*10 appname
    character*132 readopts
    character*132 filename
    integer iMOAB_InitializeFortran, iMOAB_RegisterFortranApplication, iMOAB_LoadMesh


    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, sz, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    if (rank .eq. 0) print *, "size:", sz
    call errorout(ierr, 'cannot get rank' )
    if (sz <= 3) then
      print *, " run on at least 4 processes "
      call exit(1)
    endif
    ! split the communicator in 2
    ! the mesh will be read on first n-2 tasks, and migrated on the last 2 tasks

    if (rank < sz-2) then
      color = 0
    else
      color = 1
    endif

    call MPI_Comm_split ( MPI_COMM_WORLD, color, rank, newComm, ierr )
    !    print *, "rank, color :", rank, color, " ierr:", ierr
    call errorout(ierr, 'did not split communicators' )

    call MPI_COMM_RANK(newComm, new_rank, ierr)
    call errorout(ierr, 'did not get local rank' )

    if (new_rank == 0) print *, "global rank ", rank, " of ", sz, " color ", color, " new rank :" , new_rank

    ierr = iMOAB_InitializeFortran()
    call errorout(ierr, 'did not initialize fortran' )
    if (rank == 0) print *, "initialize fortran"

    if (color .eq. 1) then
       appname='phis1'
       ierr = iMOAB_RegisterFortranApplication(appname, newComm, pid)
       print *, ' register ', appname, " on rank ", rank, " pid ", pid
    else
       appname = 'phis0'
       ierr = iMOAB_RegisterFortranApplication(appname, newComm, pid)
       print *, ' register ', appname, " on rank ", rank, " pid ", pid
    endif


    if (color .eq. 1 ) then
       filename = 'spherecube.h5m'//CHAR(0)
       readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS'//CHAR(0)
        if (new_rank == 0) print *, "loading " , trim(filename) , " with options " , trim(readopts)
       nghlay = 0

       ierr = iMOAB_LoadMesh(pid, trim(filename), trim(readopts), nghlay)
       if (new_rank.eq.0) print *, "loaded in parallel ", filename, " error: ", ierr
    endif

    call MPI_Comm_free(newComm, ierr)
    call errorout(ierr, 'did not free communicator' )

    call MPI_Finalize(ierr)
    call errorout(ierr, 'did not finalize MPI' )
contains
  SUBROUTINE errorout(ierr, message)
  integer ierr
  character*(*) message
  if (ierr.ne.0) then
    print *, message
    call exit (1)
  end if
  end

end program MigrateMesh
