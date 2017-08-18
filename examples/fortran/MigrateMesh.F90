

#define ERROR(rval) if (0 .ne. rval) call exit(1)

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
    character*10 appname
    integer iMOAB_InitializeFortran, iMOAB_RegisterFortranApplication


    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, sz, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    if (rank .eq. 0) print *, "size:", sz
    ERROR(ierr)
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
    ERROR(ierr)

    call MPI_COMM_RANK(newComm, new_rank, ierr)
    ERROR(ierr)

    if (new_rank == 0) print *, "global rank ", rank, " of ", sz, " color ", color, " new rank :" , new_rank

    ierr = iMOAB_InitializeFortran()
    ERROR(ierr)
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


end program MigrateMesh
