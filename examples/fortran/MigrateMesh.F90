
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
    integer comm1, comm2
    integer pid !  this is for physics id
    integer nghlay ! number of ghost layers for loading
    integer gr1(2)   !  2 processes in the first group
    character*10 appname
    character*132 readopts
    character*132 filename
    character*132 outfile
    character*132 wopts
    integer allgroup, group1, group2 ! Corresponding to MPI_Group in C
    integer tagcomm1, tagcomm2
    integer iMOAB_InitializeFortran, iMOAB_RegisterFortranApplication
    integer iMOAB_LoadMesh, iMOAB_SendElements, iMOAB_ReceiveElements, iMOAB_WriteMesh


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

    ! create new MPI groups for processors sz-2, sz-1 (group 1) and 0, 1, .., sz-3 (group 2)

    call MPI_COMM_GROUP (MPI_COMM_WORLD, allgroup, ierr)
    call errorout(ierr, 'cannot get group' )
    gr1(1) = sz-2
    gr1(2) = sz-1
    call MPI_Group_incl(allgroup, 2, gr1, group1, ierr)
    call errorout(ierr, 'cannot create group 1' )

    call MPI_Group_excl(allgroup, 2, gr1, group2, ierr)
    call errorout(ierr, 'cannot create group 2' )

    tagcomm1 = 1
    call MPI_Comm_create_group(MPI_COMM_WORLD, group1, tagcomm1, comm1, ierr)
    call errorout(ierr, 'cannot create communicator 1' )

    tagcomm2 = 2
    call MPI_Comm_create_group(MPI_COMM_WORLD, group2, tagcomm2, comm2, ierr)
    call errorout(ierr, 'cannot create communicator 2' )

!    if (rank < sz-2) then
!      color = 0
!    else
!      color = 1
!    endif

!    call MPI_Comm_split ( MPI_COMM_WORLD, color, rank, newComm, ierr )
!    !    print *, "rank, color :", rank, color, " ierr:", ierr
!    call errorout(ierr, 'did not split communicators' )
!
!    call MPI_COMM_RANK(newComm, new_rank, ierr)
!    call errorout(ierr, 'did not get local rank' )
!
!    if (new_rank == 0) print *, "global rank ", rank, " of ", sz, " color ", color, " new rank :" , new_rank

    ierr = iMOAB_InitializeFortran()
    call errorout(ierr, 'did not initialize fortran' )
    if (rank == 0) print *, "initialize fortran"

    if (rank >= sz-2) then
       appname='phis2'
       ierr = iMOAB_RegisterFortranApplication(appname, comm1, pid)
       print *, ' register ', appname, " on rank ", rank, " pid ", pid
    else
       appname = 'phis1'
       ierr = iMOAB_RegisterFortranApplication(appname, comm2, pid)
       print *, ' register ', appname, " on rank ", rank, " pid ", pid
    endif


    if (rank >= sz-2) then
       filename = 'spherecube.h5m'//CHAR(0)
       readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS'//CHAR(0)
       if (rank .eq. sz-2 ) print *, "loading " , trim(filename) , " with options " , trim(readopts)
       nghlay = 0

       ierr = iMOAB_LoadMesh(pid, trim(filename), trim(readopts), nghlay)
       if (rank .eq. sz-2 ) print *, "loaded in parallel ", trim(filename), " error: ", ierr
       ierr = iMOAB_SendElements(pid, comm1, MPI_COMM_WORLD, group2, pid); ! it should be different pid
       call errorout(ierr, 'cannot send elements' )
    else
       ierr = iMOAB_ReceiveElements(pid, comm2, MPI_COMM_WORLD, group1, pid); ! it should be different pid
       call errorout(ierr, 'cannot receive elements' )
       outfile = 'receivedMesh.h5m'//CHAR(0)
       wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
!      write out the mesh file to disk
       ierr = iMOAB_WriteMesh(pid, outfile, wopts)
       call errorout(ierr, 'cannot write received mesh' )
    endif

    if (MPI_COMM_NULL /= comm1) call MPI_Comm_free(comm1, ierr)
    call errorout(ierr, 'did not free comm1' )

    if (MPI_COMM_NULL /= comm2) call MPI_Comm_free(comm2, ierr)
    call errorout(ierr, 'did not free comm2' )

    call MPI_Group_free(allgroup, ierr)
    call MPI_Group_free(group1, ierr)
    call MPI_Group_free(group2, ierr)


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
