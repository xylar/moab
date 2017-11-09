
program MigrateMesh
  implicit none

#include "moab/MOABConfig.h"
#ifdef MOAB_HAVE_MPI
#  include "mpif.h"
#else
#  error "enable parallel build"
#endif

!#define NONOVERLAP

    ! init the parallel partition
    integer ierr, sz, rank, i
    integer  newComm
    integer gcomm, comm1, comm2
    integer pid1, pid2 !  this is for physics ids
    integer compid1, compid2  !  component ids are unique over all pes, and established in
                              !  advance;
    integer nghlay ! number of ghost layers for loading
    integer groupTasks(9)   !   run on at most 9 processes
    integer startG1, startG2, endG1, endG2    !   start and end for group tasks, for creation
    integer sizeG1, sizeG2           !   size of the group that gets created
    character*10 appname
    character*132 readopts
    character*132 filename
    character*132 outfile
    character*132 wopts
    integer allgroup, group1, group2 ! Corresponding to MPI_Group in C
    integer tagcomm1, tagcomm2
    integer iMOAB_InitializeFortran, iMOAB_RegisterFortranApplication
    integer iMOAB_LoadMesh, iMOAB_SendMesh, iMOAB_ReceiveMesh, iMOAB_WriteMesh
    integer iMOAB_FreeSenderBuffers
    integer iMOAB_DeregisterApplication, iMOAB_Finalize

    call MPI_INIT(ierr)
    call MPI_Comm_dup(MPI_COMM_WORLD, gcomm, ierr)
    call MPI_COMM_SIZE(gcomm, sz, ierr)
    call MPI_COMM_RANK(gcomm, rank, ierr)
    if (rank .eq. 0) print *, "size:", sz
    call errorout(ierr, 'cannot get rank' )
    if ( (0 .eq. rank) .and. (sz < 2) .and. (sz>9) ) then
      print *, "size is " , sz, ". run on at least 2 processes and at most 9 "
      call exit(1)
    endif
    ! create 2 overlapping groups, for generality
    ! create communicators for each group; 
    ! one group will represent the sender, the other group the receiver
    ! about one third of tasks will be on group 1 only,  and one fourth will be on group 2 only
    !  about (1-1./3 -1./4) will be overlapping, these tasks will be common to both groups
    ! the mesh will be read on the sender comm, will be sent to receiver comm 

    ! create new MPI groups for processors 1/3*sz, 1/3*sz+1, ..., sz-1 (group 1) and 0, 1, .., 3/4*sz-1 (group 2)
    

    call MPI_COMM_GROUP (gcomm, allgroup, ierr)
    call errorout(ierr, 'cannot get world group' )
    ! first group, sz/3 to sz-1
    startG1 = sz/2
    endG1 = sz-1
    sizeG1 = endG1 - startG1 + 1
    
    do i=1, sizeG1
      groupTasks (i) = startG1+i-1
    end do 

    call MPI_Group_incl(allgroup, sizeG1, groupTasks, group1, ierr)
    call errorout(ierr, 'cannot create group 1' )

    ! second group, 0, 1, 3/4*sz
    startG2 = 0
    endG2 = 3*sz/4 -1
#ifdef NONOVERLAP
    endG2 = startG1-1
#endif
    sizeG2 = endG2 - startG2 + 1
    do i=1, sizeG2
      groupTasks(i) = startG2+i-1
    enddo 
    
    call MPI_Group_incl(allgroup, sizeG2, groupTasks, group2, ierr)
    call errorout(ierr, 'cannot create group 2' )

    if ( (0 .eq. rank) ) then
      print *, "group 1 tasks: ", (i, i=startG1, endG1)
      print *, "group 2 tasks: ", (i, i=startG2, endG2)
    endif
    ! now create both communicators
    !  when we are not on tasks in the communicator, the MPI_Comm created will be null
    tagcomm1 = 1
    call MPI_Comm_create_group(gcomm, group1, tagcomm1, comm1, ierr)
    call errorout(ierr, 'cannot create communicator 1' )

    tagcomm2 = 2
    call MPI_Comm_create_group(gcomm, group2, tagcomm2, comm2, ierr)
    call errorout(ierr, 'cannot create communicator 2' )


    ierr = iMOAB_InitializeFortran()

    ! give some dummy values to component ids, just to differentiate between them
    ! the par comm graph is unique between components
    compid1 = 4
    compid2 = 7
    call errorout(ierr, 'did not initialize fortran' )
    if (rank == 0) print *, "initialize iMOAB fortran applications"

    if (comm1 /= MPI_COMM_NULL) then
       appname='phis1'//CHAR(0)
       ierr = iMOAB_RegisterFortranApplication(trim(appname), comm1, compid1, pid1)
       print *, ' register ', appname, " on rank ", rank, " pid1 ", pid1
    endif
    if (comm2 /= MPI_COMM_NULL) then
       appname = 'phis2'//CHAR(0)
       ierr = iMOAB_RegisterFortranApplication(trim(appname), comm2, compid2, pid2)
       print *, ' register ', appname, " on rank ", rank, " pid2 ", pid2
    endif
    

    if (comm1 /= MPI_COMM_NULL) then
       filename = 'spherecube.h5m'//CHAR(0)
       readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS'//CHAR(0)
       if (rank .eq. sz-2 ) print *, "loading " , trim(filename) , " with options " , trim(readopts)
       nghlay = 0

       ierr = iMOAB_LoadMesh(pid1, trim(filename), trim(readopts), nghlay)
       if (rank .eq. sz-2 ) print *, "loaded in parallel ", trim(filename), " error: ", ierr
       ierr = iMOAB_SendMesh(pid1, gcomm, group2, compid2); ! send to component 2
       call errorout(ierr, 'cannot send elements' )
    endif

    if (comm2 /= MPI_COMM_NULL) then
       ierr = iMOAB_ReceiveMesh(pid2, gcomm, group1, compid1); ! receive from component 1
       call errorout(ierr, 'cannot receive elements' )
       outfile = 'receivedMesh.h5m'//CHAR(0)
       if (pid2 .eq. 0) then
         wopts   = 'PARALLEL=WRITE_PART;PARALLEL_COMM=0'//CHAR(0)
       else if (pid2 .eq. 1) then
         wopts   = 'PARALLEL=WRITE_PART;PARALLEL_COMM=1'//CHAR(0)
       else
          call errorout(1, "wrong pid2")
       endif

       print *, "from ", rank, wopts, outfile
!      write out the mesh file to disk
       ierr = iMOAB_WriteMesh(pid2, trim(outfile), trim(wopts))
       call errorout(ierr, 'cannot write received mesh' )
    endif

    call MPI_Barrier(gcomm, ierr)
    call errorout(ierr, 'cannot stop at barrier' )

    ! we can now free the sender buffers
    if (comm1 /= MPI_COMM_NULL) then
       ierr = iMOAB_FreeSenderBuffers(pid1, gcomm, compid2)
    endif

    if (comm1 /= MPI_COMM_NULL) then
       ierr = iMOAB_DeregisterApplication(pid1)
         call errorout(ierr, 'cannot deregister app 1 sender' )
    endif
    if (comm2 /= MPI_COMM_NULL) then
       ierr = iMOAB_DeregisterApplication(pid2)
       call errorout(ierr, 'cannot deregister app 2 receiver' )
    endif

    ierr = iMOAB_Finalize()
    call errorout(ierr, 'did not finalize iMOAB' )

    if (MPI_COMM_NULL /= comm1) call MPI_Comm_free(comm1, ierr)
    call errorout(ierr, 'did not free comm1' )

    if (MPI_COMM_NULL /= comm2) call MPI_Comm_free(comm2, ierr)
    call errorout(ierr, 'did not free comm2' )

    call MPI_Group_free(allgroup, ierr)
    call MPI_Group_free(group1, ierr)
    call MPI_Group_free(group2, ierr)
    call MPI_Comm_free(gcomm, ierr)

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
