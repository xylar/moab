      program fdriver
C      include '/home/iulian/3rdparty/mpich-3.1/install/include/mpif.h'
      include 'mpif.h'
      integer ierr, num_procs, my_id
      integer pid
      character*10 appname

      call MPI_INIT ( ierr )

c     find out MY process ID, and how many processes were started.

      call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

      if (my_id .eq. 0) then

        print *, " I'm process ", my_id, " out of ",
     &      num_procs, " processes."

      end if

      ierr = iMOABInitializeFortran()

      appname = 'PROTEUS'

      ierr = RegisterFortranApplication(appname, MPI_COMM_WORLD, pid)

      ierr = DeregisterApplication(pid)

      ierr = iMOABFinalize()


      call MPI_FINALIZE ( ierr )

      stop
      end   
