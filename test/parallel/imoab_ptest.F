C      push parallel mesh into moab, using iMoab API, F77 version
C
C      This program shows how to push a mesh into MOAB in parallel using iMoab, with sufficient
C      information to resolve boundary sharing and exchange a layer of ghost information.
C
C       After resolving the sharing, mesh is saved to a file, then read back again, in parallel
C        the mesh is a quad mesh, 2x2 on each task,
C       by default, this test is run on 2 processors

C       Test is similar to itaps/iMeshP_unit_tests.cpp
C
C
C       Create a mesh:
C       Groups of four quads will be arranged into parts as follows:
C       +------+------+------+------+------+-----
C       |             |             |
C       |             |             |
C       +    Part 0   +    Part 2   +    Part 4
C       |             |             |
C       |             |             |
C       +------+------+------+------+------+-----
C       |             |             |
C       |             |             |
C       +    Part 1   +    Part 3   +    Part 5
C       |             |             |
C       |             |             |
C       +------+------+------+------+------+-----
C
C       Vertices will be enumerated as follows:
C       1------6-----11-----16-----21-----26----- > x  x from 0, 1, 2
C       |             |             |
C       |             |             |
C       2      7     12     17     22     27
C       |             |             |
C       |             |             |
C       3------8-----13-----18-----23-----28-----
C       |             |             |
C       |             |             |
C       4      9     14     19     24     29
C       |             |             |
C       |             |             |
C       5-----10-----15-----20-----25-----30-----
C       |
C       y varies from 0 to 4
C
C       Processor 0 will have vertices 1, 2, 3, 6, 7, 8, 11, 12, 13
C         and 4 quads, with ids from 1 to 4, and connectivity
C          1, 2, 7, 6;  2, 3, 8, 7; 6 ,7, 12, 11; 7, 8, 13, 12
C       Processor 1 will have vertices 3, 4, 5, 8, 9, 10, 13, 14, 15
C         and 4 quads, with ids from 5 to 8, and connectivity
C          3, 4, 8, 9;  4, 5, 10, 9; 8, 9, 14, 13; 9, 10, 15, 14,
C       and so on
C
C       Vertex Global IDs will be used to resolve sharing
C
C       Element IDs will be [4*rank+1,4*rank+5]
C  #define ERROR(rval) if (0 .ne. rval) call exit(1)


      SUBROUTINE errorout(ierr, message)
      integer ierr
      character*(*) message
      if (ierr.ne.0) then
        print *, message
        call exit (1)
      end if
      return
      end

#include "moab/MOABConfig.h"
      program imoabpar_test
C      implicit none

#include "mpif.h"

      integer ierr, my_id, num_procs, pid, i, ix, iy, numv, nume
      integer dime, lco, mbtype, blockid, npe
      character*10 appname
      character*100 outfile, wopts
C coordinates for 9 vertices
      real*8  coordinates(27) , coords_core(27), deltax, deltay
      integer ids(9)
      integer connec(16)
C      used for ghosting
      integer dimgh, bridge, num_layers

      data coords_core/ 0., 0., 0.,
     &                  0., 1., 0.,
     &                  0., 2., 0.,
     &                  1., 0., 0.,
     &                  1., 1., 0.,
     &                  1., 2., 0.,
     &                  2., 0., 0.,
     &                  2., 1., 0.,
     &                  2., 2., 0./
      data ids/1, 2, 3, 6, 7, 8, 11, 12, 13/
      data connec/ 1, 2, 5, 4,
     &             2, 3, 6, 5,
     &             4, 5, 8, 7,
     &             5, 6, 9, 8/


      call MPI_INIT ( ierr )
      call errorout(ierr, 'fail to initialize MPI')

c     find out MY process ID, and how many processes were started.

      call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
      call errorout(ierr, 'fail to get MPI rank')

      call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
      call errorout(ierr, 'fail to get MPI size')

      if (my_id .eq. 0) then
        print *, " I'm process ", my_id, " out of ",
     &      num_procs, " processes."

      end if

      ierr = iMOAB_InitializeFortran()
      call errorout(ierr, 'fail to initialize iMOAB')

      appname = 'IMTEST'//CHAR(0)
      ierr = iMOAB_RegisterFortranApplication(appname, MPI_COMM_WORLD,
     & pid)
      call errorout(ierr, 'fail to initialize fortran app')

C  create first 9 vertices,  in a square 3x3;
      deltax=(my_id/2) * 2.
      deltay = mod(my_id, 2) * 2.
      ix = (my_id/2) * 10
      iy = mod(my_id,2) * 2
      do i=1,9
        coordinates( 3*(i-1) + 1 ) = coords_core (3*(i-1)+1) + deltax
        coordinates( 3*(i-1) + 2 ) = coords_core (3*(i-1)+2) + deltay
        coordinates( 3*(i-1) + 3 ) = coords_core (3*(i-1)+3)
C        translate the ids too, by multiples of 10 or add 2
        ids(i) = ids(i) + ix+iy
      enddo
      numv = 9
      nume = 4
      lco = numv*3
      dime = 3
      ierr = iMOAB_CreateVertices(pid, lco, dime, coordinates);
      call errorout(ierr, 'fail to create vertices')
C create now 4 quads with those 9 vertices
      mbtype = 3
      blockid  =100
      npe = 4
      ierr = iMOAB_CreateElements(pid, nume, mbtype, npe, connec,
     & blockid)
      call errorout(ierr, 'fail to create quads')

      ierr = iMOAB_ResolveSharedEntities( pid, numv, ids );
      call errorout(ierr, 'fail to resolve shared entities')

C  see ghost elements
      dimgh = 2 ! will ghost quads, topological dim 2
      bridge = 0 ! use vertex as bridge
      num_layers = 1 ! so far, one layer only
      ierr = iMOAB_DetermineGhostEntities( pid, dimgh, num_layers,
     & bridge)
      call errorout(ierr, 'fail to determine ghosts')
C     write out the mesh file to disk, in parallel, if h5m
#ifdef MOAB_HAVE_HDF5
      outfile = 'whole.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(pid, outfile, wopts)
      call errorout(ierr, 'fail to write the mesh file')
#endif
C      all done. de-register and finalize
      ierr = iMOAB_DeregisterApplication(pid)
      call errorout(ierr, 'fail to deregister application')

      ierr = iMOAB_Finalize()
      call errorout(ierr, 'fail to finalize iMOAB')


      call MPI_FINALIZE ( ierr )
      call errorout(ierr, 'fail to finalize MPI')

      end

