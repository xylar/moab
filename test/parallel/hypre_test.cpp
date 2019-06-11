// Parts of this test have been taken from MFEM and HPCG benchmark example

#include "moab/Core.hpp"
#include "HypreSolver.hpp"

#include <iostream>
using namespace std;

#undef DEBUG

moab::ErrorCode GenerateTestMatrixAndVectors ( int nx, int ny, int nz,
    moab::HypreParMatrix &A,
    moab::HypreParVector &sol,
    moab::HypreParVector &rhs,
    moab::HypreParVector &exactsol );

int main ( int argc, char *argv[] )
{
  // double norm, d;
  int ierr = 0;
  double times[5];
  int nx=5, ny=5, nz=5;
#ifdef MOAB_HAVE_MPI
  MPI_Init ( &argc, &argv );
  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size ( MPI_COMM_WORLD, &size );
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
#endif
  MPI_Comm comm = MPI_COMM_WORLD;

  moab::Core mbCore;
  moab::Interface& mb = mbCore;

  // rval = mb.load_file(example.c_str(), &euler_set, opts.c_str());

  moab::ParallelComm* pcomm = new moab::ParallelComm(&mb, comm);

  moab::HypreParMatrix A ( pcomm );
  moab::HypreParVector x ( pcomm ), b ( pcomm ), xexact ( pcomm );
#ifdef DEBUG

  if ( rank == 0 ) {
    int junk = 0;
    cout << "Press enter to continue" << endl;
    cin >> junk;
  }

  MPI_Barrier ( MPI_COMM_WORLD );
#endif

  if ( argc > 4 ) {
    if ( rank == 0 )
      cerr << "Usage:" << endl
           << "\t" << argv[0] << " nx ny nz" << endl
           << "     where nx, ny and nz are the local sub-block dimensions, or" << endl;

    exit ( 1 );
  }
  else {
    if ( rank == 0 )
      cout << "Using command:" << "\t" << argv[0] << " nx ny nz" << endl;
  }

  times[0] = MPI_Wtime();   // Initialize it (if needed)

  if ( argc == 4 ) {
    nx = atoi ( argv[1] );
    ny = atoi ( argv[2] );
    nz = atoi ( argv[3] );
  }

  GenerateTestMatrixAndVectors ( nx, ny, nz, A, x, b, xexact );

  times[1] = MPI_Wtime() - times[0];   // operator creation time
  const bool useCG = true; // TODO make into command line option
  const bool useAMG = false; // TODO make into command line option
  int niters = 0;
  double normr = 0;
  int max_iter = 150;
  double tolerance = 1e-15; // Set tolerance to zero to make all runs do max_iter iterations
  times[2] = MPI_Wtime();   // reset
  moab::HypreSolver *lsSolver, *psSolver;

  if ( useCG ) {
    lsSolver = new moab::HyprePCG ( A );
    moab::HyprePCG *cgSolver = dynamic_cast<moab::HyprePCG *> ( lsSolver );
    cgSolver->SetTol ( tolerance );
    cgSolver->SetMaxIter ( max_iter );
    cgSolver->SetLogging ( 1 );
    cgSolver->Verbosity ( 1 );
    cgSolver->SetResidualConvergenceOptions ( 3, 0.0 );
  }
  else {
    lsSolver = new moab::HypreGMRES ( A );
    moab::HypreGMRES *gmresSolver = dynamic_cast<moab::HypreGMRES *> ( lsSolver );
    gmresSolver->SetTol ( tolerance, normr );
    gmresSolver->SetMaxIter ( max_iter );
    gmresSolver->SetLogging ( 1 );
    gmresSolver->Verbosity ( 5 );
    gmresSolver->SetKDim ( 30 );
  }

  if ( useAMG ) {
    psSolver = new moab::HypreBoomerAMG ( A );
    dynamic_cast<moab::HypreBoomerAMG *> ( psSolver )->SetSystemsOptions ( 3 );

  } else {
    psSolver = new moab::HypreParaSails ( A );
    dynamic_cast<moab::HypreParaSails *> ( psSolver )->SetSymmetry ( 1 );
  }

  if ( psSolver )
    lsSolver->SetPreconditioner ( *psSolver );

  times[2] = MPI_Wtime() - times[2];   // solver setup time
  /* Now let us solve the linear system */
  times[3] = MPI_Wtime();   // reset
  ierr = lsSolver->Solve ( b, x );

  if ( ierr ) cerr << "Error in call to CG/GMRES HYPRE Solver: " << ierr << ".\n" << endl;

  times[3] = MPI_Wtime() - times[3];   // solver time
  lsSolver->GetNumIterations ( niters );
  lsSolver->GetFinalResidualNorm ( normr );
  times[4] = MPI_Wtime() - times[0];
  // x.Print ( "solution.out" );
#ifdef MOAB_HAVE_MPI
  double t4 = times[3];
  double t4min = 0.0;
  double t4max = 0.0;
  double t4avg = 0.0;
  MPI_Allreduce ( &t4, &t4min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  MPI_Allreduce ( &t4, &t4max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  MPI_Allreduce ( &t4, &t4avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  t4avg = t4avg / ( ( double ) size );
#endif

  if ( rank == 0 ) { // Only PE 0 needs to compute and report timing results
    // double fniters = niters;
    // double fnrow = A.NNZ();
    // double fnnz = A.M();
    // double fnops_ddot = fniters * 4 * fnrow;
    // double fnops_waxpby = fniters * 6 * fnrow;
    // double fnops_sparsemv = fniters * 2 * fnnz;
    // double fnops = fnops_ddot + fnops_waxpby + fnops_sparsemv;
    std::cout << "Testing Laplace problem solver with HYPRE interface\n";
    {
      std::cout << "\tParallelism\n";
#ifdef MOAB_HAVE_MPI
      std::cout << "Number of MPI ranks: \t" << size << std::endl;
#else
      std::cout << "MPI not enabled\n";
#endif
    }
    std::cout << "Dimensions (nx*ny*nz): \t(" << nx << "*" << ny << "*" << nz << ")\n";
    std::cout << "Number of iterations: \t" << niters << std::endl;
    std::cout << "Final residual: \t" << normr << std::endl;
    std::cout << "************ Performance Summary (times in sec) *************" << std::endl;
    {
      std::cout << "\nTime Summary" << std::endl;
      std::cout << "Total             \t" << times[4] << std::endl;
      std::cout << "Operator Creation \t" << times[1] << std::endl;
      std::cout << "Solver Setup      \t" << times[2] << std::endl;
      std::cout << "HYPRE Solver      \t" << times[3] << std::endl;
      std::cout << "Avg. Solver       \t" << t4avg << std::endl;
    }
  }

  // Compute difference between known exact solution and computed solution
  // All processors are needed here.
//   double residual = 0;
//   if ((ierr = ComputeLinfError(x, xexact, &residual)))
//     cerr << "Error in call to ComputeLinfError: " << ierr << ".\n" << endl;
//   if (rank==0)
//      cout << "Difference between computed and exact  = "
//           << residual << ".\n" << endl;
  delete lsSolver;
  // Finish up
#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  return 0 ;
}


moab::ErrorCode GenerateTestMatrixAndVectors(int nx, int ny, int nz,
    moab::HypreParMatrix &A,
    moab::HypreParVector &sol,
    moab::HypreParVector &rhs,
    moab::HypreParVector &exactsol)
{
#ifdef DEBUG
  int debug = 1;
#else
  int debug = 0;
#endif
#ifdef MOAB_HAVE_MPI
  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
#endif
  // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
  bool use_7pt_stencil = false;
  const int NNZPERROW = 27;

  int local_nrow = nx * ny * nz; // This is the size of our subblock
  assert(local_nrow > 0);    // Must have something to work with
  int local_nnz = NNZPERROW * local_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)
  int total_nrow = 0;
  MPI_Allreduce(&local_nrow, &total_nrow, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // int total_nrow = local_nrow * size; // Total number of grid points in mesh
  long long total_nnz = NNZPERROW * (long long)
                        total_nrow;   // Approximately 27 nonzeros per row (except for boundary nodes)
  int start_row = local_nrow * rank; // Each processor gets a section of a chimney stack domain
  int stop_row = start_row + local_nrow - 1;

  // Allocate arrays that are of length local_nrow
  // Eigen::SparseMatrix<double, Eigen::RowMajor> diagMatrix(local_nrow, local_nrow);
  // diagMatrix.reserve(local_nnz);
  HYPRE_Int col[2] = {start_row, stop_row};
  A.resize(total_nrow, col);

  sol.resize(total_nrow, start_row, stop_row);
  rhs.resize(total_nrow, start_row, stop_row);
  exactsol.resize(total_nrow, start_row, stop_row);
  long long nnzglobal = 0;

  for (int iz = 0; iz < nz; iz++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        const int curlocalrow = iz * nx * ny + iy * nx + ix;
        const int currow = start_row + curlocalrow;
        int nnzrow = 0;
        std::vector<double> colvals;
        std::vector<int> indices;

        for (int sz = -1; sz <= 1; sz++) {
          for (int sy = -1; sy <= 1; sy++) {
            for (int sx = -1; sx <= 1; sx++) {
              const int curlocalcol = sz * nx * ny + sy * nx + sx;
              const int curcol = currow + curlocalcol;

              // Since we have a stack of nx by ny by nz domains , stacking in the z direction, we check to see
              // if sx and sy are reaching outside of the domain, while the check for the curcol being valid
              // is sufficient to check the z values
              if ((ix + sx >= 0) && (ix + sx < nx) && (iy + sy >= 0) && (iy + sy < ny) && (curcol >= 0 && curcol < total_nrow)) {
                if (!use_7pt_stencil ||
                    (sz * sz + sy * sy + sx * sx <= 1)) {     // This logic will skip over point that are not part of a 7-pt stencil
                  if (curcol == currow) {
                    colvals.push_back(27.0);
                  } else {
                    colvals.push_back(-1.0);
                  }

                  indices.push_back(curcol);
                  nnzrow++;
                }
              }
            } // end sx loop
          } // end sy loop
        } // end sz loop

        int ncols = indices.size();
        A.SetValues(1, &ncols, &currow, indices.data(), colvals.data());
        // diagMatrix.set_row_values ( curlocalrow, indices.size(), indices.data(), colvals.data() );
        nnzglobal += nnzrow;
        sol.SetValue(currow, 0.0);
        rhs.SetValue(currow, 27.0 - ((double)(nnzrow - 1)));
        exactsol.SetValue(currow, 1.0);
      } // end ix loop
    } // end iy loop
  } // end iz loop

  if (debug) cout << "Global size of the matrix: " << total_nrow << " and NNZ (estimate) = " << total_nnz << " NNZ (actual) = " << nnzglobal <<
                    endl;

  if (debug) cout << "Process " << rank << " of " << size << " has " << local_nrow << endl;

  if (debug) cout << " rows. Global rows " << start_row
                    << " through " << stop_row << endl;

  if (debug) cout << "Process " << rank << " of " << size
                    << " has " << local_nnz << " nonzeros." << endl;

  A.FinalizeAssembly();
  sol.FinalizeAssembly();
  rhs.FinalizeAssembly();
  exactsol.FinalizeAssembly();
  // A.Print("matrix.out");
  return moab::MB_SUCCESS;
}
