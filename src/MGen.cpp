/*
 * MGen.cpp
 *
 */

#include "moab/MGen.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

namespace moab {

MGen::MGen(Interface *impl, ParallelComm *comm, EntityHandle rset)
    : mbImpl(impl), pc(comm), cset(rset)
  {
    //ErrorCode error;

#ifdef MOAB_HAVE_MPI
    // Get the Parallel Comm instance to prepare all new sets to work in parallel
    // in case the user did not provide any arguments
    if (!comm)
      pc = moab::ParallelComm::get_pcomm(mbImpl, 0);
#endif

  }

ErrorCode MGen::BrickInstance(brOpts & opts)
{
  return MB_SUCCESS;
}
MGen::~MGen() {
  // TODO Auto-generated destructor stub
}

} /* namespace moab */
