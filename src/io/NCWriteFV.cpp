/*
 * NCWriteFV.cpp
 *
 *  Created on: April 9, 2014
 */

#include "NCWriteFV.hpp"
#include "moab/WriteUtilIface.hpp"

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { _writeNC->mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { _writeNC->mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

NCWriteFV::~NCWriteFV()
{
  // TODO Auto-generated destructor stub
}

ErrorCode NCWriteFV::write_values(std::vector<std::string>& var_names, EntityHandle fileSet)
{
  return MB_NOT_IMPLEMENTED;
}

} /* namespace moab */
