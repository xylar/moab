#ifndef DAGMC_PREPROC_H
#define DAGMC_PREPROC_H

#include <iostream>
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"
using namespace moab;


/// program-wide verbose output flag
extern bool verbose;

/**
 * Generic halt-on-error error checking
 */
void chkerr( Interface* mbi, ErrorCode code, int line, const char* file );
void chkerr( Interface& mbi, ErrorCode code, int line, const char* file );
void chkerr( GeomTopoTool& gtt, ErrorCode code, int line, const char* file );

#define CHECKERR(M,C) chkerr(M,C,__LINE__,__FILE__)

// features provided by obb_analysis.cpp
/* no longer provided by dagmc
ErrorCode obbvis_create( DagMC& dag, std::vector<int> &volumes, int grid, std::string& filename );
ErrorCode obbstat_write( DagMC& dag, std::vector<int> &volumes,
                         std::vector<std::string> &properties, std::ostream& out );
*/

#endif /* DAGMC_PREPROC_H */
