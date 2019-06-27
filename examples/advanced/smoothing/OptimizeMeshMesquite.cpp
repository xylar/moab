
#include "Mesquite.hpp"
#include "MsqIBase.hpp"
#include "MsqIGeom.hpp"
#include "MsqIMesh.hpp"
#include "MBiMesh.hpp"
#include "MeshImpl.hpp"
#include "moab/Core.hpp"
#include "moab/Skinner.hpp"
#include "moab/LloydSmoother.hpp"
#include "FacetModifyEngine.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "SphericalDomain.hpp"
#include "PlanarDomain.hpp"
#include "MeshWriter.hpp"

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

// algorithms
#include "IdealWeightInverseMeanRatio.hpp"
#include "TMPQualityMetric.hpp"
#include "AspectRatioGammaQualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "VertexConditionNumberQualityMetric.hpp"
#include "MultiplyQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "PMeanPTemplate.hpp"
#include "SteepestDescent.hpp"
#include "FeasibleNewton.hpp"
#include "QuasiNewton.hpp"
#include "ConjugateGradient.hpp"
#include "SmartLaplacianSmoother.hpp"
#include "Randomize.hpp"

#include "ReferenceMesh.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "TShapeB1.hpp"
#include "TQualityMetric.hpp"
#include "IdealShapeTarget.hpp"

#include <sys/stat.h>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include "iBase.h"

using namespace MBMesquite;

static int print_iGeom_error( const char* desc,
                         int err,
                         iGeom_Instance geom,
                         const char* file,
                         int line )
{
  char buffer[1024];
  iGeom_getDescription( geom, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';

  std::cerr << "ERROR: " << desc << std::endl
            << "  Error code: " << err << std::endl
            << "  Error desc: " << buffer << std::endl
            << "  At        : " << file << ':' << line << std::endl
            ;

  return -1; // must always return false or CHECK macro will break
}

static int print_iMesh_error( const char* desc,
                         int err,
                         iMesh_Instance mesh,
                         const char* file,
                         int line )
{
  char buffer[1024];
  iMesh_getDescription( mesh, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';

  std::cerr << "ERROR: " << desc << std::endl
            << "  Error code: " << err << std::endl
            << "  Error desc: " << buffer << std::endl
            << "  At        : " << file << ':' << line << std::endl
            ;

  return -1; // must always return false or CHECK macro will break
}

#define CHECK_IGEOM( STR ) if (err != iBase_SUCCESS) return print_iGeom_error( STR, err, geom, __FILE__, __LINE__ )

#define CHECK_IMESH( STR ) if (err != iBase_SUCCESS) return print_iMesh_error( STR, err, instance, __FILE__, __LINE__ )

const std::string default_file_name = std::string(MESH_DIR) + std::string("/surfrandomtris-4part.h5m");
const std::string default_geometry_file_name = std::string(MESH_DIR) + std::string("/surfrandom.facet");

std::vector<double> solution_indicator;

int write_vtk_mesh( Mesh* mesh, const char* filename);

// Construct a MeshTSTT from the file
int get_imesh_mesh( MBMesquite::Mesh** , const char* file_name, int dimension );

  // Construct a MeshImpl from the file
int get_native_mesh( MBMesquite::Mesh** , const char* file_name, int dimension );

int get_itaps_domain(MeshDomain**, const char*);
int get_mesquite_domain(MeshDomain**, const char*);

  // Run FeasibleNewton solver
int run_global_smoother( MeshDomainAssoc& mesh, MsqError& err, double OF_value=1e-4 );

  // Run SmoothLaplacian solver
int run_local_smoother( MeshDomainAssoc& mesh, MsqError& err, double OF_value=1e-3 );
int run_local_smoother2( MeshDomainAssoc& mesh_and_domain, MsqError& err, double OF_value=1e-3 );

int run_quality_optimizer( MeshDomainAssoc& mesh_and_domain, MsqError& err );

int run_solution_mesh_optimizer( MeshDomainAssoc& mesh_and_domain, MsqError& err );

bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[])
{
  MBMesquite::MsqPrintError err(cout);
  // command line arguments
  std::string file_name, geometry_file_name;
  bool use_native = false;
  int dimension = 2;

#ifdef MOAB_HAVE_MPI
//  MPI_Init(&argc, &argv);
#endif
  ProgOptions opts;

  opts.addOpt<void>(std::string("native,N"),
      std::string("Use native representation (default=false)"), &use_native);
  opts.addOpt<int>(std::string("dim,d"),
      std::string("Topological dimension of the mesh (default=2)"), &dimension);
  opts.addOpt<std::string>(std::string("input_geo,i"),
      std::string("Input file name for the geometry (pattern=*.stl, *.facet)"), &geometry_file_name);
  opts.addOpt<std::string>(std::string("input_mesh,m"),
      std::string("Input file name for the mesh (pattern=*.vtk, *.h5m)"), &file_name);

  opts.parseCommandLine(argc, argv);

  if (!geometry_file_name.length())
  {
    file_name = default_file_name;
    geometry_file_name = default_geometry_file_name;
    cout << "No file specified: Using defaults.\n";
  }
  cout << "\t Mesh filename = " << file_name << endl;
  cout << "\t Geometry filename = " << geometry_file_name << endl;

  // Vector3D pnt(0,0,0);
  // Vector3D s_norm(0,0,1);
  // PlanarDomain plane(s_norm, pnt);

//  PlanarDomain plane( PlanarDomain::XY );
  MeshDomain* plane;
  int ierr;
  if (!file_exists(geometry_file_name)) geometry_file_name = "";
  ierr = get_itaps_domain(&plane, geometry_file_name.c_str());//MB_CHK_ERR(ierr);

    // Try running a global smoother on the mesh
  Mesh* mesh=0;
  if(use_native) {
    ierr = get_native_mesh(&mesh, file_name.c_str(), dimension);//MB_CHK_ERR(ierr);
  }
  else {
    ierr = get_imesh_mesh(&mesh, file_name.c_str(), dimension);//MB_CHK_ERR(ierr);
  }

  if (!mesh) {
    std::cerr << "Failed to load input file.  Aborting." << std::endl;
    return 1;
  }

  MeshDomainAssoc mesh_and_domain(mesh, plane);

  // ierr = write_vtk_mesh( mesh, "original.vtk");MB_CHK_ERR(ierr);
  // cout << "Wrote \"original.vtk\"" << endl;

  // run_global_smoother( mesh_and_domain, err );
  // if (err) return 1;

    // Try running a local smoother on the mesh
//  Mesh* meshl=0;
//  if(use_native)
//    ierr = get_native_mesh(&meshl, file_name.c_str(), dimension);
//  else
//    ierr = get_imesh_mesh(&meshl, file_name.c_str(), dimension);
//  if (!mesh || ierr) {
//    std::cerr << "Failed to load input file.  Aborting." << std::endl;
//    return 1;
//  }

//  MeshDomainAssoc mesh_and_domain_local(meshl, plane);

  // run_solution_mesh_optimizer( mesh_and_domain, err );
  // if (err) return 1;

  run_local_smoother( mesh_and_domain, err, 1e-4 );//MB_CHK_ERR(err);
  if (err) return 1;

  double reps = 5e-2;
  for (int iter=0; iter < 5; iter++) {

    if (!(iter % 2)) {
      run_local_smoother2( mesh_and_domain, err, reps*10 );//CHECK_IMESH("local smoother2 failed");
      if (err) return 1;
    }

    // run_global_smoother( mesh_and_domain, err, reps );MB_CHK_ERR(ierr);

    run_solution_mesh_optimizer( mesh_and_domain, err );//CHECK_IMESH("solution mesh optimizer failed");
    if (err) return 1;

    reps *= 0.01;
  }

  run_local_smoother2( mesh_and_domain, err, 1e-4 );//CHECK_IMESH("local smoother2 failed");
  if (err) return 1;

  // run_quality_optimizer( mesh_and_domain, err );MB_CHK_ERR(ierr);

  // run_local_smoother( mesh_and_domain, err );MB_CHK_ERR(ierr);

  // Delete MOAB instance
  delete mesh;
  delete plane;

#ifdef MOAB_HAVE_MPI
  //MPI_Finalize();
#endif

  return 0;
}


int run_global_smoother( MeshDomainAssoc& mesh_and_domain, MsqError& err, double OF_value )
{
  // double OF_value = 1e-6;

  Mesh* mesh = mesh_and_domain.get_mesh();
  MeshDomain* domain = mesh_and_domain.get_domain();

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  // ConditionNumberQualityMetric* mean_ratio = new ConditionNumberQualityMetric();
  // TMPQualityMetric* mean_ratio = new TMPQualityMetric();

  // VertexConditionNumberQualityMetric* mean_ratio = new VertexConditionNumberQualityMetric();
  if (err) return 1;
  mean_ratio->set_averaging_method(QualityMetric::RMS, err);
  if (err) return 1;

  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  // LInfTemplate* obj_func = new LInfTemplate(mean_ratio);
  if (err) return 1;

  // creates the feas newt optimization procedures
  // ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  // FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  SteepestDescent* pass1 = new SteepestDescent( obj_func );
  pass1->use_element_on_vertex_patch();
  pass1->do_block_coordinate_descent_optimization();
  pass1->use_global_patch();
  if (err) return 1;

  QualityAssessor stop_qa( mean_ratio );

  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  if (err) return 1;
  TerminationCriterion tc_outer;
  tc_inner.add_iteration_limit( 10 );
  tc_outer.add_iteration_limit( 5 );
  tc_outer.set_debug_output_level(3);
  tc_inner.set_debug_output_level(3);
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // launches optimization on mesh_set
  if (domain) {
    queue1.run_instructions(&mesh_and_domain, err);
  }
  else {
    queue1.run_instructions(mesh, err);
  }
  if (err) return 1;


  // Construct a MeshTSTT from the file
  int ierr = write_vtk_mesh( mesh, "feasible-newton-result.vtk");
  if (ierr) return 1;
  // MeshWriter::write_vtk(mesh, "feasible-newton-result.vtk", err);
  // if (err) return 1;
  cout << "Wrote \"feasible-newton-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}

int write_vtk_mesh( Mesh* mesh, const char* filename )
{
  moab::Interface* mbi = reinterpret_cast<MBiMesh*>(dynamic_cast<MsqIMesh*>(mesh)->get_imesh_instance())->mbImpl;

  mbi->write_file(filename);

  return 0;
}

int run_local_smoother2( MeshDomainAssoc& mesh_and_domain, MsqError& err, double OF_value );
int run_local_smoother( MeshDomainAssoc& mesh_and_domain, MsqError& err, double OF_value )
{
  Mesh* mesh = mesh_and_domain.get_mesh();
  moab::Interface* mbi = reinterpret_cast<MBiMesh*>(dynamic_cast<MsqIMesh*>(mesh)->get_imesh_instance())->mbImpl;


  moab::Tag fixed;
  moab::ErrorCode rval = mbi->tag_get_handle("fixed", 1, moab::MB_TYPE_INTEGER, fixed); MB_CHK_SET_ERR(rval, "Getting tag handle failed");
  moab::Range cells;
  rval = mbi->get_entities_by_dimension(0, 2, cells); MB_CHK_SET_ERR(rval, "Querying elements failed");

  moab::LloydSmoother lloyd(mbi, 0, cells, 0, 0/*fixed*/);

  lloyd.perform_smooth();

  // run_local_smoother2(mesh_and_domain, err, OF_value);

  // Construct a MeshTSTT from the file
  int ierr = write_vtk_mesh( mesh, "smart-laplacian-result.vtk");
  if (ierr) return 1;
  // MeshWriter::write_vtk(mesh, "smart-laplacian-result.vtk", err);
  // if (err) return 1;
  cout << "Wrote \"smart-laplacian-result.vtk\"" << endl;
  return 0;
}

int run_local_smoother2( MeshDomainAssoc& mesh_and_domain, MsqError& err, double OF_value )
{
  // double OF_value = 1e-5;

  Mesh* mesh = mesh_and_domain.get_mesh();
  MeshDomain* domain = mesh_and_domain.get_domain();

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  // IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  ConditionNumberQualityMetric* mean_ratio = new ConditionNumberQualityMetric();
  // VertexConditionNumberQualityMetric* mean_ratio = new VertexConditionNumberQualityMetric();
  if (err) return 1;
  //  mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
  //   mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
  mean_ratio->set_averaging_method(QualityMetric::RMS, err);
  if (err) return 1;

  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;

  if (false)
  {
    InstructionQueue qtmp;
    Randomize rand (-0.005) ;
    TerminationCriterion sc_rand;
    sc_rand.add_iteration_limit( 2 ) ;
    rand.set_outer_termination_criterion(&sc_rand);
    qtmp.set_master_quality_improver(&rand, err);
    if (err) return 1;
    if (domain) {
      qtmp.run_instructions(&mesh_and_domain, err);
    }
    else {
      qtmp.run_instructions(mesh, err);
    }
    if (err) return 1;
  }

  // creates the smart laplacian optimization procedures
  SmartLaplacianSmoother* pass1 = new SmartLaplacianSmoother( obj_func );
  // SteepestDescent* pass1 = new SteepestDescent( obj_func );

  QualityAssessor stop_qa( mean_ratio );

  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 10 );
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // adds 1 pass of pass1 to mesh_set
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // launches optimization on mesh_set
  if (domain) {
    queue1.run_instructions(&mesh_and_domain, err);
  }
  else {
    queue1.run_instructions(mesh, err);
  }
  if (err) return 1;

  // Construct a MeshTSTT from the file
  int ierr = write_vtk_mesh( mesh, "smart-laplacian-result.vtk");
  if (ierr) return 1;
  // MeshWriter::write_vtk(mesh, "smart-laplacian-result.vtk", err);
  // if (err) return 1;
  cout << "Wrote \"smart-laplacian-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}


int run_quality_optimizer( MeshDomainAssoc& mesh_and_domain, MsqError& err )
{
  Mesh* mesh = mesh_and_domain.get_mesh();
  MeshDomain* domain = mesh_and_domain.get_domain();

  // creates an intruction queue
  InstructionQueue q;

  // do optimization
  const double eps = 0.01;
  IdealShapeTarget w;
  TShapeB1 mu;
  TQualityMetric metric( &w, &mu );
  PMeanPTemplate func( 1.0, &metric );

  SteepestDescent solver( &func );
  solver.use_element_on_vertex_patch();
  solver.do_jacobi_optimization();

  TerminationCriterion inner, outer;
  inner.add_absolute_vertex_movement( 0.5*eps );
  outer.add_absolute_vertex_movement( 0.5*eps );

  QualityAssessor qa( &metric );

  q.add_quality_assessor( &qa, err );
  if (err) return 1;
  q.set_master_quality_improver( &solver, err );
  if (err) return 1;
  q.add_quality_assessor( &qa, err );
  if (err) return 1;

// launches optimization on mesh_set
  if (domain) {
    q.run_instructions(&mesh_and_domain, err);
  }
  else {
    q.run_instructions(mesh, err);
  }
  if (err) return 1;

  // Construct a MeshTSTT from the file
  int ierr = write_vtk_mesh( mesh, "ideal-shape-result.vtk");
  if (ierr) return 1;
  // MeshWriter::write_vtk(mesh, "ideal-shape-result.vtk", err);
  // if (err) return 1;
  cout << "Wrote \"ideal-shape-result.vtk\"" << endl;

  print_timing_diagnostics( cout );
  return 0;
}


int run_solution_mesh_optimizer( MeshDomainAssoc& mesh_and_domain, MsqError& err )
{
  double OF_value = 0.01;

  Mesh* mesh = mesh_and_domain.get_mesh();
  MeshDomain* domain = mesh_and_domain.get_domain();

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  // IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  ConditionNumberQualityMetric* mean_ratio = new ConditionNumberQualityMetric();
  // VertexConditionNumberQualityMetric* mean_ratio = new VertexConditionNumberQualityMetric();
  // AspectRatioGammaQualityMetric* mean_ratio = new AspectRatioGammaQualityMetric();

  //ElementSolIndQM* solindqm = new ElementSolIndQM(solution_indicator);
  //MultiplyQualityMetric* mean_ratio = new MultiplyQualityMetric(new VertexConditionNumberQualityMetric(), solindqm, err);
  // ElementSolIndQM* mean_ratio = solindqm;

  // LocalSizeQualityMetric* mean_ratio = new LocalSizeQualityMetric();

  mean_ratio->set_averaging_method(QualityMetric::SUM_OF_RATIOS_SQUARED, err);
  if (err) return 1;

  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;

  // // creates the feas newt optimization procedures
  ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  // QuasiNewton* pass1 = new QuasiNewton( obj_func );
  // FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  pass1->use_global_patch();

  QualityAssessor stop_qa( mean_ratio );

  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  if (err) return 1;
  TerminationCriterion tc_outer;
  tc_inner.add_iteration_limit( 20 );
  tc_outer.add_iteration_limit( 5 );
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);
  pass1->set_debugging_level(3);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // launches optimization on mesh_set
  if (domain) {
    queue1.run_instructions(&mesh_and_domain, err);
  }
  else {
    queue1.run_instructions(mesh, err);
  }
  if (err) return 1;

  // Construct a MeshTSTT from the file
  int ierr = write_vtk_mesh( mesh, "solution-mesh-result.vtk");
  if (ierr) return 1;
  // MeshWriter::write_vtk(mesh, "solution-mesh-result.vtk", err);
  // if (err) return 1;
  cout << "Wrote \"solution-mesh-result.vtk\"" << endl;

  print_timing_diagnostics( cout );
  return 0;
}



int get_imesh_mesh( MBMesquite::Mesh** mesh, const char* file_name, int dimension )
{
  int err;
  iMesh_Instance instance = 0;
  iMesh_newMesh( NULL, &instance, &err, 0 ); CHECK_IMESH("Creation of mesh instance failed");

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet( instance, &root_set, &err );CHECK_IMESH("Could not get root set");

  iMesh_load( instance, root_set, file_name, 0, &err, strlen(file_name), 0 );CHECK_IMESH("Could not load mesh from file");

  iBase_TagHandle fixed_tag;
  iMesh_getTagHandle( instance, "fixed", &fixed_tag, &err, strlen("fixed") );
  // if (iBase_SUCCESS != err)
  {
    // get the skin vertices of those cells and mark them as fixed; we don't want to fix the vertices on a
    // part boundary, but since we exchanged a layer of ghost cells, those vertices aren't on the skin locally
    // ok to mark non-owned skin vertices too, I won't move those anyway
    // use MOAB's skinner class to find the skin

    // get all vertices and cells
    // make tag to specify fixed vertices, since it's input to the algorithm; use a default value of non-fixed
    // so we only need to set the fixed tag for skin vertices
    moab::Interface* mbi = reinterpret_cast<MBiMesh*>(instance)->mbImpl;
    moab::EntityHandle currset=0;
    moab::Tag fixed;
    int def_val = 0;
    err=0;
    moab::ErrorCode rval = mbi->tag_get_handle("fixed", 1, moab::MB_TYPE_INTEGER, fixed, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE, &def_val); MB_CHK_SET_ERR(rval, "Getting tag handle failed");
    moab::Range verts, cells, skin_verts;
    rval = mbi->get_entities_by_type(currset, moab::MBVERTEX, verts); MB_CHK_SET_ERR(rval, "Querying vertices failed");
    rval = mbi->get_entities_by_dimension(currset, dimension, cells); MB_CHK_SET_ERR(rval, "Querying elements failed");
    std::cout << "Found " << verts.size() << " vertices and " << cells.size() << " elements" << std::endl;

    moab::Skinner skinner(mbi);
    rval = skinner.find_skin(currset, cells, true, skin_verts); MB_CHK_SET_ERR(rval, "Finding the skin of the mesh failed"); // 'true' param indicates we want vertices back, not cells

    std::vector<int> fix_tag(skin_verts.size(), 1); // initialized to 1 to indicate fixed
    rval = mbi->tag_set_data(fixed, skin_verts, &fix_tag[0]); MB_CHK_SET_ERR(rval, "Setting tag data failed");
    std::cout << "Found " << skin_verts.size() << " vertices on the skin of the domain." << std::endl;

    // fix_tag.resize(verts.size(),0);
    // rval = mbi->tag_get_data(fixed, verts, &fix_tag[0]); MB_CHK_SET_ERR(rval, "Getting tag data failed");

    iMesh_getTagHandle( instance, "fixed", &fixed_tag, &err, strlen("fixed") );CHECK_IMESH("Getting tag handle (fixed) failed");

    // Set some arbitrary solution indicator
    moab::Tag solindTag;
    double def_val_dbl=0.0;
    rval = mbi->tag_get_handle("solution_indicator", 1, moab::MB_TYPE_DOUBLE, solindTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE, &def_val_dbl); MB_CHK_SET_ERR(rval, "Getting tag handle failed");
    solution_indicator.resize(cells.size(),0.01);
    for (unsigned i=0; i < cells.size()/4; i++)
      solution_indicator[i]=0.1;
    for (unsigned i=cells.size()/4; i < 2*cells.size()/4; i++)
      solution_indicator[i]=0.5;
    for (unsigned i=2*cells.size()/4; i < 3*cells.size()/4; i++)
      solution_indicator[i]=0.5;
    for (unsigned i=3*cells.size()/4; i < cells.size(); i++)
      solution_indicator[i]=0.5;

    rval = mbi->tag_set_data(solindTag, cells, &solution_indicator[0]); MB_CHK_SET_ERR(rval, "Setting tag data failed");

  }


  MsqError ierr;
  MBMesquite::MsqIMesh* imesh = new MBMesquite::MsqIMesh( instance, root_set, (dimension == 3 ? iBase_REGION : iBase_FACE), ierr, &fixed_tag );
  if (MSQ_CHKERR(ierr)) {
    delete imesh;
    cerr << err << endl;
    err = iBase_FAILURE;
    CHECK_IMESH("Creation of MsqIMesh instance failed");
    return 0;
  }

  *mesh = imesh;
  return iBase_SUCCESS;
}



int get_native_mesh( MBMesquite::Mesh** mesh, const char* file_name, int  )
{
  MsqError err;
  MBMesquite::MeshImpl* imesh = new MBMesquite::MeshImpl();
  imesh->read_vtk( file_name, err );
  if (err)
  {
    cerr << err << endl;
    exit(3);
  }
  *mesh = imesh;

  return iBase_SUCCESS;
}


int get_itaps_domain(MeshDomain** odomain, const char* filename)
{

  if (filename == 0 || strlen(filename) == 0) {
    *odomain=new PlanarDomain( PlanarDomain::XY );
    return 0;
  }

  int err;
  iGeom_Instance geom;
  iGeom_newGeom( "", &geom, &err, 0 ); CHECK_IGEOM("ERROR: iGeom creation failed");

#ifdef MOAB_HAVE_CGM_FACET
  FacetModifyEngine::set_modify_enabled(CUBIT_TRUE);
#endif

  iGeom_load( geom, filename, 0, &err, strlen(filename), 0 );
  CHECK_IGEOM( "Cannot load the geometry" );

  iBase_EntitySetHandle root_set;
  iGeom_getRootSet( geom, &root_set, &err );
  CHECK_IGEOM( "getRootSet failed!" );

    // print out the number of entities
  std::cout << "Model contents: " << std::endl;
  const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};
  int nents[4];
  for (int i = 0; i <= 3; ++i) {
    iGeom_getNumOfType( geom, root_set, i, &nents[i], &err );
    CHECK_IGEOM( "Error: problem getting entities after gLoad." );
    std::cout << gtype[i] << nents[i] << std::endl;
  }

  iBase_EntityHandle* hd_geom_ents;
  int csize=0, sizealloc=0;
  if (nents[3] > 0) {
    hd_geom_ents = (iBase_EntityHandle*)malloc(sizeof(iBase_EntityHandle)*nents[2]);
    csize = nents[2];
    iGeom_getEntities(geom, root_set, 2, &hd_geom_ents, &csize, &sizealloc, &err);
  }
  else {
    hd_geom_ents = (iBase_EntityHandle*)malloc(sizeof(iBase_EntityHandle)*nents[1]);
    csize = nents[1];
    iGeom_getEntities(geom, root_set, 1, &hd_geom_ents, &csize, &sizealloc, &err);
  }
  CHECK_IGEOM( "ERROR: Could not get entities" );

  *odomain = new MsqIGeom( geom, hd_geom_ents[0] );
  return iBase_SUCCESS;
}

