// 
// Usage:
// tools/mbslavepart -d 2 -m mpas/x1.2562.grid.h5m -s mpas/x1.10242.grid.h5m -o mpas_slave.h5m -e 1e-8 -b 1e-6 -O
// 
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/BVHTree.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

using namespace moab;

// A function to get the non-default value from a std::map
template <typename K, typename V>
static V get_map_value(const  std::map <K,V> & m, const K & key, const V & defval ) {
   typename std::map<K,V>::const_iterator it = m.find( key );
   if ( it == m.end() ) {
      return defval;
   }
   else {
      return it->second;
   }
}

int main(int argc, char* argv[])
{
  int proc_id = 0, size = 1, dimension=3;
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

  double tolerance=1e-6, btolerance=1;
  std::string masterfile, slavefile, outfile("slavemesh.vtk");
  bool keepsparts=false;
  ProgOptions opts;

  opts.addOpt<std::string>("master,m", "Master mesh filename", &masterfile);
  opts.addOpt<std::string>("slave,s", "Slave mesh filename", &slavefile);
  opts.addOpt<std::string>("output,o", "Output partitioned mesh filename", &outfile);
  opts.addOpt<int>("dim,d", "Dimension of entities to use for partitioning", &dimension);
  opts.addOpt<double>("eps,e", "Tolerance for the point search", &tolerance);
  opts.addOpt<double>("beps,b", "Tolerance for the bounding box search", &btolerance);
  opts.addOpt<void>("keep,K", "Keep the existing partitions in the slave mesh (use PARALLEL_PARTITION_SLAVE instead)", &keepsparts);
  opts.parseCommandLine(argc, argv);

  if (masterfile.empty() || slavefile.empty())
  {
    opts.printHelp();
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  ErrorCode error;
  Core *mbCore = new Core();

  //Set the read options for parallel file loading
  const std::string partition_set_name = "PARALLEL_PARTITION";
  const std::string global_id_name = "GLOBAL_ID";
  std::vector<std::string> read_opts, write_opts;
  std::string read_options, write_options;

  if (size > 1) {
    read_options = ";PARALLEL=READ_PART;PARTITION="+partition_set_name+";PARALLEL_RESOLVE_SHARED_ENTS";
    write_options = ";;PARALLEL=WRITE_PART";
  }

  EntityHandle masterfileset, slavefileset;
  error = mbCore->create_meshset(moab::MESHSET_SET, masterfileset);MB_CHK_ERR(error);
  error = mbCore->create_meshset(moab::MESHSET_SET, slavefileset);MB_CHK_ERR(error);

  //Load file
  error = mbCore->load_file(masterfile.c_str(), &masterfileset, read_options.c_str());MB_CHK_ERR(error);
  error = mbCore->load_file(slavefile.c_str(), &slavefileset, read_options.c_str());MB_CHK_ERR(error);

  Tag gidtag=0, parttag=0, sparttag=0;
  error = mbCore->tag_get_handle(partition_set_name.c_str(), parttag);MB_CHK_ERR(error);
  error = mbCore->tag_get_handle(global_id_name.c_str(), gidtag);MB_CHK_ERR(error);
  if (keepsparts) {
    int defpart=0;
    error = mbCore->tag_get_handle(std::string(partition_set_name+"_SLAVE").c_str(), 1, MB_TYPE_INTEGER, sparttag, MB_TAG_CREAT | MB_TAG_SPARSE, &defpart);MB_CHK_ERR(error);
  }

  Range melems, msets, selems, ssets;
  // TODO: expand and add other dimensional elements
  error = mbCore->get_entities_by_dimension(slavefileset, dimension, selems);MB_CHK_ERR(error);

  // Get the partition sets on the master mesh
  std::map<int, int> mpartvals;
  error = mbCore->get_entities_by_type_and_tag(masterfileset, MBENTITYSET, &parttag, NULL, 1, msets, moab::Interface::UNION, true);MB_CHK_ERR(error);
  if (msets.size() == 0) {
    std::cout << "No partition sets found in the master mesh. Quitting..." << std::endl;
    exit (1);
  }

  for (unsigned i=0; i < msets.size(); ++i) {
    EntityHandle mset=msets[i];

    moab::Range msetelems;
    error = mbCore->get_entities_by_dimension(mset, dimension, msetelems);MB_CHK_ERR(error);
    melems.merge(msetelems);
    
    int partID;
    error = mbCore->tag_get_data(parttag, &mset, 1, &partID);MB_CHK_ERR(error);

    // Get the global ID and use that as the indicator
    std::vector<int> gidMelems(msetelems.size());
    error = mbCore->tag_get_data(gidtag, msetelems, gidMelems.data());MB_CHK_ERR(error);

    for (unsigned j=0; j < msetelems.size(); ++j)
      mpartvals[gidMelems[j]]=partID;
      // mpartvals[msetelems[j]]=partID;
  }

  std::cout << "Found " << melems.size() << " elements in master mesh with " << msets.size() << " partition sets." << std::endl;
  msets.clear();

  std::map<int, moab::Range > spartvals;
  int npoints_notfound=0;
  {
    EntityHandle tree_root;
    moab::AdaptiveKDTree tree(mbCore);
    // moab::BVHTree tree(mbCore);
    error = tree.build_tree(melems, &tree_root);MB_CHK_ERR(error);

    for (size_t ie=0; ie < selems.size(); ie++) {
      moab::EntityHandle selem,leaf;
      double point[3];
      selem = selems[ie];

      // Get the element centroid to be queried
      error = mbCore->get_coords(&selem, 1, point);MB_CHK_ERR(error);

      // Search for the closest source element in the master mesh corresponding
      // to the target element centroid in the slave mesh 
      error = tree.point_search( point, leaf, tolerance, btolerance );MB_CHK_ERR(error);

      if (leaf == 0) {
        leaf = masterfileset; // FIXME: This is terrible -- linear search.
      }

      std::vector<moab::EntityHandle> leaf_elems;
      // We only care about the dimension that the user specified.
      // MOAB partitions are ordered by elements anyway.
      error = mbCore->get_entities_by_dimension( leaf, dimension, leaf_elems);MB_CHK_ERR(error);

      // Now get the master element centroids so that we can compute 
      // the minimum distance to the target point
      std::vector<double> centroids(leaf_elems.size()*3);
      error = mbCore->get_coords(&leaf_elems[0], leaf_elems.size(), &centroids[0]);MB_CHK_ERR(error);

      if (!leaf_elems.size())
        std::cout << ie << ": " << " No leaf elements found." << std::endl;

      double dist=1e10;
      int pinelem=-1;
      for (size_t il=0; il < leaf_elems.size(); ++il) {
        const double *centroid = &centroids[il*3];
        const double locdist = std::pow(point[0]-centroid[0],2)+std::pow(point[1]-centroid[1],2)+std::pow(point[2]-centroid[2],2);

        if (locdist < dist) {
          dist = locdist;
          pinelem = il;
        }
      }

      if (pinelem < 0) {
        std::cout << ie << ": [Error] - Could not find a minimum distance within the leaf nodes." << std::endl;
        npoints_notfound++;
      }
      else {
        int gidMelem;
        error = mbCore->tag_get_data(gidtag, &leaf_elems[pinelem], 1, &gidMelem);MB_CHK_ERR(error);

        // if (mpartvals[ gidMelem ])
        int mpartid = get_map_value(mpartvals, gidMelem, -1);
        if (mpartid < 0)
          std::cout << "[WARNING]: Part number for element " << leaf_elems[pinelem] << " with global ID = " << gidMelem << " not found.\n";

        spartvals[ mpartid ].insert(selems[ie]);
      }
    }
    error = tree.reset_tree();MB_CHK_ERR(error);
  }
  if (npoints_notfound) std::cout << "Could not find " << npoints_notfound << " points in the master mesh" << std::endl;

  // Find parallel partition sets in the slave mesh - and delete it since we are going to overwrite the sets
  if (!keepsparts) {
    error = mbCore->get_entities_by_type_and_tag(slavefileset, MBENTITYSET, &parttag, NULL, 1, ssets, moab::Interface::UNION);MB_CHK_ERR(error);
    std::cout << "Deleting " << ssets.size() << " sets in the slave mesh" << std::endl;
    error = mbCore->delete_entities(ssets);MB_CHK_ERR(error);
    ssets.clear();
  }

  for (std::map<int, moab::Range >::iterator it = spartvals.begin(); it != spartvals.end(); ++it) {
    int partID = it->first;
    moab::EntityHandle pset;
    error = mbCore->create_meshset(moab::MESHSET_SET, pset);MB_CHK_ERR(error);
    error = mbCore->add_entities(pset, it->second);MB_CHK_ERR(error);
    error = mbCore->add_parent_child(slavefileset, pset);MB_CHK_ERR(error);

    if (keepsparts) {
      error = mbCore->tag_set_data(sparttag, &pset, 1, &partID);MB_CHK_ERR(error);
    }
    else {
      error = mbCore->tag_set_data(parttag, &pset, 1, &partID);MB_CHK_ERR(error);
    }
  }

  error = mbCore->write_file("mpas_slave.vtk", NULL, write_options.c_str(), &slavefileset, 1);MB_CHK_ERR(error);
  error = mbCore->write_file(outfile.c_str(), NULL, write_options.c_str(), &slavefileset, 1);MB_CHK_ERR(error);
  
  delete mbCore;

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  exit(0);
}
