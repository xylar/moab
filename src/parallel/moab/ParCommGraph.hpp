/*
 * ParCommGraph.hpp
 *
 *  will be used to setup communication between 2 separate meshes, on
 *  different communicators;
 *  there are 3 communicators in play, one for each mesh, and one for the joined
 *  communicator, that spans both sets of processes
 *
 *  various methods should be available to partition meshes between the 2 communicators
 *  communicators are represented by their MPI groups, not by their communicators,, because
 *  the groups are always defined, irrespective of what tasks are they on.
 *
 *
 */
#include "moab_mpi.h"
#include "moab/Interface.hpp"
#include <map>


#ifndef SRC_PARALLEL_MOAB_PARCOMMGRAPH_HPP_
#define SRC_PARALLEL_MOAB_PARCOMMGRAPH_HPP_


namespace moab {

	class ParCommGraph {
	public:
	  ParCommGraph();
	  virtual ~ParCommGraph();

	  ParCommGraph(MPI_Comm joincomm, MPI_Group group1, MPI_Group group2);

	  void find_group_ranks(MPI_Group group, MPI_Comm join, std::vector<int> & ranks);
	  /*
	   * based on the number of elements on each task in group 1, partition for group 2, trivially
	   * establish how many elements are  sent from each task in group 1 to tasks in group 2
	   */
	  ErrorCode trivial_partition (std::vector<int> & numElemsPerTaskInGroup1);


	private:
	  MPI_Comm  comm;
	  MPI_Group gr1, gr2;
	  std::vector<int>  senderTasks;
	  std::vector<int>  receiverTasks;

	  // communication graph from group1 to group2;
	  //  graph[task1] = vec1; // vec1 is a stl vector of tasks in group2
	  std::map<int, std::vector<int> >  recv_graph; // to what tasks from group2 to send  (actual communication graph)
	  std::map<int, std::vector<int> >  recv_sizes; // how many elements to actually send from a sender task to receiver tasks
	  std::map<int, std::vector<int> >  sender_graph; // to what tasks from group2 to send  (actual communication graph)
	  std::map<int, std::vector<int> >  sender_sizes; // how many elements to actually send from a sender task to receiver tasks


};

} // namespace moab
#endif /* SRC_PARALLEL_MOAB_PARCOMMGRAPH_HPP_ */
