/*
 * ParCommGraph.hpp
 *
 *  will be used to setup communication between 2 separate meshes, on
 *  different communicators;
 *  there are 3 communicators in play, one for each mesh, and one for the joined
 *  communicator, that spans both sets of processes
 *
 *  various methods should be available to partition meshes between the 2 communicators
 *  communicators are represented by their MPI groups, not by their communicators, because
 *  the groups are always defined, irrespective of what tasks are they on.
 *  Some of the methods in here are executed over the sender communicator, some are over the receiver communicator
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
	  // ParCommGraph();
	  virtual ~ParCommGraph();

	  // collective constructor, will be called on all sender tasks and receiver tasks
	  ParCommGraph(MPI_Comm joincomm, MPI_Group group1, MPI_Group group2, int coid1, int coid2);

	  /**
	    \brief find ranks of a group with respect to an encompassing communicator

	    <B>Operations:</B> Local, usually called on root process of the group

      \param[in]  joincomm (MPI_Comm)
	    \param[in]  group (MPI_Group)
	    \param[out] ranks ( std::vector<int>)  ranks with respect to the joint communicator
	  */
	  void find_group_ranks(MPI_Group group, MPI_Comm join, std::vector<int> & ranks);

	  /**
	    \brief  Based on the number of elements on each task in group 1, partition for group 2, trivially

     <B>Operations:</B> Local, usually called on root process of the group

      Note:  establish how many elements are sent from each task in group 1 to tasks in group 2
            This call is usually made on a root / master process, and will construct local maps that are member data,
            which contain the communication graph, in both directions
            Also, number of elements migrated/exchanged between each sender/receiver

       \param[in]  numElemsPerTaskInGroup1 (std::vector<int> &)  number of elements on each sender task
	   */
	  ErrorCode compute_trivial_partition (std::vector<int> & numElemsPerTaskInGroup1);

	  /**
	     \brief  pack information about receivers view of the graph, for future sending to receiver root

        <B>Operations:</B> Local, usually called on root process of the group

	     \param[out] packed_recv_array
	       packed data will be sent eventually to the root of receivers, and distributed from there, and
	         will have this information, for each receiver, concatenated
	          receiver 1 task, number of senders for receiver 1, then sender tasks for receiver 1, receiver 2 task,
	            number of senders for receiver 2, sender tasks for receiver 2, etc
	   */
	  ErrorCode pack_receivers_graph(std::vector<int> & packed_recv_array );

	  /**
	   * \brief  distribute send information (as decided on root) to all senders in the group
	   * <B>Operations:</B> collective, needs to be called on all tasks on sender group
	   *
	   * sending info will contain the rank of the receivers and number of elements
	   *   (sizes) for each receiver  recv1, size1, recv2, size2, ...
	   *   size of the array will be the number of receivers times 2
	   * \param[in]      senderComm(MPI_Comm)
	   * \param[out]     sendingInfo(std::vector<int> & )
	   */
	  ErrorCode distribute_sender_graph_info(MPI_Comm senderComm, std::vector<int> &sendingInfo );

	  // get methods for private data
	  bool is_root_sender() { return rootSender;}

	  bool is_root_receiver () { return rootReceiver;}

	  int get_index_sender() { return index_sender;}
	  int get_index_receiver() { return index_receiver;}

	  int sender(int index) {return senderTasks[index];}

	  int receiver(int index) {return receiverTasks[index];}

	  // setter methods for private data
	  void set_index_sender(int ix1) {index_sender=ix1;}
	  void set_index_receiver(int ix2) {index_receiver=ix2;}

	  int get_component_id1(){return compid1;}
	  int get_component_id2(){return compid2;}

	  // return local graph for a specific root
	  ErrorCode split_owned_range (int sender_rank, Range & owned, std::map<int, Range> & split_ranges);

	private:
	  MPI_Comm  comm;
	  MPI_Group gr1, gr2;
	  std::vector<int>  senderTasks;
	  std::vector<int>  receiverTasks;
	  bool rootSender;
	  bool rootReceiver;
	  int rankInGroup1, rankInGroup2;
	  int rankInJoin, joinSize;
	  int index_sender, index_receiver; // indices in the list of local graphs referred by this application (not used yet)
	  int compid1, compid2;

	  // communication graph from group1 to group2;
	  //  graph[task1] = vec1; // vec1 is a stl vector of tasks in group2
	  std::map<int, std::vector<int> >  recv_graph; // to what tasks from group2 to send  (actual communication graph)
	  std::map<int, std::vector<int> >  recv_sizes; // how many elements to actually send from a sender task to receiver tasks
	  std::map<int, std::vector<int> >  sender_graph; // to what tasks from group2 to send  (actual communication graph)
	  std::map<int, std::vector<int> >  sender_sizes; // how many elements to actually send from a sender task to receiver tasks


};

} // namespace moab
#endif /* SRC_PARALLEL_MOAB_PARCOMMGRAPH_HPP_ */
