/*
 * ParCommGraph.cpp
 *
 */

#include "moab/ParCommGraph.hpp"

namespace moab {

ParCommGraph::ParCommGraph(MPI_Comm joincomm, MPI_Group group1, MPI_Group group2):
comm(joincomm), gr1(group1), gr2(group2) {
  // find out the tasks from each group, in the joint communicator
  find_group_ranks(gr1, comm, senderTasks);
  find_group_ranks(gr2, comm, receiverTasks);

  rootSender = rootReceiver = false;
  rankInGroup1 =  rankInGroup2 = rankInJoin = -1; // not initialized, or not part of the group
  int mpierr = MPI_Group_rank(gr1, &rankInGroup1);
  if (MPI_SUCCESS != mpierr)
    rankInGroup1 = -1;
  mpierr = MPI_Group_rank(gr2, &rankInGroup2);
  if (MPI_SUCCESS != mpierr)
    rankInGroup2 = -1;
  mpierr = MPI_Comm_rank(comm, &rankInJoin);
  if (MPI_SUCCESS != mpierr) // it should be a fatal error
    rankInJoin = -1;

  if (0==rankInGroup1)rootSender=true;
  if (0==rankInGroup2)rootReceiver=true;
}

ParCommGraph::~ParCommGraph() {
  // TODO Auto-generated destructor stub
}

// utility to find out the ranks of the processes of a group, with respect to a joint comm,
// which spans for sure the group
// it is used locally (in the constructor), but it can be used as a utility
void ParCommGraph::find_group_ranks(MPI_Group group, MPI_Comm joincomm, std::vector<int> & ranks)
{
   MPI_Group global_grp;
   MPI_Comm_group(joincomm, &global_grp);

   int grp_size;

   MPI_Group_size(group, &grp_size);
   std::vector<int> rks (grp_size);
   ranks.resize(grp_size);

   for (int i = 0; i < grp_size; i++)
     rks[i] = i;

   MPI_Group_translate_ranks(group, grp_size, &rks[0], global_grp, &ranks[0]);
   MPI_Group_free(&global_grp);
   return;
}

ErrorCode ParCommGraph::trivial_partition (std::vector<int> & numElemsPerTaskInGroup1)
{

  recv_graph.clear(); recv_sizes.clear();
  sender_graph.clear(); sender_sizes.clear();

  if (numElemsPerTaskInGroup1.size() != senderTasks.size())
    return MB_FAILURE; // each sender has a number of elements that it owns

  // first find out total number of elements to be sent from all senders
  int total_elems=0;
  std::vector<int> accum;
  accum.push_back(0);

  int num_senders = (int) senderTasks.size();

  for (size_t k=0; k<numElemsPerTaskInGroup1.size(); k++)
  {
    total_elems+=numElemsPerTaskInGroup1[k];
    accum.push_back(total_elems);
  }

  int num_recv  =  ((int)receiverTasks.size());
  // in trivial partition, every receiver should get about total_elems/num_receivers elements
  int num_per_receiver = (int)(total_elems/num_recv);
  int leftover = total_elems - num_per_receiver*num_recv;

  // so receiver k will receive  [starts[k], starts[k+1] ) interval
  std::vector<int> starts;
  starts.resize(num_recv+1);
  starts[0]=0;
  for (int k=0; k<num_recv; k++)
  {
    starts[k+1] = starts[k]+  num_per_receiver;
    if (k<leftover) starts[k+1]++;
  }

  // each sender will send to a number of receivers, based on how the
  // arrays starts[0:num_recv] and accum[0:sendr] overlap
  int lastUsedReceiverRank = 0; // first receiver was not treated yet
  for (int j = 0; j < num_senders; j++ )
  {
    // we could start the receiver loop with the latest receiver that received from previous sender
    for (int k = lastUsedReceiverRank; k<num_recv; k++ )
    {
      // if overlap:
      if (starts[k]<accum[j+1] && starts[k+1]> accum[j] )
      {
        recv_graph[receiverTasks[k]].push_back(senderTasks[j]);
        sender_graph[senderTasks[j]].push_back(receiverTasks[k]);

        // we still need to decide what is the overlap
        int sizeOverlap = 1; // at least 1, for sure
        //1
        if ( starts[k] >= accum[j]) // one end is starts[k]
        {
          if (starts[k+1] >= accum[j+1]) // the other end is accum[j+1]
            sizeOverlap = accum[j+1]-starts[k];
          else //
            sizeOverlap = starts[k+1]-starts[k];
        }
        else // one end is accum[j]
        {
          if (starts[k+1] >= accum[j+1]) // the other end is accum[j+1]
            sizeOverlap = accum[j+1]-accum[j];
          else
            sizeOverlap = starts[k+1]-accum[j];
        }
        recv_sizes[receiverTasks[k]].push_back(sizeOverlap); // basically, task k will receive from
                                                        //   sender j, sizeOverlap elems
        sender_sizes[senderTasks[j]].push_back(sizeOverlap);
        if (starts[k] > accum[j+1])
        {
          lastUsedReceiverRank = k-1; // so next k loop will start a little higher, we probably
                                            // finished with first few receivers (up to receiver lastUsedReceiverRank)
          break ; // break the k loop, we distributed all elements from sender j to some receivers
        }
      }
    }
  }

  return MB_SUCCESS;
}
ErrorCode ParCommGraph::pack_receivers_graph(std::vector<int> & packed_recv_array )
{
  // it will basically look at local data, to pack communication graph, each receiver task will have to post receives
  // for each sender task that will send data to it;
  // the array will be communicated to root receiver, and eventually distributed to receiver tasks

  /*
   * packed_array will have receiver, number of senders, then senders, etc
   */
  for (std::map<int, std::vector<int> >::iterator it=recv_graph.begin(); it!=recv_graph.end(); it++ )
  {
    int recv = it->first;
    std::vector<int> & senders = it->second;
    packed_recv_array.push_back(recv);
    packed_recv_array.push_back( (int) senders.size() );

    for (int k = 0; k<(int)senders.size(); k++ )
      packed_recv_array.push_back( senders[k] );
  }

  return MB_SUCCESS;
}

ErrorCode ParCommGraph::distribute_sender_graph_info(MPI_Comm senderComm, std::vector<int> &sendingInfo )
{
  // only the root of the sender has all info
  // it will use direct send/receives, for number of elems to be sent to each receiver, in an array
  if (rootSender)
  {
    // just append the local info for the array , and send some data for each receiver
    int nrecv0 = (int) recv_graph[senderTasks[0]].size();
    for (int k=0; k<nrecv0; k++)
    {
      sendingInfo.push_back(recv_graph[senderTasks[0]][k]);
      sendingInfo.push_back(recv_sizes[senderTasks[0]][k]);
    }
    // the rest of the info will be sent for each sender task
    for (int j=1; j< (int) senderTasks.size(); j++)
    {
      std::vector<int> array_to_send;
      // in the sender comm and sender group , rank to send to is j
      for (int k=0; k<(int) recv_graph[senderTasks[j]].size(); k++)
      {
        array_to_send.push_back( recv_graph[senderTasks[j]][k]);
        array_to_send.push_back( recv_sizes[senderTasks[j]][k]);
      }

      int ierr = MPI_Send(&array_to_send[0], (int)array_to_send.size(), MPI_INT, j, 11, senderComm);
      if (MPI_SUCCESS != ierr)
        return MB_FAILURE;
    }
  }
  else
  {
    // expect to receive some info about where to send local data
    // local array it is max the size of 2 times nreceivers;
    int sizeBuffMax = 2 * (int) receiverTasks.size();
    std::vector<int> data(sizeBuffMax);
    MPI_Status status;
    int ierr = MPI_Recv(&data[0], sizeBuffMax, MPI_INT, 0, 11, senderComm, &status);
    if (MPI_SUCCESS != ierr)
      return MB_FAILURE;
    int count;
    // find out how much data we actually got
    MPI_Get_count(&status, MPI_INT, &count);
    for (int k=0; k<count; k++)
      sendingInfo.push_back(data[k]);

  }
  return MB_SUCCESS;
}

} // namespace moab
