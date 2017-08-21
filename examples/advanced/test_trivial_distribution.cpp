#include <iostream>
#include "moab/Range.hpp"
#include <map>
#include <vector>
#include "assert.h"

using namespace moab;

void trivial_partition (int sender_rank, std::vector<int> & recvRanks,  std::vector<int> & number_elems_per_part,
    Range & owned, std::map<int, Range> & ranges_to_send)
{
  // first find out total number of elements to be sent from all senders
  int total_elems=0;
  for (size_t k=0; k<number_elems_per_part.size(); k++)
    total_elems+=number_elems_per_part[k];

  assert ( (int) (owned.size()) == number_elems_per_part[sender_rank]);

  int num_recv  =  ((int)recvRanks.size());
  // in trivial partition, every receiver should get about total_elems/num_receivers elements
  int num_per_receiver = (int)(total_elems/num_recv);
  int leftover = total_elems - num_per_receiver*num_recv;

  // so receiver k will receive  [starts[k], starts[k+1] ) interval
  std::vector<int> starts; starts.resize(num_recv+1);
  starts[0]=0;
  for (int k=0; k<num_recv; k++)
  {
    starts[k+1] = starts[k]+  num_per_receiver;
    if (k<leftover) starts[k+1]++;
  }

  std::cout << " intervals partitions: " ;
  for (int k=0; k< (int)starts.size(); k++)
  {
    std::cout<<" " << starts[k];
  }
  std::cout << "\n";

  int elems_before_sender_rank = 0;
  for (int k=0; k<sender_rank; k++)
    elems_before_sender_rank += number_elems_per_part[k];


  // receiver at index_recv will start receiving from owned range
  // find rank j that will receive the first subrange
  int j = 0;

  while (starts[j+1] <= elems_before_sender_rank )
    j++;

  /// so first element in range will go to receiver j interval:  [ starts[j], starts[j+1] )
  Range current = owned; // get the full range first, then we will subtract stuff, fot
  // the following ranges
  int current_index = elems_before_sender_rank;
  Range rleftover=current;
  while (elems_before_sender_rank + (int)owned.size() > starts[j+1] )
  {

    int upper_number = starts[j+1]-current_index;
    Range newr;
    newr.insert(current.begin(), current.begin() +upper_number);
    ranges_to_send[ recvRanks[j] ] = newr;

    rleftover = subtract(current, newr );
    current = rleftover;
    current_index = starts[j+1];
    j++;
  }
  if (!rleftover.empty())
    ranges_to_send[ recvRanks[j] ] = rleftover;

  return;
}

// will find each who is sending data to each receiver
// will be sent from rank 0 of senders to rank 0 of communicators, so the receivers are informed
// who is about to send elements to them; each receiver need to post receives from the senders for
// their data
void senders_trivial_partition(std::vector<int> & sndrRanks, std::vector<int> & number_elems_per_part,
    std::vector<int> & recvRanks, std::map<int, std::vector<int> > & recv_sets)
{
  // first find out total number of elements to be sent from all senders
  int total_elems=0;
  std::vector<int> accum;
  accum.push_back(0);

  int num_senders = (int) sndrRanks.size();

  for (size_t k=0; k<number_elems_per_part.size(); k++)
  {
    total_elems+=number_elems_per_part[k];
    accum.push_back(total_elems);
  }

  int num_recv  =  ((int)recvRanks.size());
  // in trivial partition, every receiver should get about total_elems/num_receivers elements
  int num_per_receiver = (int)(total_elems/num_recv);
  int leftover = total_elems - num_per_receiver*num_recv;

  // so receiver k will receive  [starts[k], starts[k+1] ) interval
  std::vector<int> starts; starts.resize(num_recv+1);
  starts[0]=0;
  for (int k=0; k<num_recv; k++)
  {
    starts[k+1] = starts[k]+  num_per_receiver;
    if (k<leftover) starts[k+1]++;
  }

  for (int j = 0; j < num_senders; j++ )
  {
    for (int k = 0; k<num_recv; k++ )
    {
      // if overlap:
      if (starts[k]<accum[j+1] && starts[k+1]> accum[j] )
      {
        recv_sets[ recvRanks[k]].push_back(sndrRanks[j]);
        if (starts[k] > accum[j+1])
          break ; // break the k loop
      }
    }
  }

}


int main(int argc, char * argv[])
{
  int sender_rank=1;
  std::vector<int> recvRanks;
  std::vector<int>  number_elems_per_part;
  recvRanks.push_back(3);
  recvRanks.push_back(4);
  recvRanks.push_back(5);
  recvRanks.push_back(6);
  recvRanks.push_back(7);
  Range owned(7,7+10-1);
  std::map<int, Range> ranges_to_send;
  number_elems_per_part.push_back(6);  number_elems_per_part.push_back(10); number_elems_per_part.push_back(6);

  std::cout<<" send sizes " ;
  for (int k=0; k< (int)number_elems_per_part.size(); k++)
  {
    std::cout<<" " << number_elems_per_part[k];
  }
  std::cout << "\n";
  std::cout<< " sender id: " << sender_rank << "  owned.size() " << owned.size() << " \n";
  std::cout << " receivers ranks : " << recvRanks.size() << " :" ;
  for (int k=0; k< (int)recvRanks.size(); k++)
  {
    std::cout<<" " << recvRanks[k];
  }
  std::cout << "\n";
  trivial_partition (sender_rank, recvRanks,  number_elems_per_part, owned, ranges_to_send);

  for (std::map<int, Range>::iterator it = ranges_to_send.begin(); it!=ranges_to_send.end(); it++ )
  {
    Range & ran = it->second;
    std::cout<< " receiver " << it->first << " receive range: [" << ran[0] << ", " << ran[ran.size()-1]  << "] \n";
  }

  // build communication matrix, each receiver will receive from what sender
  std::map<int, std::vector<int> > recv_sets; // map from receiver to its senders (in global space)

  std::vector<int> sndrRanks;
  sndrRanks.push_back(0); sndrRanks.push_back(1); sndrRanks.push_back(2);
  senders_trivial_partition(sndrRanks, number_elems_per_part, recvRanks, recv_sets);

  for (std::map<int, std::vector<int> >::iterator it=recv_sets.begin(); it!=recv_sets.end(); it++ )
  {
    int recv = it->first;
    std::vector<int> & senders = it->second;
    std::cout << " receiver: " << recv << " will receive from senders: " ;
    for (int k = 0; k<(int)senders.size(); k++ )
      std::cout << " " << senders[k];
    std::cout << "\n";
  }

  return 0;
}
