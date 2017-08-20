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

  int elems_before_sender_rank = 0;
  for (int k=0; k<sender_rank; k++)
    elems_before_sender_rank += number_elems_per_part[k];

  // receiver at index_recv will start receiving from owned range
  // find rank j that will receive the first subrange
  int j = 0;

  while (starts[j+1] < elems_before_sender_rank )
    j++;

  /// so first element in range will go to receiver j interval:  [ starts[j], starts[j+1] )
  ranges_to_send[ recvRanks[j] ] = owned; // get the full range first, then we will subtract stuff, fot
  // the following ranges
  while (elems_before_sender_rank + owned.size() > starts[j+1] )
  {
    Range & current = ranges_to_send[ recvRanks[j] ] ; // previous range that needs to be chopped off

  }

  return;
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
  Range owned(7,7+7-1);
  std::map<int, Range> ranges_to_send;
  number_elems_per_part.push_back(6);  number_elems_per_part.push_back(7); number_elems_per_part.push_back(6);
  trivial_partition (sender_rank, recvRanks,  number_elems_per_part, owned, ranges_to_send);

  return 0;
}
