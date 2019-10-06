#ifndef APART_NODE_WBUBBLE_HPP
#define APART_NODE_WBUBBLE_HPP
// ==================================================================
// APart_Node_wBubble.hpp
//
// Adopt the partitioned node indices and enrich them with cell
// bubbles.
// 
// The bubble index is appended after the original nodal indices.
// In element ee, the bubble node indices are
//
//   nFunc + ee x num_bubble_per_cell + 0
//   nFunc + ee x num_bubble_per_cell + 1
//   ...
//   nFunc + ee x num_bubble_per_cell + num_bubble_per_cell - 1
//
// The number of bubble nodes in this partition is
// num_bubble_per_cell x num_local_cell
// 
// In the local_to_global array, the bubble nodes are appended after
// the ghost nodes, and hence has the structure
//       {node_local} + {node_ghost} + {node_bubble}.
// 
// The ntotalnode = local_to_global.size()
//                = nlocalnode + nghostnode + num_total_bubbles
// 
// The nlocalnode = nlocalnode + num_total_bubbles
//
// Author: Ju Liu
// Date: Nov. 27 2017
// ==================================================================
#include "APart_Node.hpp"
#include "ALocal_Elem.hpp"

class APart_Node_wBubble : public APart_Node
{
  public:
    APart_Node_wBubble(
        const std::string &fileBaseName, const int &rank,
        const int &nbub_per_cell );

    virtual ~APart_Node_wBubble();

    virtual void print_info() const;

    virtual int get_nbubblenode() const {return num_total_bubbles;}

    // 0 <= index < num_total_bubbles
    virtual int get_bubble_node(const int &index) const
    {return node_bubble[index];}

  private:
    int num_bubble_per_cell, num_local_cell, num_total_bubbles;
    
    std::vector<int> node_bubble;
};

#endif
