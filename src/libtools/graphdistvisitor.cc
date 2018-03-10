

bool GraphDistVisitor::queEmpty(){
  return edge_que_.empty(); 
}

// Add the distance to the node that has not yet been explored
void GraphDistVisitor::exploreNode_(pair<int,GraphNode> p_gn,Graph g, Edge ed){
  // Determine if the node has already been explored
  int vertex = p_gn.first;
  if(vertex==startingVertex_){
    GraphNode gn = p_gn.second;
    gn.int_vals["Dist"] = 0;
  }else{
    // Node has not been explored
    if(explored_.count(vertex)==0){ 
      int prev_vertex = ed.getOtherV(vertex);
      GraphNode gn_prev = g.getNode(prev_vertex); 
      gn.int_vals["Dist"]=gn_prev["Dist"]+1;
    }
  }
}

Edge GraphDistVisitor::getEdge_(Graph g){
  Edge ed = edge_que_.at(0).pop();
  if(edge_que.at(0).size()==0){
    egdge_que.pop_front();
  }
}

// Add edges to be explored
void GraphVisitor::addEdges_(Graph g, int vertex){
  auto eds = g.getNeighEdges(vertex);
  // Proceed to add them to queue if the vertices
  // they refer to are not already explored
 
  // If first edges to be added
  if(edge_que_.empty()){
    queue<Edge> first_que;
    for(auto ed : eds ){
      first_que.push_back(ed);
    }
    edge_que_.push_back(first_que);
  }else{

    if(edge_que.size()==1){
      queue<edge> next_que;
      for(auto ed : eds ){
        int neigh_vert = ed.getOther(vertex);
        if(explored_.count(neigh_vert)==0){
          next_que.push_back(ed);
        }
      }
      edge_que_.push_back(next_que);
    }else{
      for(auto ed : eds ){
        int neigh_vert = ed.getOther(vertex);
        if(explored_.count(neigh_vert)==0){
          // Add the edges to the next highest distance queue    
          edge_que_.at(1).push_back(ed);
        }
      }
    }
  }
}

