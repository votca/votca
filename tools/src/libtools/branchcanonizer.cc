/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writingraph, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

namespace votca {
namespace tools {
/*
void BranchCanonizer::canonize() {

  std::vector<Index> starting_vertices  = findUniqueStartingVertex_();

  // For the first iteration we will assign distances to starting nodes, but
  // that is only done on the first itteration
  std::map<Index, Graph> chosen_graphs = assignDistances_(starting_vertices);

  // For memory reasons cycle graphs one at a time
  for ( std::pair<Index,Graph> vert_graph_temp : chosen_graphs){

    Graph g = vert_graph_temp.second;
    // Continue the cycle until there are no more parallel branches
    do {
      Index starting_vertex = vert_graph_temp.first;
      ReducedGraph reduced_graph = reduceGraph_(g);
      assignLevelsToBranches_(reduced_graphs,starting_vertex);
      std::map<Index,std::vector<Branch>> branch_levels =
getBranchLevels(reduced_graph) BranchGroup branch_groups =
separateBranches(branch_levels); BranchTree tree =
calculateBranchTrees_(branch_groups); g = newGraph_(tree, branch_groups); }
while (parallel_branches.size()>0);

    // Now compare the string id of the tree with the existing string ids the
largest
    // string id will be chosen
  }
}

// Return the graph with the starting vertex
static std::map<Index, Graph> assignDistances_(Index
starting_vertex,std::vector<Graph> & chosen_graphs) {

  std::map<Index, Graph> chosen_graphs;
  for( const Index & starting_vertex : starting_vertices ) {
    GraphDistanceVisitor graph_visitor;
    graph_visitor.setStartingVertex(starting_vertex);
    Graph graph_temp = graph;
    exploreGraph(graph_temp, graph_visitor);
    std::string temp_struct_id = graph_temp.getId();
    if (chosenId.compare(temp_struct_id) < 0) {
      // If the id is the smallest encountered thus far the vector of graphs is
      // cleared
      chosen_graphs.clear();
      chosenId = temp_struct_id;
      chosen_graphs[starting_vertex] = graph_temp;
    } else if ( chosenId.compare(temp_struct_id) == 0){
      // If the string ids are the same then all the graphs with the same string
      // ids must be further analized
      chosen_graphs[starting_vertex] = graph_temp;
    }
  }
  return chosen_graphs;
}

static Graph newGraph_(BranchTree tree,BranchGroup branch_groups){
  // So we need to create a new graph
}

static void assignLevelsToBranches_( Graph * graph, Index starting_vertex ) {
  GraphBranchLevelVisitor gv;
  gv.setStartingVertex(starting_vertex);
  exploreGraph( graph, gv);
}

// Return all the edges sorted into groups based on the "Level" property
static std::map<Index,std::vector<Branch>> getBranchLevels(reduced_graph) {
  std::map<Index,std::vector<Branch>> branch_levels;
  std::vector<Edge> edges;
  for ( Edge & ed : edges ) {
    branch_levels[ed.getInt("Level")].push_back(ed);
  }
  return branch_levels;
}

static BranchTree calculateBranchTrees_(BranchGroup branch_groups){
  // Note that merges do not share branches of the level above
  //           splits will share branches of the level above
  BranchTree tree(branch_groups.max_level);

  for ( int level = branch_groups.max_level; level>=0; level--){
    for( Branch & branch : perpendicular_branches[level] ){
      tree.addBranch(level,branch);
    }
    tree.stepLevel();
  }
  return tree;
}

struct BranchGroups {
  int max_level = 0;
  std::map<Index,std::vector<Branch>> parallel_branches;
  std::map<Index,std::vector<Branch>> perpendicular_branches;
};


// Of the branches that are in the same level get the ones that are connected
static BranchGroups separateBranches(std::map<Index,std::vector<Branch>>
branch_levels){

  BranchGroups branch_group;

  for ( std::pair<Index, vector<Branch> > & level : branch_levels ){
    if(level.first > branch_group.max_level ) {
      branch_group.max_level = level.first;
    }
    std::map<size_t,int> tagged_branch_indices;
    for ( size_t ind = 0; ind < level.second.size(); ++ind ){
      Branch branch1 = level.second.at(ind);
      for ( size_t ind2 = ind+1; ind2 < level.second.size(); ++ind2){
        Branch branch2 = level.second.at(ind2);
        if( branch1.getSource() == branch2.getSource() ){
          tagged_branch_indices[ind]++;
          tagged_branch_indices[ind2]++;
        }
        if (branch1.getSource() == branch2.getTerminal() ){
          tagged_branch_indices[ind]++;
          tagged_branch_indices[ind2]++;
        }
        if (branch1.getTerminal() == branch2.getSource() ) {
          tagged_branch_indices[ind]++;
          tagged_branch_indices[ind2]++;
        }
        if (branch1.getTerminal() == branch2.getTerminal() ) {
          tagged_branch_indices[ind]++;
          tagged_branch_indices[ind2]++;
        }
      }
    } // for ind

    // Any branch tagged twice or more is a parallel branch
    for ( pair<size_t, int> & ind_and_count : tagged_branch_indices ){
      size_t ind = ind_and_count.first;
      if( ind_and_count.second > 1 ){
        branch_group.parallel_branches[level.first].push_back(level.second.at(ind));
      } else {
        branch_group.perpendicular_branches[level.first].push_back(level.second.at(ind));

      }
    }
  } // for level
  return branch_group;
}

static std::vector<Index> findUniqueStartingVertex_(Graph & graph) {

  std::vector<Index> starting_vertices =
getPotentialStartingVerticesBasedOnDegree_(graph); starting_vertices =
filterPotentialVerticesByNodeContent_(starting_vertices); return
starting_vertices;
}

static std::vector<Index> getPotentialStartingVerticesBasedOnDegree_(Graph &
graph) { Index maxD = graph.getMaxDegree();
  // Get the vertices with this degree
  std::vector<Index> vertices = graph.getVerticesDegree(maxD);
  return vertices;
}

static std::vector<Index> filterPotentialVerticesByNodeContent_(Graph & graph,
std::vector<Index> potential_starting_vertices) {
  // Get the nodes and determine which node has the greatest stringID
  // When compared using compare function
  std::string str_id = "";
  std::vector<Index> graph_node_ids;
  for (const Index& vertex : potential_starting_vertices) {
    GraphNode graph_node = graph.getNode(vertex);
    Index comp_int = str_id.compare(graph_node.getStringId());
    if (comp_int > 0) {
      str_id = graph_node.getStringId();
      graph_node_ids.clear();
      graph_node_ids.push_back(vertex);
    } else if (comp_int == 0) {
      graph_node_ids.push_back(vertex);
    }
  }

  if (graph_node_ids.size() == 0) {
    graph_node_ids = vertices;
  }

  return graph_node_ids;

}
*/
}
}  // namespace votca
