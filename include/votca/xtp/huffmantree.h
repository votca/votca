#ifndef HUFFMANTREE_H
#define HUFFMANTREE_H
#include <vector>
#include <votca/xtp/glink.h>
using namespace std;

namespace votca { namespace xtp {




class huffmanTree
{

    struct huffmanNode{
    //huffmanNode * for the inner nodes, GLink * for the nodes on the last level before the leafs (The GLinks themselves represent the "leaf" level)
    huffmanNode * leftChild;
    huffmanNode * rightChild;
    GLink * rightLeaf;
    GLink * leftLeaf;
    double probability;
    bool isOnLastLevel=false;
};

public:
    void makeTree();
    GLink *findHoppingDestination(double p);
    void setEvents(std::vector<GLink> * v);

private:
    void addProbabilityFromRightSubtreeToLeftSubtree(huffmanNode *n, double add);
    void moveProbabilitiesFromRightSubtreesOneLevelUp(huffmanNode *n);
    vector <huffmanNode> htree;
    bool treeIsMade=false;
    double escape_rate;
    vector <GLink> * events=nullptr;

};

}}
#endif // HUFFMANTREE_H
