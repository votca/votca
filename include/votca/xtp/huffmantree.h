#ifndef HUFFMANTREE_H
#define HUFFMANTREE_H
#include <vector>
#include <votca/xtp/glink.h>
using namespace std;

namespace votca { namespace xtp {


struct huffmanNode{
    //huffmanNode * for the inner nodes, GLink * for the nodes on the last level before the leafs (The GLinks themselves represent the "leaf" level)
    huffmanNode * leftChild;
    huffmanNode * rightChild;
    GLink * rightLeaf;
    GLink * leftLeaf;
    double probability;
    bool isOnLastLevel=false;
};

class huffmanTree
{
private:
    void organizeProbabilities(huffmanNode *n, double add);
    void moveProbabilities(huffmanNode *n);
    vector <huffmanNode> htree;
    bool treeIsMade=false;
public:
    double escape_rate=0.0;
    vector <GLink> * events=nullptr;
    void makeTree();
    GLink *findHoppingDestination(double p);
    huffmanTree();
};

}}
#endif // HUFFMANTREE_H
