#include <vector>
#include <list>
#include <stdlib.h>
#include <queue>
#include <votca/xtp/huffmantree.h>
#include <votca/xtp/glink.h>
using namespace std;



namespace votca {namespace xtp {

huffmanTree::huffmanTree()
{
}

//https://en.wikipedia.org/wiki/Huffman_coding
void huffmanTree::makeTree(){
    if (!events) throw runtime_error("Error in Huffmantree::makeTree : Pointer to Events not set!");

    //queue of the nodes, sorted by probability
    auto compare = [](huffmanNode * n1, huffmanNode * n2)
    { return n1->probability>n2->probability;};
    priority_queue<huffmanNode *,vector<huffmanNode *>, decltype(compare)> queue(compare);

    htree=vector<huffmanNode>(events->size()%2?events->size():events->size()-1);

    auto comp2 = [](GLink * e1, GLink * e2){
        return e1->rate>e2->rate;
    };
    priority_queue<GLink *, vector<GLink *>, decltype(comp2)> eventQueue(comp2);

    int i=0;
    for (GLink &e:*events){
        eventQueue.push(&e);
     }
    while (eventQueue.size()>1){
        htree[i].isOnLastLevel=true;
        htree[i].leftLeaf=eventQueue.top();
        eventQueue.pop();
        htree[i].rightLeaf=eventQueue.top();
        eventQueue.pop();
        htree[i].probability=(htree[i].leftLeaf->rate+htree[i].rightLeaf->rate)/escape_rate;
        queue.push(&(htree[i]));
        i++;
    }
    if (!eventQueue.empty()){
        htree[i].isOnLastLevel=true;
        htree[i].rightLeaf=eventQueue.top();
        htree[i].leftLeaf=eventQueue.top();
        htree[i].probability=2*htree[i].leftLeaf->rate/escape_rate;
        queue.push(&(htree[i]));
        i++;
    }

    //now connect the hnodes, making a new one for every connection:
    //always take the two nodes with the smallest probability and "combine" them, repeat, until just one node (the root) is left.
    huffmanNode * h1;
    huffmanNode * h2;
    while (queue.size()>1){
        h1=queue.top();
        queue.pop();
        h2=queue.top();
        queue.pop();
        htree[i].probability=h1->probability+h2->probability;
        htree[i].leftChild=h1;
        htree[i].rightChild=h2;
        queue.push(&(htree[i]));
        i++;
    }
    //reorganize the probabilities: in every node, add the probability of one subtree ("small")
    //to all nodes of the other subtree.
    organizeProbabilities(&htree[htree.size()-1],0);
    moveProbabilities(&htree[htree.size()-1]);
    treeIsMade=true;
}



void huffmanTree::organizeProbabilities(huffmanNode * n,double add){

    //adds "add" to the probability, then calls itself recursively.
    //this calculates the probabilities needed to traverse the tree quickly
    n->probability+=add;
    //if leftId=-1 (=> node is leaf), returns
    if (n->isOnLastLevel){
        return;
    }

    organizeProbabilities(n->leftChild,add+n->rightChild->probability);
    organizeProbabilities(n->rightChild,add);
}



void huffmanTree::moveProbabilities(huffmanNode * n){
    if (n->isOnLastLevel){
        n->probability-=n->leftLeaf->rate/escape_rate;
    }
    else{
        n->probability=n->rightChild->probability;
        moveProbabilities(n->rightChild);
        moveProbabilities(n->leftChild);
    }
}



GLink * huffmanTree::findHoppingDestination(double p){
    if (!treeIsMade) throw runtime_error("Tried to find Hopping Destination without initializing the Huffmantree first!");
    huffmanNode * node=&(htree[htree.size()-1]);
    while (!node->isOnLastLevel){
        if (p>node->probability) node=node->leftChild;
        else node=node->rightChild;
    }
    return (p>node->probability?node->leftLeaf:node->rightLeaf);
}

}}
