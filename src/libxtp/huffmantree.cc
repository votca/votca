#include <vector>
#include <list>
#include <stdlib.h>
#include <queue>
#include <votca/xtp/huffmantree.h>
#include <votca/xtp/glink.h>
using namespace std;



namespace votca {namespace xtp {


void huffmanTree::setEvents(vector<GLink> * v){
    this->events=v;
}

//https://en.wikipedia.org/wiki/Huffman_coding
void huffmanTree::makeTree(){
    if (!events) throw runtime_error("Error in Huffmantree::makeTree : Pointer to Events not set!");

    //queue of the nodes, sorted by probability
    auto compare = [](huffmanNode * n1, huffmanNode * n2)
    { return n1->probability>n2->probability;};

    //priority queues, because the algorithm always needs the element with the smallest probability. Also, it keep adding nodes to it, so it would we very inefficient to sort it in every iteration.
    priority_queue<huffmanNode *,vector<huffmanNode *>, decltype(compare)> queue(compare);

    htree=vector<huffmanNode>(events->size()%2?events->size():events->size()-1);

    auto comp2 = [](GLink * e1, GLink * e2){
        return e1->rate>e2->rate;
    };
    priority_queue<GLink *, vector<GLink *>, decltype(comp2)> eventQueue(comp2);
    escape_rate=0.0;

    int firstEmptyFieldIndex=0;
    for (GLink &e:*events){
        eventQueue.push(&e);
        escape_rate+=e.rate;
     } 
    while (eventQueue.size()>1){
        htree[firstEmptyFieldIndex].isOnLastLevel=true;
        htree[firstEmptyFieldIndex].leftLeaf=eventQueue.top();
        eventQueue.pop();
        htree[firstEmptyFieldIndex].rightLeaf=eventQueue.top();
        eventQueue.pop();
        htree[firstEmptyFieldIndex].probability=(htree[firstEmptyFieldIndex].leftLeaf->rate+htree[firstEmptyFieldIndex].rightLeaf->rate)/escape_rate;
        queue.push(&(htree[firstEmptyFieldIndex]));
        firstEmptyFieldIndex++;
    }
    if (!eventQueue.empty()){
        htree[firstEmptyFieldIndex].isOnLastLevel=true;
        htree[firstEmptyFieldIndex].rightLeaf=eventQueue.top();
        htree[firstEmptyFieldIndex].leftLeaf=eventQueue.top();
        htree[firstEmptyFieldIndex].probability=htree[firstEmptyFieldIndex].leftLeaf->rate/escape_rate;
        queue.push(&(htree[firstEmptyFieldIndex]));
        firstEmptyFieldIndex++;
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
        htree[firstEmptyFieldIndex].probability=h1->probability+h2->probability;
        htree[firstEmptyFieldIndex].leftChild=h1;
        htree[firstEmptyFieldIndex].rightChild=h2;
        queue.push(&(htree[firstEmptyFieldIndex]));
        firstEmptyFieldIndex++;
    }
    //reorganize the probabilities: in every node, add the probability of one subtree ("small")
    //to all nodes of the other subtree.
    addProbabilityFromRightSubtreeToLeftSubtree(&htree[htree.size()-1],0);
    moveProbabilitiesFromRightSubtreesOneLevelUp(&htree[htree.size()-1]);
    treeIsMade=true;

}



void huffmanTree::addProbabilityFromRightSubtreeToLeftSubtree(huffmanNode * n,double add){
    //for each node, adds the probability of the right childnode to the left childnode and every node under it. 
    //if the Tree would look like this (with the Numbers representing the probability of every node) before calling this function

    //           1.0
    //       ____||____
    //      |          |
    //     0.4        0.6
    //     _||_      _||_ 
    //    |    |    |    |
    //  0.25 0.15  0.35 0.25
    //        _||_
    //       |    |
    //      0.1  0.05
    //then it would look like this after calling it
    //           1.0
    //       ____||____
    //      |          |
    //     1.0        0.6
    //     _||_      _||_ 
    //    |    |    |    |
    //   1.0 0.75  0.6  0.25
    //        _||_
    //       |    |
    //      0.75 0.65
    //now the tree could be traversed with "while (!n.isLeaf()) n=p>n.right.p?n.left:n.right"
    //so in the function moveProbabilitiesFromRightSubtreesOneLevelUp the numbers are moved one level up to call n.p instead of n.right.p


    //adds "add" to the probability, then calls itself recursively.
    //this calculates the probabilities needed to traverse the tree quickly
    n->probability+=add;
    //if leftId=-1 (=> node is leaf), returns
    if (n->isOnLastLevel){
        return;
    }

    addProbabilityFromRightSubtreeToLeftSubtree(n->leftChild,add+n->rightChild->probability);
    addProbabilityFromRightSubtreeToLeftSubtree(n->rightChild,add);
}



void huffmanTree::moveProbabilitiesFromRightSubtreesOneLevelUp(huffmanNode * n){
    //moves the Probabilities on the right subtrees one level up.
    //if the Tree would look like this (with the Numbers representing the probability of every node) before calling this function
    //           1.0
    //       ____||____
    //      |          |
    //     1.0        0.6
    //     _||_      _||_ 
    //    |    |    |    |
    //   1.0 0.75  0.6  0.25
    //        _||_
    //       |    |
    //      0.75 0.65 
    //then it would look like this after calling it
    //           0.6
    //       ____||____
    //      |          |
    //     0.75      0.25
    //     _||_      _||_ 
    //    |    |    |    |
    //   1.0 0.65  0.6  0.25
    //        _||_
    //       |    |
    //     0.75  0.65
    //note, that now the probabilities on the leaf level are not needed anymore to traverse the tree; the algorithm now is "while (!n.isLeaf()) n=p>n.p?n.left:n.right"
    if (n->isOnLastLevel){
        n->probability-=n->leftLeaf->rate/escape_rate;
    }
    else{
        n->probability=n->rightChild->probability;
        moveProbabilitiesFromRightSubtreesOneLevelUp(n->rightChild);
        moveProbabilitiesFromRightSubtreesOneLevelUp(n->leftChild);
    }
}



GLink * huffmanTree::findHoppingDestination(double p){
    if (!treeIsMade) throw runtime_error("Tried to find Hopping Destination without initializing the Huffmantree first!");
    huffmanNode * node=&htree.back();
    while (!node->isOnLastLevel){
        if (p>node->probability) node=node->leftChild;
        else node=node->rightChild;
    }
    return (p>node->probability?node->leftLeaf:node->rightLeaf);
}

}}
