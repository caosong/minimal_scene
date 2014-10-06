/*
 *
 *  Created on: Feb 8, 2011
 *  Author: Tomasz Trzcinski (tomasz.trzcinski@epfl.ch)
 */

#ifndef PARCINDEX_H_
#define PARCINDEX_H_

#include "Utils.h"
#include "ParcTree.h"
#include "Accumulator.h"

using namespace std;

class ParcIndex {
 private:
  // number of trees
  unsigned int nr_trees;
  // branching factor of the trees
  unsigned int branch;

 public:
  // forest of the trees
  vector<ParcTree<Bin128, unsigned int>* > trees;

  // default constructor
  ParcIndex() {};

  // parametrized constructor
  ParcIndex(const vector<Bin128*>& dataset,
            unsigned int in_nr_trees,
            unsigned int in_branch) :
             nr_trees(in_nr_trees),
             branch(in_branch)
  {
   trees.resize(nr_trees);
   for (unsigned int i = 0; i < nr_trees; ++i)
    trees[i] = new ParcTree<Bin128, unsigned int>(dataset, branch);
  }

  // default destructor
  ~ParcIndex()
  {
   for (unsigned int i = 0; i < trees.size(); ++i)
    if (trees[i])
     delete trees[i];
   trees.clear();
  };

  // method searches the nearest neighbor of a query
  // (passed as a reference) and returns n nearest
  // neighbors (unsorted, stored in the accumulator
  int search(const Bin128& query, Accumulator<Bin128,unsigned int>& acc);
};

int ParcIndex::search(const Bin128& query,
                      Accumulator<Bin128,unsigned int>& acc)
{
 vector<Bin128*> anns;

 for (unsigned int i = 0; i < trees.size(); ++i)
 {
  anns.clear();
  trees[i]->search_tree(query,
                        trees[i]->root,
                        anns,
                        acc.distances[acc.greatest_distance_index]);
  for (unsigned int j = 0; j < anns.size(); ++j)
  {
   if (!contains(acc.elements, *(anns[j])) && (query != *anns[j]) )
   {
    unsigned int distance = query^(*(anns[j]));
    if (distance < acc.distances[acc.greatest_distance_index])
    {
     acc.elements[acc.greatest_distance_index] = *(anns[j]);
     acc.distances[acc.greatest_distance_index] = distance;
     acc.update_greatest_distance();
    }
   }
  }
 }
 return EXIT_SUCCESS;
}

#endif /* PARCINDEX_H_ */
