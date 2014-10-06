/*
 *  Created on: Feb 8, 2011
 *  Author: Tomasz Trzcinski (tomasz.trzcinski@epfl.ch)
 */

#ifndef INDEX_H_
#define INDEX_H_

#include <vector>
#include <iostream>

#include "ParcIndex.h"
#include "LshIndex.h"
#include "Utils.h"
#include "Accumulator.h"

using namespace std;

class Index {
 private:

  // pointers to the subsequent data structures
  ParcIndex* parcIndex;
  LshIndex* lshIndex;

  // number of nearest neighbors to be retrieved
  unsigned int nn;

  // accumulator used to store and retrieve nearest neighbors
  Accumulator<Bin128, unsigned int> acc;

 public:
  // default constructor
  Index() {};

  // parametrized constructor

  // index parameters passed using params vector
  // for LSH:
  // # params[0]: 0 (index type)
  // # params[1]: number of keys used
  // # params[2]: number of bits used in each key
  // # Decent results (0.9 precision @ first place) should be obtained for 40 keys of 14 bits
  // for PARC trees:
  // # params[0]: 1 (index type)
  // # params[1]: number of trees
  // # params[2]: branching factor of the trees
  // # Decent results (0.9 precision @ first place) should be obtained for 30 trees with 64 branching factor
  Index(vector<Bin128>& dataset,
        const vector<unsigned int>& params,
        const unsigned int& in_nn);

  ~Index()
  {
   if (lshIndex)
    delete lshIndex;
   if (parcIndex)
    delete parcIndex;
  };

  // method to retrieve nearest neighbors from the accumulator.
  // anns is a vector of sorted anns (from the nearest to the
  // furthest) and distances is the list of the corresponding
  // distances
  int search(const Bin128& query,
             vector<Bin128>& anns,
             vector<unsigned int>& distances);

};

Index::Index(vector<Bin128>& dataset,
             const vector<unsigned int>& params,
             const unsigned int& in_nn) :
             nn(in_nn)
{
 // create pointers to the data to save space
 vector<Bin128*> dataset_ptrs;
 dataset_ptrs.resize(dataset.size());
 for (unsigned int i = 0; i < dataset.size(); ++i)
  dataset_ptrs[i] = &dataset[i];

 parcIndex = NULL;
 lshIndex = NULL;
 switch(params[0])
 {
  // create LSH index
  case 0:
   lshIndex = new LshIndex(dataset_ptrs, params[1], params[2]); break;
  // create Parc trees index
  case 1:
   parcIndex = new ParcIndex(dataset_ptrs, params[1], params[2]); break;
 }
}

int Index::search(const Bin128& query,
                  vector<Bin128>& anns,
                  vector<unsigned int>& distances)
{
 acc.init(nn);
 if (lshIndex)
  lshIndex->search(query, acc);
 else if (parcIndex)
  parcIndex->search(query, acc);

 acc.retrieve_nn(anns, distances);

 return EXIT_SUCCESS;
}

#endif /* INDEX_H_ */
