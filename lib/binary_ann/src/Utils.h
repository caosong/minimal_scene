/*
 *  Created on: Jan 21, 2011
 *  Author: Tomasz Trzcinski (tomasz.trzcinski@epfl.ch)
 */


#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <climits>
#include <string.h>
#include <time.h>
#include <math.h>

#include "Bin128.h"
#include "SIFT.h"

#include "Accumulator.h"

using namespace std;

// COUNTERS
template <class T>
class Counter
{
 public:
  static unsigned int no_dist_computations;
};

template <class T> unsigned int Counter<T>::no_dist_computations = 0;

// READ/WRITE METHODS

int readFloat(string& sInputFile, vector<SIFT>& dataset, vector<SIFT>& query)
{
 ifstream ifData;
 ifData.open(sInputFile.c_str(), ios::in | ios::binary);

 unsigned int nrDim;
 unsigned int iNrDescriptors;
 float  w[128];

 if (ifData.is_open())
 {
  ifData.read((char*)&iNrDescriptors, sizeof(unsigned int));
  ifData.read((char*)&nrDim, sizeof(unsigned int));
#ifdef DEBUG
  cout << "There are " << iNrDescriptors << " float descriptors ";
  cout << "of " << nrDim << " dimensions" << endl;
#endif

  // Creating the dataset
  dataset.resize(iNrDescriptors);
  for(unsigned int i=0; i < iNrDescriptors; ++i)
  {
   float* dsc = new float[nrDim];
   ifData.read((char*)w, sizeof(float)*128);
   for (unsigned int j=0; j < nrDim; j++)
    dsc[j] = w[j];
   dataset[i].set(dsc, nrDim, i);
  }

  // Creating the query
  query.assign(dataset.begin(), dataset.end());

  ifData.close();

  return EXIT_SUCCESS;
 }
 else
 {
  fprintf(stderr,"Unable to open file");
  return EXIT_FAILURE;
 }
}

int readBin(string& sInputFile, vector<Bin128>& dataset, vector<Bin128>& query)
{
 ifstream ifData;
 ifData.open(sInputFile.c_str(), ios::in | ios::binary);
 unsigned int nrDim, iNrDescriptors;
 BIN_WORD w1, w2;

 if (ifData.is_open())
 {
  ifData.read((char*)&iNrDescriptors, sizeof(unsigned int));
  ifData.read((char*)&nrDim, sizeof(unsigned int));
#ifdef DEBUG
  cout << "There are " << iNrDescriptors << " binary descriptors ";
  cout << "of " << nrDim << " dimensions" << endl;
#endif

  dataset.resize(iNrDescriptors);
  for(unsigned int i=0; i < iNrDescriptors; i++)
  {
   ifData.read((char*)&w1, sizeof(BIN_WORD)); //read first 64 bits of the binary descriptor
   ifData.read((char*)&w2, sizeof(BIN_WORD)); //read second 64 bits of the binary descriptor

   dataset[i].set(w1, w2, i);
  }

  // Creating the query
  query.assign(dataset.begin(), dataset.end());

  ifData.close();

  return EXIT_SUCCESS;
 }
 else
 {
  fprintf(stderr,"Unable to open file");
  return EXIT_FAILURE;
 }
}

int create_ground_truth (string sGroundTruth, vector<vector<unsigned int> >& vvGT)
{
 ifstream ifInputFile;
 ifInputFile.open(sGroundTruth.c_str(), ios::in);
 unsigned int iDescriptor;
 string sLine;
 if (ifInputFile.is_open())
 {
  while (ifInputFile.good())
  {
   getline(ifInputFile, sLine);

   char str[200];
   strcpy(str,sLine.c_str());

   // Extraction of Descriptor index
   char *pch;
   pch = strtok (str," ");
   if (pch)
   {
    iDescriptor = atoi(pch);
    if (vvGT.size() < iDescriptor)
    {
     fprintf(stderr, "\nError. Vector vvGT size %zu too small for iDescriptor %d\n", vvGT.size(), iDescriptor);
     exit(EXIT_FAILURE);
    }
    pch = strtok (NULL, " ");
   }
   // Extraction of the neighbours
   while (pch != NULL)
   {
    vvGT[iDescriptor].push_back(atoi(pch));
    pch = strtok (NULL, " ");
   }
  }
 }

 ifInputFile.close();

 return EXIT_SUCCESS;
}

// END OF READ/WRITE METHODS BLOCK

// VECTOR MANIPULATIONS
template <class T>
static bool contains(const vector<T>& v, const T& element)
{
 for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); ++it)
  if (*it == element) return true;
 return false;
}

inline int generate_random_indices(unsigned int max, unsigned int size, vector<unsigned int>& random_indices)
{
 if (size == 0)
  return EXIT_FAILURE;

 random_indices.clear();
 unsigned int num;
 for (unsigned int i = 0; i < size; i++)
 {
  while( contains( random_indices, num = rand() % max ) );
  random_indices.push_back(num);
 }
 sort(random_indices.begin(), random_indices.end());
 return EXIT_SUCCESS;
}

inline unsigned int median(vector<unsigned int> v)
{
 if (v.empty())
  return EXIT_FAILURE;

  sort(v.begin(), v.end());
  unsigned int median_value = (v.size() % 2) ?
                               v[((v.size() + 1) / 2) - 1] :
//                       (float) (v[(v.size() / 2) - 1] + v[(v.size() / 2)]) / 2.0f;
                               v[((v.size() + 1) / 2)];
  return median_value;
}

template <typename T>
inline void print_vector(const vector<T> v)
{
 for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); ++it)
  cout << *it << " ";
}
// END OF VECTOR MANIPULATIONS BLOCK

// SEARCH METHODS (and associated)

template <class T>
inline unsigned int compute_distance(const T& a, const T& b)
{
 ++(Counter<T>::no_dist_computations);
 return a^b;
// return (a > b) ? a - b : b - a;
}

int linearSearch(const Bin128* query,
                 const vector<Bin128*>& dataset,
                 unsigned int nn,
                 vector<Bin128>& nns)
{
 Accumulator<Bin128, unsigned int> acc;
 acc.init(nn);
 for (unsigned int i = 0; i < dataset.size(); ++i)
 {
  Bin128* datapoint = dataset[i];
  if (!contains(acc.elements, *datapoint))
  {
   unsigned int distance = (*query)^(*datapoint);
   if (distance < acc.distances[acc.greatest_distance_index])
   {
    acc.elements[acc.greatest_distance_index] = *datapoint;
    acc.distances[acc.greatest_distance_index] = distance;
    acc.update_greatest_distance();
   }
  }
 }
 vector<unsigned int> dists;
 acc.retrieve_nn(nns, dists);

 return EXIT_SUCCESS;
}
// END OF SEARCH METHODS BLOCK

#endif /* UTILS_H_ */
