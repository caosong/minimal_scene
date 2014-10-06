/*
 *  Created on: Jan 25, 2011
 *  Author: Tomasz Trzcinski (tomasz.trzcinski@epfl.ch)
 */

#ifndef ACCUMULATOR_H_
#define ACCUMULATOR_H_

#include <vector>

using namespace std;

template <class T, class D>
class Accumulator {
public:
 Accumulator() {};
 ~Accumulator() {};

 vector<T> elements;
 vector<D> distances;
 unsigned int greatest_distance_index;

 int init(unsigned int nn)
 {
  elements.clear();
  distances.clear();
  elements.resize(nn);
  for (unsigned int i = 0; i < nn; ++i)
   distances.push_back(UINT_MAX);
  update_greatest_distance();

  return EXIT_SUCCESS;
 }

 int update_greatest_distance()
 {
  D greatest_temp_distance = 0;
  for (unsigned int it = 0; it < distances.size(); ++it)
  {
   if (distances[it] >= greatest_temp_distance )
   {
    greatest_distance_index = it;
    greatest_temp_distance = distances[it];
   }
  }
  return EXIT_SUCCESS;
 }

 int retrieve_nn(vector<T>& anns,
                 vector<D>& dists) const
 {
  anns.assign(elements.begin(), elements.end());
  dists.assign(distances.begin(), distances.end());

  // sort anns
  bool swapped;
  do
  {
   swapped = false;
   for (unsigned int it = 1; it < dists.size(); ++it)
   {
    if (dists[it - 1] > dists[it])
    {
     D temp = dists[it];
     dists[it] = dists[it - 1];
     dists[it - 1] = temp;
     T temp2 = anns[it];
     anns[it] = anns[it - 1];
     anns[it - 1] = temp2;
     swapped = true;
    }
   }
  }
  while (swapped);

  // cut off the wrong anns
  unsigned int anns_end = 0;
  while (anns_end < dists.size() && dists[anns_end] < UINT_MAX) ++anns_end;
  anns.erase(anns.begin() + anns_end, anns.end());
  dists.erase(dists.begin() + anns_end, dists.end());

  return EXIT_SUCCESS;
 }

};

#endif /* ACCUMULATOR_H_ */
