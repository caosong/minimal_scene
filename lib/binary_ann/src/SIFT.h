/*
 * SIFT.h
 *
 *  Created on: Feb 17, 2011
 *      Author: tt
 */

#ifndef SIFT_H_
#define SIFT_H_

#include <iostream>
#include <vector>
#include <cmath>

class SIFT
{
 public:
  unsigned int size;
  float* descriptor;
  int id;
  static int id_counter;

  SIFT() :
   size(0)
  {
   id = SIFT::id_counter++;
   descriptor = NULL;
  }

  SIFT(float* &dsc, unsigned int s) :
   size(s)
  {
   id = SIFT::id_counter++;
   descriptor = dsc;
  }

  SIFT(float* &dsc, unsigned int s, unsigned int i) :
   id(i)
  {
   size = s;
   descriptor = dsc;
  }

  ~SIFT() {};

  void clear()
  {
   delete descriptor;
   descriptor = NULL;
  }

  void set(float* &dsc, const unsigned int& s, const unsigned int& i)
  {
   size = s;
   descriptor = dsc;
   id = i;
  }

  friend ostream& operator<<(ostream& output, const SIFT);

  float operator^(SIFT& sift) const
  {
   float result = 0,
         diff = 0;

   for (unsigned int i = 0; i < size && i < sift.size; ++i)
   {
    diff = descriptor[i] - sift.descriptor[i];
    result += diff * diff;
   }

   return result;
  }

  bool operator==(const SIFT& sift) const
  {
   return ( id == sift.id );
  }

  bool operator!=(const SIFT& sift) const
  {
   return (  id != sift.id );
  }

  static vector<SIFT> generate_random_samples(unsigned int number_of_random_samples)
  {
   vector<SIFT> dataset;
   unsigned int random_size = 3;
   float* dsc = new float[random_size];
   for (unsigned int i = 0; i < number_of_random_samples; ++i)
   {
    for (unsigned int j = 0; j < random_size; ++j)
     dsc[j] = (float)rand()/(float)RAND_MAX;
    dataset.push_back(SIFT(dsc, random_size));
    for (unsigned int j = 0; j < random_size; ++j)
     dsc[j] = 0.0f;
   }
   delete dsc;
   return dataset;
  }

  static SIFT find_gravity_point(const vector<SIFT>& data)
  {
   return data[rand() % data.size()];
  }


  unsigned int operator&(SIFT sift) const
  {
   return 0;
  }

  int shift(unsigned int i)
  {
   return 0;
  }


};

int SIFT::id_counter = 0;

inline ostream& operator<<(ostream& output, const SIFT bin)
{
 output << "ID: " << bin.id << " ";
 for (unsigned int i = 0; i < bin.size; ++i)
  output << bin.descriptor[i] << " ";
 return output;
}

#endif /* SIFT_H_ */
