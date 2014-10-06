/*
 * Bin128.h
 *
 *  Created on: Jan 25, 2011
 *      Author: tt
 */

#ifndef BIN128_H_
#define BIN128_H_

#define ONE  ((BIN_WORD)(1))
#define BIN_WORD unsigned long long
#define DIM 128
#define DIM_2 64

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Bin128
{
 public:
  BIN_WORD w1, w2;
  int id;
  static int id_counter;

  Bin128() : w1(0), w2(0)
  {
   id = Bin128::id_counter++;
  }

  Bin128(BIN_WORD word1, BIN_WORD word2) :
   w1(word1), w2(word2)
  {
   id = Bin128::id_counter++;
  }

  Bin128(BIN_WORD word1, BIN_WORD word2, int identifier) :
   w1(word1), w2(word2), id(identifier) { }

  ~Bin128() {};

  void clear() {};

  void set(BIN_WORD word1, BIN_WORD word2, int identifier)
  {
   w1 = word1;
   w2 = word2;
   id = identifier;
  }

  friend ostream& operator<<(ostream& output, const Bin128);
  friend Bin128 operator~(const Bin128& binary);

  Bin128 XOR(const Bin128& binary) const
  {
   return Bin128(w1^binary.w1, w2^binary.w2);
  }

  unsigned int operator^(const Bin128& binary) const
  {
   return ( __builtin_popcountll(w1^binary.w1) + __builtin_popcountll(w2^binary.w2));
  }

  unsigned int operator&(const Bin128& binary) const
  {
   return ( (w1 & binary.w1) + (w2 & binary.w2));
  }

  Bin128 operator|(const Bin128& binary) const
  {
   return Bin128( (w1 | binary.w1), (w2 | binary.w2));
  }

  bool operator==(const Bin128& binary) const
  {
   return ( id == binary.id );
  }

  bool operator!=(const Bin128& binary) const
  {
   return (  id != binary.id );
  }

  static Bin128 find_gravity_point(const vector<Bin128>& data)
  {
   Bin128 gravity_point;
   unsigned int gravity[128] = {0};
   BIN_WORD b;
   for (unsigned int i = 0; i < data.size(); ++i)
   {
    for(int l=0; l < 64; ++l)
    {
        b = ONE << l;
        if((data[i].w1 & b) == b) ++gravity[l];
        if((data[i].w2 & b) == b) ++gravity[l+64];
    }
   }

   for (unsigned int i = 0; i < 128; ++i)
    gravity[i] = (unsigned int) round((float) gravity[i] / (float) data.size());

   for(int l=0; l < 64; ++l)
   {
    b = ONE << l;
    if (gravity[l])
     gravity_point.w1 += b;
    if (gravity[l+64])
     gravity_point.w2 += b;
   }

   return gravity_point;
  }

  static vector<Bin128> generate_random_samples(unsigned int number_of_random_samples)
  {
   vector<Bin128> dataset;
   for (unsigned int i = 0; i < number_of_random_samples; ++i)
   {
    BIN_WORD w1( rand() ),
             w2( rand() );
    w1 <<= 32; w2 <<= 32;
    w1 += rand(); w2 += rand();
    dataset.push_back(Bin128(w1, w2));
   }
   return dataset;
  }

  int shift(unsigned int i)
  {
   w1 = 0;
   w2 = 0;
   if (i <= 64)
   {
    w1 += 1;
    w1 <<= i;
   }
   else
   {
    w2 += 1;
    w2 <<= (i-64);
   }
   return EXIT_SUCCESS;
  }
};

int Bin128::id_counter = 0;

inline ostream& operator<<(ostream& output, const Bin128 bin)
{
 BIN_WORD b;
 for(int l=63; l >= 0; --l)
 {
  b = ONE << l;
  if((bin.w1 & b) == b) output << "1";
  else              output << "0";
 }
 output << " (" << bin.w1 << ")\t";
 for(int l=63; l >= 0; --l)
 {
  b = ONE << l;
  if((bin.w2 & b) == b) output << "1";
  else              output << "0";
 }
 output << " (" << bin.w2 << ")";
 return output;
}

Bin128 operator~(const Bin128& binary)
{
 return Bin128(~binary.w1,~binary.w2);
}

#endif /* BIN128_H_ */
