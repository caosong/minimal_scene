/*
 *  Created on: May 4, 2011
 *  Author: Tomasz Trzcinski (tomasz.trzcinski@epfl.ch)
 */

// TODO: Hashing keys selected so that they do not have overlapping bits
// TODO: No empty hash-buckets as a good-mask condition
// TODO: Have a look at the vectors that are not retrieved correctly

using namespace std;

#include "Utils.h"
#include "Accumulator.h"
#include "BucketInfo.h"

class LshIndex
{
private:
 vector<Bin128*> masks;
 vector<vector<vector<Bin128*> > > data;
 unsigned int nr_datapoints;
 
#ifdef DEBUG
 vector<BucketInfo> buckets;
#endif
  
 // different methods for masks creation
 int createMasks(unsigned int nr_mask_bits, unsigned int nr_masks, const vector<Bin128*>& dataset);
 // random masks initialization
 int createMasks_random(unsigned int nr_mask_bits, unsigned int nr_masks);
 // random masks but with uniform distribution of bits
 int createMasks_overlap(unsigned int nr_mask_bits, unsigned int nr_masks);
 // masks created so that the distribution of 1s resembles distribution of bits in the dataset
 int createMasks_dist(unsigned int nr_mask_bits, unsigned int nr_masks, const vector<Bin128*>& dataset);

 int analyseDataset(const vector<Bin128*>& dataset, vector<float>& distro);
#ifdef DEBUG
 int analyseMasks();
 int analyseData();
 int analyseCorel(const vector<Bin128*>& dataset, unsigned int nr_mask_bits);
#endif

 // method to create the data structure
 int indexData(unsigned int nr_mask_bits, unsigned int nr_masks, const vector<Bin128*>& dataset);

 // method computes index for a given point and a given mask
 unsigned int computeIndex(Bin128* const mask, const Bin128* const datapoint);

public:
 // default constructor
 LshIndex() {};
 
 // constructing data structure from the dataset
 LshIndex(const vector<Bin128*>& dataset,
            unsigned int nr_masks,
            unsigned int nr_mask_bits);
 
 // destructor
 ~LshIndex();
 
 // search method
 // arg0: query item
 // arg1: accumulator
 int search(const Bin128& query, Accumulator<Bin128, unsigned int>& acc);


#ifdef DEBUG
 int print_indices(Bin128 datapoint);
#endif

};			

LshIndex::LshIndex(const vector<Bin128*>& dataset,
                       unsigned int nr_masks,
                       unsigned int nr_mask_bits)
{
 nr_datapoints = dataset.size();
#ifdef DEBUG
 printf("\nAnalysing dataset...\n");
 vector<float> distro;
 analyseDataset(dataset, distro);
#endif
 masks.resize(nr_masks);
 createMasks(nr_mask_bits, nr_masks, dataset);

#ifdef DEBUG
 printf("\nAnalysing masks...\n");
 analyseMasks();
#endif

 data.resize(nr_masks);
 indexData(nr_mask_bits, nr_masks, dataset);

#ifdef DEBUG
 printf("\nAnalysing data structure...\n");
 analyseData();
// analyseCorel(dataset, nr_mask_bits);
#endif
}

LshIndex::~LshIndex()
{
 for (vector<Bin128*>::const_iterator it = masks.begin();
      it != masks.end();
      ++it)
  delete *it;
 masks.clear();

 for (unsigned int i = 0; i < data.size(); ++i)
 {
  for (unsigned int j = 0; j < data[i].size(); ++j)
  {
   data[i][j].clear();
  }
  data[i].clear();
 }
 data.clear();

#ifdef DEBUG
 buckets.clear();
#endif
}

int LshIndex::createMasks(unsigned int nr_mask_bits,
                            unsigned int nr_masks,
                            const vector<Bin128*>& dataset)
{

 unsigned int choose = 1;
 switch (choose)
 {
     case 0:
      createMasks_random(nr_mask_bits, nr_masks); break;
     case 1:
      createMasks_overlap(nr_mask_bits, nr_masks); break;
     case 2:
      createMasks_dist(nr_mask_bits, nr_masks, dataset); break;
 }
 return EXIT_SUCCESS;
}

int LshIndex::createMasks_random(unsigned int nr_mask_bits,
                                   unsigned int nr_masks)
{
 Bin128* mask;
 vector<unsigned int> random_indices,
                      masks_stats(DIM, 0);
 unsigned int shift;
 for (unsigned int i = 0; i < nr_masks; ++i)
 {
  generate_random_indices(DIM - 1,
                          nr_mask_bits,
                          random_indices);
  mask = new Bin128();
  vector<bool> mask_bit (DIM, false);
  for (unsigned int j = 0; j < nr_mask_bits; ++j)
  {
   shift = random_indices[j];
   if (shift < (DIM_2))
    mask->w1 |= (ONE << shift);
   else
    mask->w2 |= (ONE << (shift - (DIM_2)));
  }
  masks[i] = mask;
#ifdef DEBUG
    cout << *mask << "\tnr of 1s: " << __builtin_popcountll(mask->w1) + __builtin_popcountll(mask->w2) << endl;
#endif
 }
 return EXIT_SUCCESS;
}

int LshIndex::createMasks_overlap(unsigned int nr_mask_bits,
                                       unsigned int nr_masks)
{
 Bin128* mask;
 vector<unsigned int> random_indices,
                      masks_stats(DIM, 0);
 float avg_bit = (float) (nr_mask_bits*nr_masks) / (float) DIM;;
 unsigned int shift;
 for (unsigned int i = 0; i < nr_masks; ++i)
 {
  generate_random_indices(DIM - 1,
                          nr_mask_bits,
                          random_indices);
  mask = new Bin128();

  vector<bool> mask_bit (DIM, false);
  for (unsigned int j = 0; j < nr_mask_bits; ++j)
  {
   float max_bit = avg_bit;
   shift = random_indices[j];
   if ( (float) masks_stats[shift] > max_bit || mask_bit[shift])
   {
    unsigned int min = nr_masks,
                 min_id = DIM;
    for (unsigned int m = 0; m < DIM; ++m)
     if (masks_stats[m] < min)
     {
      min = masks_stats[m];
      min_id = m;
     }
    shift = (min_id != DIM) ? min_id : (shift + 1) % DIM;
   }
   ++masks_stats[shift];
   mask_bit[shift] = true;
   if (shift < (DIM_2))
    mask->w1 |= (ONE << shift);
   else
    mask->w2 |= (ONE << (shift - (DIM_2)));
  }
  masks[i] = mask;
#ifdef DEBUG
  //  cerr << "Mask creation took: " <<  (clock() - begin_mask) / (float) CLOCKS_PER_SEC << " secs" << endl;
    cout << *mask << "\tnr of 1s: " << __builtin_popcountll(mask->w1) + __builtin_popcountll(mask->w2) << endl;
#endif
 }
 return EXIT_SUCCESS;
}

int LshIndex::createMasks_dist(unsigned int nr_mask_bits,
                                 unsigned int nr_masks,
                                 const vector<Bin128*>& dataset)
{
 vector<float> distro;
 analyseDataset(dataset, distro);

 Bin128* mask;
 vector<unsigned int> random_indices,
                      masks_stats(DIM, 0);
 unsigned int shift;
 for (unsigned int i = 0; i < nr_masks; ++i)
 {
  generate_random_indices(DIM - 1,
                          nr_mask_bits,
                          random_indices);
  mask = new Bin128();

  vector<bool> mask_bit (DIM, false);
  for (unsigned int j = 0; j < nr_mask_bits; ++j)
  {
   float max_bit = distro[j] * (float) nr_mask_bits;
   shift = random_indices[j];
   if ( (float) masks_stats[shift] > max_bit || mask_bit[shift])
   {
    unsigned int min = nr_masks,
                 min_id = DIM;
    for (unsigned int m = 0; m < DIM; ++m)
     if (masks_stats[m] < min)
     {
      min = masks_stats[m];
      min_id = m;
     }
    shift = (min_id != DIM) ? min_id : (shift + 1) % DIM;
   }
   ++masks_stats[shift];
   mask_bit[shift] = true;
   if (shift < (DIM_2))
    mask->w1 |= (ONE << shift);
   else
    mask->w2 |= (ONE << (shift - (DIM_2)));
  }
  masks[i] = mask;
#ifdef DEBUG
    cout << *mask << "\tnr of 1s: " << __builtin_popcountll(mask->w1) + __builtin_popcountll(mask->w2) << endl;
#endif
 }
 return EXIT_SUCCESS;
}

int LshIndex::indexData(unsigned int nr_mask_bits, unsigned int nr_masks, const vector<Bin128*>& dataset)
{
 Bin128* mask;
 for (unsigned int i = 0; i < nr_masks; ++i)
 {
  mask = masks[i];
  BIN_WORD max_index = ONE << (__builtin_popcountll(mask->w1) + __builtin_popcountll(mask->w2));
  data[i].resize(max_index);
  for (unsigned int d = 0;
       d < dataset.size();
       ++d)
  {
   data[i][computeIndex(mask,dataset[d])].push_back(dataset[d]);
  }
 }

 return EXIT_SUCCESS;
}

int LshIndex::analyseDataset(const vector<Bin128*>& dataset,
                               vector<float>& distro)
{
 // iterate over all masks
 vector<unsigned int> dataset_stats(DIM, 0);
 distro.resize(DIM);
 BIN_WORD b,
           temp = (BIN_WORD) 0;
 for (unsigned int i = 0; i < DIM; ++i)
 {
  b = ONE << (i % DIM_2);
  for (unsigned int j = 0; j < dataset.size(); ++j)
  {
   temp = (i < 64) ? (b & dataset[j]->w1) : (b & dataset[j]->w2);
   if (temp)
    ++dataset_stats[i];
  }
 }

 float sum = 0.0f;
 printf("Distribution of 1's along the dataset\n");
 for (int i = DIM_2 - 1; i >= 0; --i)
 {
  distro[i] = (float) dataset_stats[i] * 100.0f / (float) nr_datapoints;
  printf("%3.2f%% ", distro[i]);
  sum += (float) dataset_stats[i];
 }
 for (int i = DIM - 1; i >= DIM_2; --i)
 {
  distro[i] = (float) dataset_stats[i] * 100.0f / (float) nr_datapoints;
  printf("%3.2f%% ",distro[i]);
  sum += (float) dataset_stats[i];
 }
 printf("\tAverage: %3.2f\n", sum * 100.0f / (float) (DIM*nr_datapoints));

 return EXIT_SUCCESS;
}

unsigned int LshIndex::computeIndex(Bin128* const mask, const Bin128* const datapoint)
{
 unsigned int index = 0;
 int mask_counter = 0;
 BIN_WORD b;
 for (int i = 0; i < DIM; ++i)
 {
  if (i < DIM_2)
  {
   b = ONE << i;
   if (mask->w1 & b)
   {
    if (datapoint->w1 & b)
     index += ONE << mask_counter;
    ++mask_counter;
   }
  }
  else
  {
   b = ONE << (i - DIM_2);
   if (mask->w2 & b)
   {
    if (datapoint->w2 & b)
     index += ONE << mask_counter;
    ++mask_counter;
   }
  }
 }
 return index;
}

int LshIndex::search(const Bin128& query,
                     Accumulator<Bin128, unsigned int>& acc)
{
 for (unsigned int i = 0; i < masks.size(); ++i)
 {
  unsigned int mask_index = i;
  unsigned int query_index = computeIndex(masks[mask_index], &query);
  for (unsigned int j = 0; j < data[mask_index][query_index].size(); ++j)
  {
   Bin128* datapoint = data[mask_index][query_index][j];
   if (!contains(acc.elements, *datapoint) && (*datapoint != query))
   {
    unsigned int distance = query^(*datapoint);
    if (distance < acc.distances[acc.greatest_distance_index])
    {
     acc.elements[acc.greatest_distance_index] = *datapoint;
     acc.distances[acc.greatest_distance_index] = distance;
     acc.update_greatest_distance();
    }
   }
  }
 }

 return EXIT_SUCCESS;
}

#ifdef DEBUG

int LshIndex::analyseMasks()
{
 vector<unsigned int> stats (DIM, 0);
 BIN_WORD b,
          temp = (BIN_WORD) 0;
 for (unsigned int i = 0; i < DIM; ++i)
 {
  b = ONE << (i % DIM_2);
  for (unsigned int j = 0; j < masks.size(); ++j)
  {
   temp = (i < 64) ? (b & masks[j]->w1) : (b & masks[j]->w2);
   if (temp)
    ++stats[i];
  }
 }

 float sum = 0.0f;
 printf("Distribution of 1's along the masks\n");
 for (int i = DIM_2 - 1; i >= 0; --i)
 {
  printf("%u",stats[i]);
  sum += (float) stats[i];
 }
 printf("\t\t\t");
 for (int i = DIM - 1; i >= DIM_2; --i)
 {
  printf("%u",stats[i]);
  sum += (float) stats[i];
 }
 printf("\tAverage: %3.2f\n", sum / (float) DIM);

 return EXIT_SUCCESS;
}

int LshIndex::analyseData()
{
 // iterate over all masks
 float avg_bucket_size = 0.0f;
 for (unsigned int i = 0; i < data.size(); ++i)
 {
  cout << "Mask " << i << " with " << data[i].size() << " indices: ";
  // iterate over all indices
  int empty_buckets = 0;
//  float avg_bucket_size = nr_datapoints / data[i].size();
//  float variance = 0.0f;
  for (unsigned int j = 0; j < data[i].size(); ++j)
  {
   if (data[i][j].empty()) ++empty_buckets;
//   float temp = (float) data[i][j].size() - avg_bucket_size;
//   variance += temp * temp;
  }
  buckets.push_back(BucketInfo(empty_buckets, i));
  avg_bucket_size += (float) nr_datapoints / (float) (data[i].size() - empty_buckets);
  cout << empty_buckets << " empty buckets, avg_bucket_size = " << (float) nr_datapoints / (float) (data[i].size() - empty_buckets) << endl; // variance = " << sqrt(variance) << endl;
 }

 printf("Sorting buckets...\n");
 sort(buckets.begin(), buckets.end());

 for (unsigned int i = 0; i < buckets.size(); ++i)
  printf("Mask %u has %u empty buckets.\n", buckets[i].mask_index, buckets[i].empty_buckets);

// printf("Average bucket size = %3.2f, avg_number_distance_computations = %3.2f\n", avg_bucket_size / (float) data.size(), avg_bucket_size);
// fprintf(stderr,"%3.2f\t%3.2f\t%3.2f\n", avg_bucket_size / (float) data.size(), avg_bucket_size, (float) nr_datapoints / avg_bucket_size);

 return EXIT_SUCCESS;
}

int LshIndex::analyseCorel(const vector<Bin128*>& dataset, unsigned int nr_mask_bits)
{
 ofstream output("output.txt", ios::out | ios::app);

 unsigned int* sigma = new unsigned int[nr_datapoints];

 float total_common_dsc = 0.0f;
 float total_common_bck = 0.0f;

// fprintf(stderr,"Created correlation matrix %d x %d\t", nr_datapoints, nr_datapoints);

 for (unsigned int i = 0; i < dataset.size(); ++i)
 {
  // clearing the array
  for (unsigned int j = 0; j < nr_datapoints; ++j)
   sigma[j] = 0;
  fprintf(stderr,"\r %d/%d descriptors analysed.", i+1, nr_datapoints);
  for (unsigned int j = 0; j < masks.size(); ++j)
  {
//   fprintf(stderr,"j = %d ", j);
//   cerr << "Mask[" << j << "] = " << *masks[j];
//   if (dataset[i])
//    cerr << " dataset[" << i << "] = " << *dataset[i];
   unsigned int index = computeIndex(masks[j], dataset[i]);
//   fprintf(stderr, "Computed index: %d\t",index);
   for (unsigned int k = 0; k < data[j][index].size(); ++k)
   {
////    fprintf(stderr,"k = %d", k);
//    cerr << data[j][index][k]->id << " ";
//    if ((i*nr_datapoints + data[j][index][k]->id) > (nr_datapoints * nr_datapoints))
//     cerr << "i = " << i << " id = " << data[j][index][k]->id << "\t";
    ++sigma[data[j][index][k]->id];
   }
  }
  float common_dsc = 0.0f;
  float common_bck = 0.0f;
  for (unsigned int col = 0; col < nr_datapoints; ++col)
  {
   if ( i != col && sigma[col] ) // && ((float) sigma[col] / (float) sigma[i]) > 0.2f )
   {
    common_dsc += 1.0f;
    common_bck += (float) sigma[col];
//    output << col << "_" << sigma[col] << " ";
   }
  }
  total_common_dsc += common_dsc;
  total_common_bck += (common_bck / common_dsc);
//  output << endl;
 }

 fprintf(stderr, "Avg. common descriptors: %3.2f\t Avg. common buckets: %3.2f\n",
         total_common_dsc / (float) nr_datapoints,
         total_common_bck / (float) nr_datapoints);

 output << masks.size() << "\t" << nr_mask_bits << "\t";
 output <<  total_common_dsc / (float) nr_datapoints << " ";
 output << total_common_bck / (float) nr_datapoints << endl;

// ofstream output("output.txt", ios::out);
// for (unsigned int row = 0; row < nr_datapoints; ++row)
// {
//  for (unsigned int col = 0; col < nr_datapoints; ++col)
//  {
//   if (sigma[row][col])
//    output << col << "_" << sigma[row][col] << " ";
//  }
//  output << endl;
// }
 output.close();

 delete [] sigma;
 return EXIT_SUCCESS;
}

int LshIndex::print_indices(Bin128 datapoint)
{
 vector<unsigned int> indices;
 for (unsigned int i = 0; i < masks.size(); ++i)
  indices.push_back(computeIndex(masks[i], &datapoint));

 sort(indices.begin(), indices.end());

 for (unsigned int i = 0; i < indices.size(); ++i)
  printf("%u ", indices[i]);
 return EXIT_SUCCESS;
}
#endif
