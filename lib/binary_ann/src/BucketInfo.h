/*
 * BucketInfo.h
 *
 *  Created on: May 11, 2011
 *      Author: tt
 */

#ifndef BUCKETINFO_H_
#define BUCKETINFO_H_

class BucketInfo
{
public:
 unsigned int empty_buckets;
 unsigned int mask_index;

 BucketInfo(unsigned int e, unsigned int m) :
  empty_buckets (e), mask_index(m)
 {

 }

 ~BucketInfo() {};

 bool operator<(const BucketInfo& bucket) const
 {
  return ( empty_buckets < bucket.empty_buckets );
 }

};

#endif /* BUCKETINFO_H_ */
