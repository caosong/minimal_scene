//============================================================================
// Name        : extract_sift.cpp
// Author      : Tomasz Trzcinski (tomasz.trzcinski@epfl.ch)
//============================================================================

#define EXIT_SUCCESS 0
#define EXIT_FAILURE -1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <vector>

#include "hashpro.h"

using namespace std;

int readSIFTdescriptors(string inputDesc, vector<vector<float> >& descriptors)
{
 string line;
 ifstream inputData(inputDesc.c_str());
 if (inputData.is_open())
 {
  getline (inputData, line);
  unsigned int maxIndex = atoi(line.c_str());
 //  cout << "Max index: " << maxIndex << endl;
  descriptors.resize(maxIndex+1);
  size_t found;
  while (inputData.good())
  {
   getline (inputData, line);
   if (strcmp(line.c_str(),"") == 0) continue;
   found = line.find(" ");
   unsigned int index = atoi(line.substr(0,found).c_str());
   string desc = line.substr(found+1);
//    cout << "\nIndex: " << index << " desc: ";
   char * pch = strtok ((char*) desc.c_str()," ");
   while (pch != NULL)
   {
    descriptors[index].push_back( atof(pch));
//     cout << atoi(pch) << " ";
    pch = strtok (NULL, " ");
   }
  }
  inputData.close();
 }
 else
 {
  fprintf(stderr, "Unable to open file");
  return EXIT_FAILURE;
 }

 return EXIT_SUCCESS;
}

//void sseg_matrix_vector_mul(const float* A, int ar, int ac, int ald, const float* b, float* c)
//{
//  for( int r=0; r<ar; r++ )
//  {
//    c[r] = 0.0;
//    for(int j=0; j < ac; j++)
//    {
//      c[r] += A[r*ald+j]*b[j];
//    }
//  }
//}

void multiplySIFTbyMatrix(vector<float>& sift, vector<float>& siftMultiplied)
{
 for (unsigned int r = 0; r < 128; ++r)
 {
  siftMultiplied[r] = 0.0f;
  for (unsigned int c = 0; c < 128; ++c)
  {
   siftMultiplied[r] += Adif128[r*128+c]*sift[c];
  }
 }
}

int computeLDAdescriptors(vector<vector<float> >& siftDescriptors, vector<vector<bool> >& ldaDescriptors)
{
 for (unsigned int i = 0; i < siftDescriptors.size(); ++i)
 {
  if (!siftDescriptors[i].empty())
  {
   vector<float> siftMultiplied, siftDescriptor;
   siftMultiplied.resize(siftDescriptors[i].size());
   siftDescriptor.assign(siftDescriptors[i].begin(), siftDescriptors[i].end());
   for (unsigned int d = 0; d < siftDescriptor.size(); ++d)
    siftDescriptor[d] /= 512.f;
   multiplySIFTbyMatrix(siftDescriptor, siftMultiplied);
   for (unsigned int j = 0; j < siftMultiplied.size(); ++j)
    if(siftMultiplied[j] + tdif128[j] <= 0.0)
     ldaDescriptors[i].push_back(true);
    else
     ldaDescriptors[i].push_back(false);
  }
 }

 return EXIT_SUCCESS;
}

int saveLDAdescriptors(vector<vector<bool> >& ldaDescriptors, const string& output)
{
 ofstream outputData(output.c_str());
 outputData << ldaDescriptors.size() - 1;
 for (unsigned int i = 0; i < ldaDescriptors.size(); ++i)
 {
  if (!ldaDescriptors[i].empty())
  {
   outputData << endl << i << " ";
   for (unsigned int j = 0; j < ldaDescriptors[i].size(); ++j)
    if (ldaDescriptors[i][j])
     outputData << "1";
    else
     outputData << "0";
  }
 }
 outputData.close();

 return EXIT_SUCCESS;
}

int main(int argc, char** argv)
{
 if( argc != 3 )
 {
  cout << "Wrong number of parameters." << endl;
  cout << "SYNTAX: " << argv[0] << " InputSIFTdescriptorsFile OutputLDAHashDescriptorsFile" << endl;
  return EXIT_FAILURE;
 }

 vector<vector<float> > siftDescriptors;
 readSIFTdescriptors((string) argv[1], siftDescriptors);

 vector<vector<bool> > ldaDescriptors;
 ldaDescriptors.resize(siftDescriptors.size());
 computeLDAdescriptors(siftDescriptors, ldaDescriptors);

 saveLDAdescriptors(ldaDescriptors, (string) argv[2]);

 return EXIT_SUCCESS;
}
