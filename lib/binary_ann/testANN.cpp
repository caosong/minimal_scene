/*
 *  Created on: Feb 8, 2011
 *  Author: Tomasz Trzcinski (tomasz.trzcinski@epfl.ch)
 */

#define EXIT_SUCCESS 0
#define EXIT_FAILURE -1

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <math.h>
#include <time.h>

#include "Index.h"

using namespace std;

int print_results(unsigned int iCorrect,
                  unsigned int iIncorrect,
                  float fRecall,
                  unsigned int nn)
{

#ifdef DEBUG
 printf("\n==============================================\n\n");
 printf("               R E S U L T S :\n");
 printf("\nCorrect 1st NN: %u\tIncorrect: %u\n", iCorrect, iIncorrect);
 printf("Correct 1st NN\tIncorrect\tRecall\n");
#endif
 fprintf(stdout, "%3.2f%%\t", 100.0f * ( (float) iCorrect / (float) (iCorrect + iIncorrect) ));
#ifdef DEBUG

 fprintf(stdout,"%3.2f%%\t", 100.0f * ( (float) iIncorrect / (float) (iCorrect + iIncorrect) ));
#endif
 fprintf(stdout,"%3.5f\t", fRecall);
#ifdef DEBUG
 printf("\n==============================================\n\n");
#endif

 return EXIT_SUCCESS;
}

template <class T>
int check_results(vector<vector<unsigned int> >& ground_truth,
                  vector<vector<T> >& anns,
                  unsigned int nn)
{

 if (ground_truth.empty() || anns.empty())
  return EXIT_FAILURE;

 unsigned int iIncorrect = 0,
              iCorrect = 0;
 float fRecall = 0.0f;
 for (unsigned int i = 0; i < ground_truth.size(); ++i)
 {
  if (!ground_truth[i].empty())
  {
   int iRecall = 0;
   for (unsigned int j = 0; j < ground_truth[i].size() && j < nn; ++j)
   {
    unsigned int k = 0;
    while ( (k < nn) &&
            (ground_truth[i][j] != (unsigned) anns[i][k].id)) ++k;
    if (j == 0)
    {
     if (k == 0)
      ++iCorrect;
     else
      ++iIncorrect;
    }
    if (k < nn)
    {
     ++iRecall;
    }
   }
   if (ground_truth[i].size() < nn)
    fRecall += ((float) iRecall ) / ((float) ground_truth[i].size());
   else
    fRecall += ((float) iRecall ) / ((float) nn);
  }
 }
 fRecall /= ((float) ground_truth.size());

 print_results(iCorrect, iIncorrect, fRecall, nn);

 return EXIT_SUCCESS;
}

void perform_search(vector<Bin128>& dataset,
                    vector<Bin128>& query,
                    string& ground_truth_file,
                    const vector<unsigned int>& params,
                    unsigned int nn)
{
 vector<vector<unsigned int> > ground_truth;
 vector<vector<Bin128> > anns;
 unsigned int nr_descriptors = dataset.size();

 // structure creation
 clock_t begin_createTree = clock();
 Index index(dataset, params, nn);
 float create_time = (clock() - begin_createTree) / (float) CLOCKS_PER_SEC;
#ifdef DEBUG
 printf("\nIndex_hash (%u point(s)) created in %3.2f s\n",
        nr_descriptors,
        create_time);
#endif

//exit(-1);

 // create ground truth
#ifdef DEBUG
 clock_t begin_createGT = clock();
#endif
 ground_truth.resize(nr_descriptors);
 create_ground_truth(ground_truth_file, ground_truth);
#ifdef DEBUG
 printf("Ground Truth read from %s created in %3.2f s\n",
        ground_truth_file.c_str(),
        (clock() - begin_createGT) / (float) CLOCKS_PER_SEC );
#endif

 // search the trees
 vector<unsigned int> adist;
 anns.resize(query.size());
 clock_t begin_query = clock();

 for (unsigned int i = 0; i < query.size(); ++i)
 {
  adist.clear();
  index.search(query[i], anns[i], adist);
 }
 float query_time =  (clock() - begin_query) / (float) CLOCKS_PER_SEC;
#ifdef DEBUG
 printf("\nQuerying took %3.2f\n",
        query_time );
#endif
 // check results against the ground truth
#ifdef DEBUG
 clock_t begin_check = clock();
#endif
 check_results(ground_truth, anns, nn);
#ifdef DEBUG
 printf("Checking results took %3.2f s\n", (clock() - begin_check) / (float) CLOCKS_PER_SEC );
#endif

 fprintf(stdout, "%3.2f\t%3.2f\n", create_time, (float) ( (float) query_time * 1000000.0f / (float) nr_descriptors ) );

 // clearing the elements
 for (unsigned int i = 0; i < dataset.size(); ++i)
  dataset[i].clear();
 dataset.clear();
 query.clear();
 ground_truth.clear();
 for (unsigned int i = 0; i < anns.size(); ++i)
  anns[i].clear();
 anns.clear();
}


int main(int argc, char** argv)
{
 // Check input parameter(s)
 int number_of_arguments = argc - 1;
 if ( number_of_arguments != 5 )
 {
  printf("Wrong number of parameters.\n");
  exit(EXIT_SUCCESS);
 }

#ifdef DEBUG
 // Info screen
 printf("\n================================================================================\n\n");
 printf("     Binary Approximate Nearest-Neighbour Search     \n\n");
 printf("\tAuthor: Tomasz Trzcinski\n\tE-mail: tomasz.trzcinski@epfl.ch\n\n");
 printf("\n================================================================================\n\n");
#endif

 // setting seed random
// srand(time(NULL));
 srand(-21);

 string sInput = (string) (argv[1]);
 string sGroundTruth = (string) (argv[2]);

 // Finding the extension of the input file
 size_t found = sInput.find_last_of(".");
 string sExtension = sInput.substr(found+1);

 unsigned int nn = 2;
 vector<Bin128> dataset, query;
 readBin(sInput, dataset, query);
 vector<unsigned int> params;
 for (int i = 3; i < argc; ++i)
  params.push_back(atoi(argv[i]));
 perform_search(dataset,
                query,
                sGroundTruth,
                params,
                nn);

	return EXIT_SUCCESS;
}
