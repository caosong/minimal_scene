/*
 *  Created on: Feb 8, 2011
 *  Author: Tomasz Trzcinski (tomasz.trzcinski@epfl.ch)
 */

#ifndef PARCTREE_H_
#define PARCTREE_H_

#include <vector>
#include <limits.h>

#include "Node.h"
#include "Utils.h"

using namespace std;

template <class T, class D>
class ParcTree {

 private:
  unsigned int branching;

  int balanced_split(const vector<T*>& dataset,
                     vector<unsigned int>& center_indices) const
  {
   // choose dataset image
   vector<unsigned int> image_indices;
   unsigned int image_size = (unsigned int) round(1.0f * (float) dataset.size());
   image_size = (image_size < branching) ? branching : image_size;
   generate_random_indices(dataset.size(), image_size, image_indices);

   // compute average cluster size
   unsigned int avg_cluster_size = round((float) image_indices.size() / (float) branching);
   unsigned int max_disbalance = UINT_MAX,
                temp_disbalance;

   for (unsigned int iteration = 0; iteration < 20; ++iteration)
   {
    // choose random centers
    vector<unsigned int> temp_center_indices;
    generate_random_indices(dataset.size(), branching, temp_center_indices);

    // assign data to the centers
    vector<vector<T> > image_clustered;
    image_clustered.resize(branching);

    unsigned int closest_point,
                 closest_dist,
                 temp_dist;
    for (unsigned int i = 0; i < image_indices.size(); ++i)
    {
     // find the nearest cluster center using exact search
     closest_point = branching;
     closest_dist = UINT_MAX;
     T image_point = *dataset[image_indices[i]];
     for (unsigned int j = 0; j < temp_center_indices.size(); ++j)
     {
      if (image_indices[i] == temp_center_indices[j])
      {
       closest_point = branching;
       break;
      }
      temp_dist = image_point^(*dataset[temp_center_indices[j]]);
      if (temp_dist < closest_dist)
      {
       closest_point = j;
       closest_dist = temp_dist;
      }
     }

     if (closest_point < branching)
      image_clustered[closest_point].push_back(image_point);
    }

    // compute disbalance measure as a variance of cluster sizes
    temp_disbalance = 0;
    for (unsigned int i = 0; i < image_clustered.size(); ++i)
    {
     temp_disbalance += (image_clustered[i].size() - avg_cluster_size)*(image_clustered[i].size() - avg_cluster_size);
    }
    if (temp_disbalance < max_disbalance)
    {
     center_indices.clear();
     center_indices.assign(temp_center_indices.begin(), temp_center_indices.end());
     max_disbalance = temp_disbalance;
    }
   }

   return EXIT_SUCCESS;
  }

  int choose_cluster_centers(const vector<T*>& dataset,
                             vector<unsigned int>& chosen_center_indices) const
  {
   int choose_center = 0;

   switch (choose_center)
   {
    case 0:
     generate_random_indices(dataset.size(), branching, chosen_center_indices);
     break;
    case 1:
     balanced_split(dataset, chosen_center_indices);
     break;
   }

   return EXIT_SUCCESS;
  }

  Node<T>* cluster_data(const vector<T*>& dataset, unsigned int level) const
  {
   if (dataset.empty())
    return NULL;

   Node<T>* temp_node = new Node<T>();
   temp_node->level = level;

   if (dataset.size() < branching)
   {

    temp_node->leaf = true;

    for (unsigned int i = 0; i < dataset.size(); ++i)
     temp_node->bucket.push_back(dataset[i]);

    return temp_node;
   }
   else
   {
    temp_node->leaf = false;

    vector<unsigned int> chosen_indices, dataset_indices;

    choose_cluster_centers(dataset, chosen_indices);

    vector<vector<T*> > dataset_clustered;
    dataset_clustered.resize(branching);

    D closest_point,
      closest_dist,
      temp_dist;
    for (unsigned int i = 0; i < dataset.size(); ++i)
    {
     // find the nearest cluster center using exact search
     closest_point = branching;
     closest_dist = UINT_MAX;

     for (unsigned int j = 0; j < chosen_indices.size(); ++j)
     {
      if (i == chosen_indices[j])
      {
       closest_point = branching;
       break;
      }
      temp_dist = (*dataset[i])^(*dataset[chosen_indices[j]]);
      if (temp_dist < closest_dist)
      {
       closest_point = j;
       closest_dist = temp_dist;
      }
     }

     if (closest_point < branching)
     {
      dataset_clustered[closest_point].push_back(dataset[i]);
     }
    }
    for (unsigned int i = 0; i < branching; ++i)
    {
     if (dataset_clustered[i].empty())
      temp_node->bucket.push_back(dataset[chosen_indices[i]]);
     else
     {
      temp_node->children.push_back(cluster_data(dataset_clustered[i], level + 1));
      temp_node->children.back()->p = *dataset[chosen_indices[i]];
     }
    }

    return temp_node;
   }
  }

  void print_tree() const
  {
   vector<vector<unsigned int> > sorted_nodes;
   sort_nodes(root, sorted_nodes);
   for (unsigned int i = 0; i < sorted_nodes.size(); ++i)
   {
    for (unsigned int j = 0; j < sorted_nodes[i].size(); ++j)
     cout << sorted_nodes[i][j] << " ";
    cout << endl;
   }
  }

 public:

  Node<T>* root;

  ParcTree() :
   root(NULL), branching(0) {};

  ParcTree(const vector<T*>& dataset, unsigned int branching_factor) :
   branching(branching_factor), root(NULL)
  {
   root = cluster_data(dataset, 0);
  }

  ~ParcTree()
  {
   if (root) delete root;
  };

  int search_tree(const T& query, Node<T>* node, vector<T*>& elements, const D& max_dist)
  {
   if (!node)
    return EXIT_SUCCESS;

   D dist_to_node = query^node->p;

   if (node->level != 0 && dist_to_node < max_dist)
     elements.push_back(&(node->p));

   if (node->children.empty())
   {
    if (node->level == 0 && !node->bucket.empty())
         for (unsigned int i = 0; i < node->bucket.size(); ++i)
           elements.push_back((node->bucket[i]));
    return EXIT_SUCCESS;
   }

   vector<Node<T>* > children = node->children;

   D closest_dist = UINT_MAX,
                temp_dist;
   unsigned int closest_point = children.size();

   for (unsigned int i = 0; i < children.size(); ++i)
   {
    temp_dist = query^(children[i]->p);
    if (temp_dist < closest_dist)
    {
     closest_dist = temp_dist;
     closest_point = i;
    }
   }

   if (closest_point < children.size())
   {
    Node<T>* closest_child = children[closest_point];
    if ((closest_dist < max_dist))
     elements.push_back(&(closest_child->p));
    if (!closest_child->bucket.empty())
     for (unsigned int i = 0; i < closest_child->bucket.size(); ++i)
       elements.push_back((closest_child->bucket[i]));

    search_tree(query, closest_child, elements, max_dist);
   }
   return EXIT_SUCCESS;
  }


  friend ostream& operator<<(ostream& output, const ParcTree<T, D>& hrtree)
  {

   if (hrtree.root)
   {
    hrtree.print_tree();
    output << *(hrtree.root) << endl;
   }
   else
    output << "\nEmpty tree." << endl;
   return output;
  }
};

#endif /* PARCTREE_H_ */
