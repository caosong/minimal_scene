/*
 * Node.h
 *
 *  Created on: Jan 24, 2011
 *      Author: tt
 */

#ifndef NODE_H_
#define NODE_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

template <class T>
class Node
{
 public:
  T p;
  vector<T*> bucket;
  vector<Node<T>* > children;
  int id;
  unsigned int level;
  bool leaf;

  static unsigned int id_counter;

  Node() :
   p(T()), level(0), leaf(false)
  {
   id = Node<T>::id_counter++;
  }

  ~Node()
  {
   for (unsigned int i = 0; i < children.size(); ++i)
    if (children[i])
     delete children[i];
   bucket.clear();
  }

  friend ostream& operator<<(ostream& output, const Node<T>& node)
  {
   output << endl;
   if (node.leaf)
    output << "Leaf ";
   else
    output << "Node ";
   output << "ID: " << node.id << " Level: " << node.level << " ";
   output << "Pivot: " << node.p << " ";

   if (!node.bucket.empty())
   {
    output << "In bucket: ";
    for (unsigned int i = 0; i < node.bucket.size(); ++i)
    {
     output << node.bucket[i] << " ";
    }
   }

   if (!node.children.empty())
   {
    output << "Children: ";
    for (unsigned int i = 0; i < node.children.size(); ++i)
     output << node.children[i]->id << " ";
    for (unsigned int i = 0; i < node.children.size(); ++i)
     output << *(node.children[i]) << " ";
   }

   return output;
  }

};

template <class T>
unsigned int Node<T>::id_counter = 0;

#endif /* NODE_H_ */
