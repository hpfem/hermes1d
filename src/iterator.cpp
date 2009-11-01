// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "iterator.h"

void Iterator::reset() {
  current_coarse_elem_index = -1;
  // clear the stack
  while(!S.empty()) S.pop();
}

Element* Iterator::first_active_element()
{
  return this->mesh->first_active_element();
}

Element* Iterator::next_active_element()
{
  Element *e;
  // find the first element that has not been visited
  if(current_coarse_elem_index == -1) { //either it is the first coarse mesh element
    e = this->mesh->get_base_elems();
    current_coarse_elem_index = 0; 
  }
  else { // or we take it from the stack
    if(S.empty()) {
      if(current_coarse_elem_index == this->mesh->get_n_base_elem()-1) { // there is no new element to visit
        return NULL;
      }
      else { // we take the next coarse mesh element
        e = this->mesh->get_base_elems() + current_coarse_elem_index +1;
        current_coarse_elem_index++; 
      }
    }
    else { // if stack is not empty
      e = S.top();     // take the top element
      S.pop(); // remove it from stack
    }
  }
  // if this element is active, return it
  if(e->is_active()) return e;

  // descending to first active element
  while(!e->is_active()) {
    S.push(e->sons[1]);
    e = e->sons[0];
  }
  return e;
}

Element* Iterator::last_active_element()
{
  return this->mesh->last_active_element();
}




