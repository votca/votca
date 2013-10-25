/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_KMC_STORE_H_
#define __VOTCA_KMC_STORE_H_

#include<votca/tools/vec.h>

template<class C> class Store {
public:
  Store();
  unsigned int Buy();
  void Sell(unsigned int item_nr);
  void Grow(unsigned int nr_items);
  C* Get_item_by_index(unsigned int index); // Needed for efficient iteration over sold items
  unsigned int Get_itemnr_by_index(unsigned int index); // Needed for efficient iteration over sold items
  C* Get_item(int item_nr);
private:
  std::vector<C> items;
  std::vector<unsigned int> free_numbers;
  std::vector<unsigned int> itemnr_to_index; // Needed for efficient iteration over sold items
  std::vector<unsigned int> index_to_itemnr; // Needed for efficient iteration over sold items
  unsigned int sold_items;
  unsigned int total_items;
};

template<class C> C* Store<C>::Get_item(int item_nr) {
  return &items[item_nr];
}

template<class C> C* Store<C>::Get_item_by_index(unsigned int index) {
  return &items[index_to_itemnr[index]];
}

template<class C> unsigned int Store<C>::Get_itemnr_by_index(unsigned int index) {
  return index_to_itemnr[index];
}

template<class C> Store<C>::Store() {
  sold_items = 0;
  total_items = 0;
}

template<class C> unsigned int Store<C>::Buy() {
  unsigned int free_number = free_numbers[sold_items];
  index_to_itemnr[sold_items] = free_number;
  itemnr_to_index[free_number] = sold_items;
  sold_items++;
  return free_number;
}

template<class C> void Store<C>::Sell(unsigned int item_nr) {
  sold_items--;
  free_numbers[sold_items] = item_nr;
  unsigned int temp1 = itemnr_to_index[item_nr];
  unsigned int temp2 = index_to_itemnr[sold_items];
  index_to_itemnr[temp1] = temp2;
  itemnr_to_index[temp2] = temp1;
}

template<class C> void Store<C>::Grow(unsigned int nr_items) {
  unsigned int new_nr_items = total_items + nr_items;
  free_numbers.resize(new_nr_items);
  items.resize(new_nr_items);
  itemnr_to_index.resize(new_nr_items);
  index_to_itemnr.resize(new_nr_items);
  for (unsigned int i=total_items; i<new_nr_items; i++) {
    free_numbers[i] = i;
  }
  total_items = new_nr_items;
}

#endif