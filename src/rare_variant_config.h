#ifndef RARE_IBS0_H
#define RARE_IBS0_H

#include "cramore.h"
#include "bcf_ordered_reader.h"
#include "bp2cm.h"
#include "utils.h"

struct Node {

  Node *left, *right;
  std::vector<int32_t> leaves;
  int32_t height = 0;

  Node(int32_t x) {
    left = NULL; right = NULL;
    leaves.push_back(x);
  }

  Node(Node *_l, Node *_r, int32_t _h) : height(_h) {
    if (_l->height > _r->height) { // Left always has lower height
      left = _r; right = _l;
    } else {
      left = _l; right = _r;
    }
    leaves.reserve(_l->leaves.size() + _r->leaves.size());
    leaves.insert(leaves.end(), _l->leaves.begin(), _l->leaves.end());
    leaves.insert(leaves.end(), _r->leaves.begin(), _r->leaves.end());
  }

  void Print(std::vector<int32_t> &rec) {
    if (height == 0) {
      rec.push_back(leaves[0]);
    }
    else{
      left->Print(rec);
      right->Print(rec);
    }
  }
};

inline void Delete(Node* rt) {
  if (rt->height == 0) {return;}
  if (rt->left != nullptr) {
    Delete(rt->left);
    delete rt->left;
  }
  if (rt->right !=nullptr) {
    Delete(rt->right);
    delete rt->right;
  }
};

class RareVariant {
public:
  bcf1_t* iv = NULL;
  int32_t ac, pos;
  int32_t **ibs0mat; // ac x ac matrix of ibs0 pos; upper-right; lower-left
  float *sorted_cm;  // ac x ac matrix in planer order, flatten by row as vector
  std::map<int32_t, uint32_t> id_index; // abs id -> index in ibs0mat
  std::vector<int32_t> id_list;         // Will be updated according to planer_order after building tree
  std::vector<int32_t> planer_order;
  int32_t left_size;
  int32_t AvgDist = 0;
  int32_t ovst = -1, oved = -1; // st & ed of the overlap region
  int32_t finished = 0;         // Count number of pairs with both ibs0 found
  float   AvgDist_cm = 0.0, MedDist_cm = 0.0;

RareVariant(bcf1_t* _iv, int32_t _ac) : iv(_iv), ac(_ac) {
  pos = iv->pos;
  sorted_cm = new float[ac*ac]{0.0};
  ibs0mat = new int32_t*[ac];
  for (int32_t i = 0; i < ac; ++i) {
    ibs0mat[i] = new int32_t[ac]{-1};
  }
}
~RareVariant() {
  for (int32_t i = 0; i < ac; ++i) {
    delete ibs0mat[i];
  }
  if (iv) {
    bcf_destroy(iv);
  }
}

void Add_id(int32_t id) { // Not robust!
  id_list.push_back(id);
  id_index[id] = id_list.size()-1;
}
void Add_id(std::vector<int32_t>& idvec) {
  int32_t ct = 0;
  id_list.clear();
  for (auto & v : idvec) {
    id_list.push_back(v);
    id_index[v] = ct;
    ct++;
  }
}

bool Add(int32_t id1, int32_t id2, int32_t st, int32_t ed);
bool Add_half(int32_t id1, int32_t id2, int32_t pt, int32_t direction);

bool IfDone() {
  return ( finished >= (ac*(ac-1)/2) );
}

// Build a tree
void Organize(bp2cmMap& pgmap);

};

#endif
