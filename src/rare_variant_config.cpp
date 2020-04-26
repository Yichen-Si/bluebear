#include "rare_variant_config.h"


// Build a tree. TODO: this is slow
void RareVariant::Organize(bp2cmMap& pgmap) {

  int32_t pt1=0, pt2=1, m=ac;
  std::vector<Node*> nodelist(m);
  for (int32_t it = 0; it < m; ++it) {
    nodelist[it] = new Node(it);
  }
  for (int32_t i = 0; i < m - 2; i++) { // m-2 merging events
    int32_t maxsim = 0;
    // Find closest branches to join
    for (int32_t p1 = 0; p1 < m-i-1; ++p1) {
      for (int32_t p2 = p1+1; p2 < m-i; ++p2) {
        int32_t avgsim = 0;
        int32_t ctsim = 0;
        for (auto & u : (nodelist[p1]->leaves) ) {
          for (auto & v : (nodelist[p2]->leaves) ) {
            int32_t irow = std::min(u,v), icol = std::max(u,v);
            avgsim += (ibs0mat[irow][icol] - ibs0mat[icol][irow]);
            ctsim += 1;
          }
        }
        avgsim /= ctsim;
        if (avgsim > maxsim) {
          maxsim = avgsim;
          pt1=p1; pt2=p2;
        }
      }
    }
    if (pt2 == m-i-1) {
      pt2 = pt1;
      pt1 = m-i-1;
    }
    else {
      Node* tmp = nodelist[m-i-1];
      nodelist[m-i-1] = nodelist[pt1];
      nodelist[pt1] = tmp;
    }
    // Merge two nodes
    Node* tmp = new Node(nodelist[pt2], nodelist[m-i-1], i+1);
    nodelist[pt2] = tmp;
  } // End of m-2 merging operations

  Node *root = new Node(nodelist[0], nodelist[1], m-1);
  root->Print(planer_order);
  left_size = root->left->leaves.size();

  // Re-order the ibd matrix according to the planer order
  for (int32_t i = 0; i < ac-1; ++i) {
    for (int32_t j = i+1; j < ac; ++j) {
      int32_t irow = std::min(planer_order[i], planer_order[j]);
      int32_t icol = std::max(planer_order[i], planer_order[j]);
      sorted_cm[i*ac+j] = pgmap.bpinterval2cm(pos, ibs0mat[irow][icol]);
      sorted_cm[j*ac+i] = pgmap.bpinterval2cm(ibs0mat[icol][irow], pos);
    }
  }

  std::vector<float> btwDist;
  for (auto & u : nodelist[0]->leaves ) {
    for (auto & v : nodelist[1]->leaves ) {
      int32_t irow = std::min(u,v), icol = std::max(u,v);
      AvgDist += (ibs0mat[irow][icol] - ibs0mat[icol][irow]);
      float d = pgmap.bpinterval2cm(ibs0mat[icol][irow], ibs0mat[irow][icol]);
      btwDist.push_back(d);
      AvgDist_cm += d;
    }
  }

  std::vector<int32_t> tmpid(id_list.begin(), id_list.end());
  for ( int32_t i = 0; i < ac; ++i ) {
    id_list[i] = tmpid[planer_order[i]];
  }

  uint32_t ctbtw = btwDist.size();
  std::sort(btwDist.begin(), btwDist.end());
  MedDist_cm = btwDist[ ctbtw/2 ];
  AvgDist /= ctbtw;
  AvgDist_cm /= ctbtw;

  Delete(root);

}

bool RareVariant::Add_half(int32_t id1, int32_t id2, int32_t pt, int32_t direction) {
  int32_t i = 0, j = 0, tmp = 0;
  try {
    i = id_index.at(id1);
    j = id_index.at(id2);
  }
  catch (const std::out_of_range& e) {
    error("RareVariant::Add_half ID out of range");
  }
  if (i > j) {
    tmp=i; i=j; j=tmp;
  }
  if (direction == 1) {  // Left
    ibs0mat[j][i] = pt;
    if (ovst < 0)
      ovst = pt;
    else
      ovst = std::max(ovst, pt);
    if (ibs0mat[i][j] != 0)
      finished ++;
  } else {              // Right
    ibs0mat[i][j] = pt;
    if (oved < 0)
      oved = pt;
    else
      oved = std::min(oved, pt);
    if (ibs0mat[j][i] != 0)
      finished ++;
  }

  if ( finished >= (ac*(ac-1)/2) ) {
    return 1;
  } else {
    return 0;
  }
}


bool RareVariant::Add(int32_t id1, int32_t id2, int32_t st, int32_t ed) {
  int32_t i=0, j=0, tmp=0;
  try {
    i = id_index.at(id1);
    j = id_index.at(id2);
  }
  catch (const std::out_of_range& e) {
    error("RareVariant::Add ID out of range");
  }
  if (i > j) {
    tmp=i; i=j; j=tmp;
  }
  ibs0mat[i][j] = ed;
  ibs0mat[j][i] = st;
  if (ovst < 0 || oved < 0) {
    ovst = st; oved = ed;
  } else {
    ovst = std::max(ovst, st);
    oved = std::min(oved, ed);
  }
  finished ++;
  if ( finished >= (ac*(ac-1)/2) ) {
    return 1;
  } else {
    return 0;
  }
}




