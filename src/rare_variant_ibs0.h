#ifndef RARE_IBS0_H
#define RARE_IBS0_H

class RareVariant {
public:
  bcf1_t* iv;
  int32_t ac;
  int32_t **ibs0mat; // ac x ac matrix of ibs0 pos; upper-right; lower-left
  std::map<int32_t, uint32_t> id_index; // abs id -> index in ibs0mat
  std::vector<int32_t> id_list;
  std::vector<int32_t> subset1,subset2;
  int32_t AvgDist = 0;
  float AvgDist_cm = 0.0;
  int32_t ovst = -1, oved = -1; // st & ed of the overlap region
  int32_t finished = 0; // Count number of pairs with both ibs0 found

RareVariant(bcf1_t* _iv, int32_t _ac) : iv(_iv), ac(_ac) {
  ibs0mat = new int32_t*[ac];
  for (int32_t i = 0; i < ac; ++i) {
    ibs0mat[i] = new int32_t[ac]{0};
  }
}
~RareVariant() {
  for (int32_t i = 0; i < ac; ++i) {
    delete ibs0mat[i];
  }
  bcf_destroy(iv);
}

void Add_id(int32_t id) { // Not robust!
  id_list.push_back(id);
  id_index[id] = id_index.size();
}
void Add_id(std::vector<int32_t>& idvec) {
  int32_t ct = 0;
  for (auto & v : idvec) {
    id_list.push_back(v);
    id_index[v] = ct;
    ct++;
  }
}

bool Add(int32_t id1, int32_t id2, int32_t st, int32_t ed) {
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

bool Add_half(int32_t id1, int32_t id2, int32_t pt, int32_t direction) {
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

// Build a tree. TODO: this is slow
void Organize(bp2cmMap& pgmap) {
  int32_t pt1=0, pt2=1, m=ac;
  std::vector<std::vector<int32_t>* > nodelist(m);
  for (int32_t it = 0; it < m; ++it) {
    nodelist[it] = new std::vector<int32_t> {it};
  }
  for (int32_t i = 0; i < m - 2; i++) { // m-2 merging events
    int32_t maxsim = 0;
    // Find closest branches to join
    for (int32_t p1 = 0; p1 < m-i-1; ++p1) {
      for (int32_t p2 = p1+1; p2 < m-i; ++p2) {
        int32_t avgsim = 0;
        int32_t ctsim = 0;
        for (auto & u : (*nodelist[p1])) {
          for (auto & v : (*nodelist[p2])) {
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
      std::vector<int32_t>* tmp = nodelist[m-i-1];
      nodelist[m-i-1] = nodelist[pt1];
      nodelist[pt1] = tmp;
    }
    // Merge two nodes
    for (auto & v : *nodelist[m-i-1])
      (*nodelist[pt2]).push_back(v);
  } // End of m-2 merging operations

  // Ends up with two clusters.
  int32_t it = ((*nodelist[0]).size() < (*nodelist[1]).size()) ? 0 : 1;
  for (auto & v : *nodelist[it])
    subset1.push_back(id_list[v]); // Id in the smaller cluster
  for (auto & v : *nodelist[1-it])
    subset2.push_back(id_list[v]); // Id in the larger cluster

  int32_t ctbtw = 0;
  for (auto & u : *nodelist[0]) {
    for (auto & v : *nodelist[1]) {
      int32_t irow = std::min(u,v), icol = std::max(u,v);
      AvgDist += (ibs0mat[irow][icol] - ibs0mat[icol][irow]);
      AvgDist_cm += (pgmap.bp2cm(ibs0mat[irow][icol]) - pgmap.bp2cm(ibs0mat[icol][irow]));
      ctbtw += 1;
    }
  }
  AvgDist /= ctbtw;
  AvgDist_cm /= ctbtw;
  for (int32_t it = 0; it < m; ++it) {
    delete nodelist[it];
  }
}

};

#endif
