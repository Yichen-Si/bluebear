#ifndef PBWT_BUILD_H
#define PBWT_BUILD_H

class pbwtCursor
{

protected:
  bool *y;
  int32_t* b; // temporary storage for a;
  int32_t* e; // temporary storage for d;

public:
  int32_t M;  // Number of haplotypes
  int32_t* a; // a[i]: original index [M] of the haplotype for i-th prefix
  int32_t* d; // Start of prefix matching between a[i] and a[i-1]

  /**
   * For prefix, k is the first (leftmost) position
   * For suffix, k is the last (rightmost) position
   */
  pbwtCursor(int32_t _M, int32_t k) {
    M = _M;
    a = new int32_t[M];
    b = new int32_t[M];
    d = new int32_t[M];
    e = new int32_t[M];
    for (int32_t i = 0; i < M; ++i) {
      d[i] = k; a[i] = i;
    }
  }

  // Initialize from snapshot
  pbwtCursor(int32_t _M, int32_t* _a, int32_t* _d) {
    M = _M;
    a = new int32_t[M];
    b = new int32_t[M];
    d = new int32_t[M];
    e = new int32_t[M];
    memcpy (a, _a, M*sizeof(int32_t));
    memcpy (d, _d, M*sizeof(int32_t));
  }

  ~pbwtCursor() {
    delete [] a;
    delete [] b;
    delete [] d;
    delete [] e;
    d = NULL; a = NULL;
  }

  void Reset(int32_t _M, int32_t k) {
    M = _M;
    for (int32_t i = 0; i < M; ++i) {
      d[i] = k; a[i] = i;
      b[i] = 0; e[i] = 0;
    }
  }

  void Reset(int32_t _M, int32_t* _a, int32_t* _d) {
    M = _M;
    memcpy (a, _a, M*sizeof(int32_t));
    memcpy (d, _d, M*sizeof(int32_t));
    for (int32_t i = 0; i < M; ++i) {
      b[i] = 0; e[i] = 0;
    }
  }

  /**
   * Should be called backward - from larger positions to smaller positions
   */
  void ForwardsAD_suffix(bool *y, int32_t k) {
    int32_t u = 0, v = 0;
    int32_t p = k-1, q = k-1;
    for (int32_t i = 0; i < M; ++i) {
      if (d[i] < p) p = d[i];
      if (d[i] < q) q = d[i];
      if (y[a[i]] == 0) {
        a[u] = a[i];
        d[u] = p;
        ++u; p = INT_MAX;
      }
      else {
        b[v] = a[i];
        e[v] = q;
        ++v; q = INT_MAX;
      }
    }
    memcpy(a+u, b, v*sizeof(int32_t));
    memcpy(d+u, e, v*sizeof(int32_t));
  }

  void ForwardsAD_prefix(bool *y, int32_t k) {
    int32_t u = 0, v = 0;
    int32_t p = k+1, q = k+1;
    for (int32_t i = 0; i < M; ++i) {
      if (d[i] > p) p = d[i];
      if (d[i] > q) q = d[i];
      if (y[a[i]] == 0) {
        a[u] = a[i];
        d[u] = p;
        ++u; p = 0;
      }
      else {
        b[v] = a[i];
        e[v] = q;
        ++v; q = 0;
      }
    }
    memcpy(a+u, b, v*sizeof(int32_t));
    memcpy(d+u, e, v*sizeof(int32_t));
  }

  void ReverseA(int32_t *r) {
    for (int32_t i = 0; i < M; ++i) {
      r[a[i]] = i;
    }
  }

  void SwitchIndex(std::vector<bool>& indi) {
    for (int32_t i = 0; i < M; ++i) {
      if (indi[a[i]/2]) {
        a[i] += 1 - 2 * (a[i]%2);
      }
    }
  }

  /**
   * Start position (inclusive) of the prefix matching between i and j
   */
  int32_t Dist_pref(int32_t i, int32_t j) {
    if (i > j) {
      int32_t tmp = j; j = i; i = tmp;
    }
    int32_t st = d[i+1];
    for (int32_t it = i+1; it <= j; ++it) {
      if (d[it] > st) {
        st = d[it];
      }
    }
    return st;
  }

  /**
   *
  */
  int32_t Dist_suff(int32_t i, int32_t j) {
    if (i > j) {
      int32_t tmp = j; j = i; i = tmp;
    }
    int32_t ed = d[i+1];
    for (int32_t it = i+1; it <= j; ++it) {
      if (d[it] < ed) {
        ed = d[it];
      }
    }
    return ed;
  }

  void SwitchIndex_OnePair(int32_t x, int32_t y) {
    for (int32_t i = 0; i < M; ++i) {
      if (a[i] == x) {
        a[i] = y;
      } else if (a[i] == y) {
        a[i] = x;
      }
    }
  }

  int32_t FindIndex(int32_t x) {
    for (int32_t i = 0; i < M; ++i) {
      if (a[i] == x)
        return i;
    }
    return -1;
  }

};

#endif
