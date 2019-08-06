class pbwtCursor
{

protected:
  bool *y;
  int32_t *b, *e;

public:
  int32_t M;
  int32_t *a, *d;

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

  void SwitchIndex(std::set<int32_t>& idset) {
    int32_t ct = 2*idset.size();
    for (int32_t i = 0; i < M; ++i) {
      if (idset.find(a[i]/2) != idset.end()) {
        a[i] += 1 - 2 * (a[i]%2);
        ct--;
        if (ct ==0) break;
      }
    }
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



