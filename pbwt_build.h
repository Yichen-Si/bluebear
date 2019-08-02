class pbwtCursor
{

private:
  bool *y;
  int32_t *b, *e;

public:
  int32_t M;
  int32_t *a, *d;

  pbwtCursor(int32_t _M, int32_t k = INT_MAX) {
    M = _M;
    a = new int32_t[M];
    b = new int32_t[M];
    d = new int32_t[M];
    e = new int32_t[M];
    for (int32_t i = 0; i < M; ++i) {
      d[i] = k; a[i] = i;
    }
  }
  ~pbwtCursor() {
    delete [] a;
    delete [] b;
    delete [] d;
    delete [] e;
    d = NULL;
  }

  void ForwardsAD(bool *y, int32_t k) {
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

  void ReverseA(int32_t *r) {
    for (int32_t i = 0; i < M; ++i) {
      r[a[i]] = i;
    }
  }

};
