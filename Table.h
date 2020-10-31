// Simple -*- C++ -*- tables
//	$Id$

#ifndef Table_H
#define Table_H

#define Max_Dim 6

#ifdef CHECK_OUT_OF_BOUNDS
# define ch_b(x,i)  if (x < 0 || x >= L[i]) out_of_bounds(x,L[i]);
# define ch_bn(x,i) if (x < 0 || x >= i )   out_of_bounds(x,i);
# define ch_dim(i)  if (i != d)             out_of_bounds_dim(i);
// extern "C" void exit(int);
#else
# define ch_b(x,i)
# define ch_bn(x,i)
# define ch_dim(i)
#endif

template <class T>
class Table
{
public:
  T *v;
  int N;
  int d;
  int L[Max_Dim];
  int own_;

  Table() : own_(0), N(0), d(0) {}
  Table(int x): own_(0) { init(x); }
  Table(int x, int y): own_(0) { init(x,y); }
  Table(int x, int y, int z): own_(0) { init(x,y,z); }
  Table(int i, int j, int k, int l): own_(0) { init(i,j,k,l); }
  Table(int i, int j, int k, int l, int m): own_(0) { init(i,j,k,l,m); }
  Table(int i, int j, int k, int l, int m, int n):own_(0){ init(i,j,k,l,m,n); }

  ~Table() { clear(); }

//   Table(const Table<T>& table) { inherit((Table<T>& ) table); }
//   Table(Table<T>& table) { inherit(table); }

  T* begin() { return v; }
  T* end()   { return v+N; }

  void init(int x) {
    if (own_) clear();
    L[0] = N = x; d=1; own_=1;
    v = new T[x];
  }

  void init(int x, int y) {
    if (own_) clear();
    L[0] = x; L[1] = y; d=2; own_=1;
    N = x*y; v = new T[N];
  }

  void init(int x, int y, int z) {
    if (own_) clear();
    L[0] = x; L[1] = y; L[2] = z; d=3; own_=1;
    N = x*y*z; v = new T[N];
  }

  void init(int i, int j, int k, int l) {
    if (own_) clear();
    L[0] = i; L[1] = j; L[2] = k, L[3] = l; d=4; own_=1;
    N = i*j*k*l; v = new T[N];
  }

  void init(int i, int j, int k, int l, int m) {
    if (own_) clear();
    L[0] = i; L[1] = j; L[2] = k, L[3] = l, L[4] = m; d=5; own_=1;
    N = i*j*k*l*m; v = new T[N];
  }

  void init(int i, int j, int k, int l, int m, int n) {
    if (own_) clear();
    L[0] = i; L[1] = j; L[2] = k, L[3] = l, L[4] = m, L[5] = n; d=6; own_=1;
    N = i*j*k*l*m*n; v = new T[N];
  }

  void clear() {
    if (own_) { delete [] v; own_=0; }
  }

  void ref(Table<T>& o) {
    if (own_) clear();
    v = o.v;
    N = o.N;
    d = o.d;
    L[0] = o.L[0]; L[1] = o.L[1]; L[2] = o.L[2];
    L[3] = o.L[3]; L[4] = o.L[4]; L[5] = o.L[5]; 
  }

  void copy(Table<T>& o) { // Copy contents from o.
    if (own_) clear();
    N = o.N;
    d = o.d;
    L[0] = o.L[0]; L[1] = o.L[1]; L[2] = o.L[2];
    L[3] = o.L[3]; L[4] = o.L[4]; L[5] = o.L[5]; 
    own_ = 1;
    v = new T[N];
    for (int i = 0; i < N; i++) v[i] = o.v[i];
  }

  void inherit(Table<T>& o) { // Transfer ownership from o.
    ref(o);
    own_ = o.own_; o.own_ = 0;
  }

  T& operator()(int x) { ch_dim(1) ch_b(x,0) return v[x]; }
  T& operator()(int x, int y) 
  { ch_dim(2) ch_b(x,0) ch_b(y,1) return v[x + L[0]*y]; }
  T& operator()(int x, int y, int z)
  { ch_dim(3) ch_b(x,0) ch_b(y,1) ch_b(z,2) return v[x + L[0]*(y + L[1]*z)]; }
  T& operator()(int i, int j, int k, int l)
  { ch_dim(4) ch_b(i,0) ch_b(j,1) ch_b(k,2) ch_b(l,3)
      return v[i + L[0]*(j + L[1]*(k + L[2]*l))]; }
  T& operator()(int i, int j, int k, int l, int m)
  { ch_dim(5) ch_b(i,0) ch_b(j,1) ch_b(k,2) ch_b(l,3) ch_b(m,4)
      return v[i + L[0]*(j + L[1]*(k + L[2]*(l + L[3]*m)))]; }
  T& operator()(int i, int j, int k, int l, int m, int n)
  { ch_dim(6) ch_b(i,0) ch_b(j,1) ch_b(k,2) ch_b(l,3) ch_b(m,4) ch_b(n,5)
      return v[i + L[0]*(j + L[1]*(k + L[2]*(l + L[3]*(m + L[4]*n))))]; }

  T& operator[](int m) { ch_bn(m,N) return v[m]; }

  // operator T*() { return v; } These are dangerous as they can appear as L-values.
  // operator void*() { return v; }

  void set(const T x) { for (int i = 0; i < N; v[i++] = x) {} }
  void set(const T x, int n) { for (int i = 0; i < n; v[i++] = x) {} }
  T operator=(const T x) { for (int i = 0; i < N; v[i++] = x) {} return x; }

  int index(int x) { return x; }
  int index(int x, int y) { return x + L[0]*y; }
  int index(int x, int y, int z) { return x + L[0]*(y + L[1]*z); }
  int index(int i, int j, int k, int l)
    { return i + L[0]*(j + L[1]*(k + L[2]*l)); }
  int index(int i, int j, int k, int l, int m)
    { return i + L[0]*(j + L[1]*(k + L[2]*(l + L[3]*m))); }
  int index(int i, int j, int k, int l, int m, int n)
    { return i + L[0]*(j + L[1]*(k + L[2]*(l + L[3]*(m + L[4]*n)))); }

  int dim() { return d; }

  void save(ostream & fil) {
    fil << N _ d << endl;
    // fil << L[0] _ L[1] _ L[2] _ L[3] _ L[4] _ L[5] << endl;
    for (int i = 0; i < d; i++) fil _ L[i];
    fil << endl;
    fil.write((char*) v, N*sizeof(T));
  }
  void load(istream & fil) {
    fil >> N >> d >> newln;
    // fil >> L[0] >> L[1] >> L[2] >> L[3] >> L[4] >> L[5] >> newln;
    for (int i = 0; i < d; i++) fil >> L[i];
    fil >> newln;
    if (own_) clear(); own_=1;
    v = new T[N];
    fil.read((char*) v, N*sizeof(T));
  }

  void write(ostream & fil) {	// Print contents in a column
    for (int i = 0; i < N; i++)
      fil << v[i] << endl;
  }

#ifdef CHECK_OUT_OF_BOUNDS
  void out_of_bounds(int x, int i) {
    cerr << "Out of bounds error: "
	 << x << " not in [ 0 , " << i-1 << " ]" << endl;
#if CHECK_OUT_OF_BOUNDS+100 == 2+100
    abort();
#else
    exit(1);
#endif
  }

  void out_of_bounds_dim(int i) {
    cerr << "Dimension mismatch: " << i << " != " << d << endl;
#if CHECK_OUT_OF_BOUNDS+100 == 2+100
    abort();
#else
    exit(1);
#endif
  }
#endif

};

template <class T>
int search(T v, T *table, int n) {
  int i, l = 0;	     // Binary search of ordered table.
  while (n > l) {
    i = (l+n)/2;
    if (v < table[i]) n = i - 1; else l = i;
  }
  return l;
}

# undef ch_b
# undef ch_bn
# undef ch_dim

#endif
