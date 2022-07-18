// Estimate.h  -*-C++-*- class for error estimation.
//	$Id: Estimate.h,v 1.2 2009/08/31 09:20:35 jack Exp $	

// Estimate errors using jacknife method.

#ifndef Estimate_H
#define Estimate_H

#define MEMORY 64

class Estimate {
public:
  const int M;
  int n, i, j;

  double *mi, m, v, m2;

  bool jack;

  Estimate(int memory = MEMORY) : M(2*(memory/2)) {
    mi = new double[M];
    reset();
  }

  virtual ~Estimate() {
    delete [] mi;
  }

  void reset() {
    n = 1; i = j = 0;
    m = v = m2 = 0;
    for (int k = 0; k < M; mi[k++]=0);
    jack = false;
  }

  void operator+=(const double x)  {
    if (j == M) rearr();
    mi[j] += x;
    m2    += x*x;
    i++;
    if (i == n) { i = 0; j++; }
  }

  void rearr() {
    j = 0;
    for (int k = 0; k < M; j++, k++, k++)
      mi[j] = mi[k] + mi[k+1];
    for (int k = j; k < M; k++)
      mi[k] = 0;
    n *= 2;
  }

  void jackknife() {		// Form Jackknife bins.
    if (jack) return;
    int N = n*j;		// Skipping the last i samples.
//     if (i != 0)
//       cerr << "jackknife: i=" << i // For optimal performance i should be 0.
// 	   << " j=" << j
// 	   << " n=" << n << endl;
    if (N == 0) return;
    m = v = 0;
    for (int k = 0; k < j; k++) {
      m += mi[k];
    }
    double mean = m/N;
    for (int k = 0; k < j; k++) {
      mi[k] = m - mi[k];	// Jackknife bins.
      v += sqr(mi[k]/(N-n)-mean);
    }
    m = mean;			// Average
    v *= j-1;
    v /= j; // already done...
    v = sqrt(v);		// Standard error
    // v = sqrt(v/(j+1));	// Standard error
    m2 = m2/(N+i) - m*m;	// Variance

    jack = true;
  }

  Estimate& operator /=(const Estimate& w) {
    if (!w.jack) {
      cerr << "Please apply jackknife to argument first!" << endl;
      abort();
    }
    jackknife();		// First jackknife

    double m_J = 0;
    for (int k = 0; k < j; k++) {
      mi[k] /= w.mi[k];
      m_J += mi[k];
    }
    m_J /= j;			// Jackknife average.

    m /= w.m;			// Ordinary average.

    m += (j-1)*(m - m_J);	// Correct for bias.

    double v_J = 0;
    for (int k = 0; k < j; k++)
      v_J += sqr( mi[k] - m_J );

    v_J *= (j-1.0)/j;		// Variance of m_J.

    v_J = sqrt(v_J);		// Standard error.
    v = v_J;			// Standard error.

    // m2 = ?

    return *this;
  }

  // Usefull for calculating cummulants, eg., cv, chi.
  Estimate& subsqr(const Estimate& A) {
    if (!A.jack) {
      cerr << "Please apply jackknife to argument first!" << endl;
      abort();
    }
    jackknife();		// First jackknife

    double m_J = 0;
    for (int k = 0; k < j; k++) {
      mi[k] -= sqr(A.mi[k]);
      m_J += mi[k];
    }
    m_J /= j;			// Jackknife average.

    m -= sqr(A.m);		// Ordinary average.

    m += (j-1)*(m - m_J);	// Correct for bias.

    double v_J = 0;
    for (int k = 0; k < j; k++)
      v_J += sqr( mi[k] - m_J );

    v_J *= (j-1.0)/j;		// Variance of m_J.

    v_J = sqrt(v_J);		// Standard error.
    v = v_J;			// Standard error.

    // m2 = ?

    return *this;
  }

  void results() {
    int N = n*j + i;
    if (N == 0) return;
    m = v = 0;
    for (int k = 0; k < j; k++) {
      m += mi[k];
    }
    if (i > 0) m += mi[j];
    m /= N;			// Average
    for (int k = 0; k < j; k++) {
      v += sqr(mi[k]/n-m);
    }
    if (i > 0)
      v += (sqr(mi[j]/i-m)*i)/n; // if the last bin is not yet full.

    v /= j-1;
    v = sqrt(v/j);		// Standard error
    m2 = m2/N - m*m;		// Variance
  }

  void bootstrap(int nboot = 200) {
    int N = n*j + i;
    if (N == 0) return;
    m = v = 0;
    for (int k = 0; k <= j; k++) // Use ordinary average
      m += mi[k];
    m /= N;			// Average
    int kmax = j - (i!=0);
    for (int t = 1; t <= nboot; t++) {
      double mm = 0;
      for (int k = 0; k < kmax; k++) {
	int l = int(rnd()*kmax);
	mm += mi[l];
      }
      mm /= n*kmax;
      mm -= m;
      v += mm*mm;
    }
    v /= nboot;
    v = sqrt(v);		// Standard error
    m2 = m2/N - m*m;		// Variance
  }

  void bootstrap2(int nboot = 200) {
    int N = n*j + i;
    if (N == 0) return;
    m = v = 0;
    int kmax = j - (i!=0);
    for (int t = 1; t <= nboot; t++) {
      double mm = 0;
      for (int k = 0; k < kmax; k++) {
	int l = int(rnd()*kmax);
	mm += mi[l];
      }
      mm /= n*kmax;

      m += mm;
      v += mm*mm;
    }
    m /= nboot;
    v /= nboot;
    v -= m*m;


//       mm -= m;

//       m += mm/t;
//       v += (t-1)*mm*mm/t;
//     }
//     v /= nboot;
    v = sqrt(v);		// Standard error
    m2 = m2/N - m*m;		// Variance
  }

  friend ostream& operator << (ostream& out, const Estimate& o);

  friend obstream& operator << (obstream& out, const Estimate& o);

  void write(ostream& out) {
    for (int k = 0; k < j; k++)
      out << k*n _ mi[k]/n << endl;
    out   << j*n _ mi[j]/i << endl << '&' << endl;
  }

  void save(ostream& out) {
    out << j _ n _ i _ m _ v _ m2 << endl;
    out.write((char*) mi,(j+1)*sizeof(*mi));
  }

  void load(istream& in) {
    in >> j >> n >> i >> m >> v >> m2 >> newln;
    in.read((char*) mi,(j+1)*sizeof(*mi));
  }

  // double time() { return (m2 > 0 ? v*v/m2*(n*j+1) : 0); }
  double time() { return (m2 > 0 ? v*v/m2*(n*j+i) : 0); }

  // Integrated autocorrelation function:
  double icorr() { return v*v*(n*j+i); }
};

ostream& operator << (ostream& out, const Estimate& o)
{ return out << o.m _ o.v; }

obstream& operator << (obstream& out, const Estimate& o)
{ return out << o.m << o.v; }

#undef MEMORY

#endif
