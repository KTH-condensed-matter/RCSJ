// Average.h -  -*-C++-*- class for averages and errorbars.
//	$Id$

// Assume samples are independent and errorbars reliable.

#ifndef Average_H
#define Average_H

class Est {
public:
  double m;			// value
  double err;			// standard error
};

class Average {
public:

  double m;			// the average
  double w;			// weight

  Average() {
    reset();
  }

  virtual ~Average() {}

  void reset() {
    m = w = 0;
  }

  void operator += (const Average a) {
    w += a.w;
    double d = a.m - m;
    d *= a.w;
    d /= w;
    m += d;
  }

  void samp(double am, double aw) {
    w += aw;
    double d = am - m;
    d *= aw;
    d /= w;
    m += d;
  }

  void operator += (const Est e) {
    double ew = 1/sqr(e.err);
    if (e.err == 0) ew = 100;
    w += ew;
    double d = e.m - m;
    d *= ew;
    d /= w;
    m += d;
  }

  // Assumes each sample has equal weight 1:
  void operator += (const double x) {
    w ++;
    double d = x - m;
    d /= w;
    m += d;
  }

  double& mean() { return m; }
  double err()  { return 1.0/sqrt(w); }

  void results() {}

  // Allow simple transformations:
  void operator*=(const double z) { m *= z; w /= z*z; }
  void operator/=(const double z) { m /= z; w *= z*z; }

  friend ostream& operator << (ostream& out, const Average& o);

  friend obstream& operator << (obstream& out, const Average& o);

  void save(ostream& out) {
    out.write((char*) this, sizeof(*this));
  }

  void load(istream& in) {
    in.read((char*) this,sizeof(*this));
  }

};

ostream& operator << (ostream& out, const Average& a)
{ return out << a.m _ 1.0/sqrt(a.w); }

obstream& operator << (obstream& out, const Average& a)
{ return out << a.m << 1.0/sqrt(a.w); }

#endif
