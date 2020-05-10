// analyze.cc - Extract data from data
// $Id$

#define FINITE_SIZE_SCALING

#define BINARY
#define ACCURACY 1e-12

#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "util.h"
#include "Table.h"

#include <algorithm>
#include <string>

#include "Average.h"

class observable : public Average {
public:
  ofstream file;

  observable(const char filnamn[]) : file(filnamn) {
    file.precision(8);
//    file << "@ legend on\n@ legend box fill on\n"
//      "@ legend x1 0.65\n@ legend y1 0.8\n";
  }

  void write(double T) {
    file << T _ *this << endl;
    reset();
  }

  void newset(int nr, string s) {
//    file << '&' << endl << "@ legend string " << nr << " \""
//	 << s << "\"" << endl;
  }
};

class Data {
public:

  double T, U;
  int Lx;
  Est curr;
  Est resistivity;
  Est voltage;
  int time;
  double step;
  int sweeps;

  Data() {}

#ifdef BINARY
  istream& read(istream& in) {
    in.Bread(&T,1);
    in.Bread(&U,1);
    in.Bread(&Lx,1);
    in.Bread(&curr,1);
    in.Bread(&resistivity,1);
    in.Bread(&voltage,1);
    in.Bread(&time,1);
    in.Bread(&step,1);
    in.Bread(&sweeps,1);
    return in; }
  //  istream& read(istream& in) { in.Bread(this,1); /*in.ignore(1);*/ return in; }
#else
  istream& read(istream& in) { 
    cerr << "No thanks!" << endl;
    return in;
  }
#endif

};

int compare(const void *a, const void *b) {
  const Data *x = (Data*) a, *y = (Data*) b;
  const double dT = x->T - y->T;
  const int    Lx = x->Lx, Ly = y->Lx;
  const double dU = x->U - y->U;
  const int    dt = x->time - y->time;
  if (Lx < Ly) return -1; else if (Lx > Ly) return +1;
  if (dT < -ACCURACY) return -1; else if (dT > +ACCURACY) return +1;
  if (dU < -ACCURACY) return -1; else if (dU > +ACCURACY) return +1;
  if (dt < 0) return -1; else if (dt > 0) return +1;
  return 0;
};

inline bool operator<(const Data& x, const Data& y)
{ return (compare(&x, &y) == -1); }

inline bool operator==(const Data& x, const Data& y)
{ return (compare(&x, &y) == 0); }

class Averages {
public:

  double T, U;
  int Lx;
  observable curr;
  observable resistivity;
  observable voltage;
  int time;
  double step;
  int sweeps;

  int n;			// Number of observations written to file

  Averages() : 
    curr("J"),
    resistivity("rho"),
    voltage("IV")
  { reset(); }

  void reset() {
    curr.reset();
    resistivity.reset();
    voltage.reset();

    n = 0;
    time = 0;
  }

  // Not used in disordered case:
  void operator +=(Data& d) {
    n++;

    T = d.T;
    U = d.U;

    Lx = d.Lx;
    
    time += d.time;
    step = d.step;
    sweeps = d.sweeps;

    curr += d.curr;
    resistivity += d.resistivity;
    voltage += d.voltage;
  }

  void write() {	// Not using replicas:
    int Vol = int(Lx);

    // curr	/= Vol*sweeps*step;
    curr	/= sweeps*step;
    resistivity	/= 2*sweeps*step*T/Vol;
    // voltage	is good as it is.

    // double xx = T;
    double xx = U;

#ifdef FINITE_SIZE_SCALING
    // ...

#endif

    double I = curr.mean();
    curr.write(xx);
    resistivity.write(xx);
    voltage.write(I); // IV

    cerr << n << " simulation runs"
	 << " at T = " << T << ", t = " << time << endl;
  }

  void newset(int nr) {
    char s[32];
    // sprintf(s,"%dx%d, t = %d",Lx,Ly, time);
    sprintf(s,"L = %d, T = %g",Lx, T);

    curr.newset(nr,s);
    resistivity.newset(nr,s);
    voltage.newset(nr,s);
  }
};

int main(int argc, char *argv[]) {

  istream& in = cin;

  cout.precision(8);

  int minsweep = 400;
  if (argc > 1)
    minsweep = atoi(argv[1]);

  int n = 0, nr = 0;

  const int Max = 40000;	// Max number of lines in raw

  Table<Data> d(Max);

  // First read everything:
  while (in.good() && d[n].read(in))
    n++;

  if (n > Max) {
    cerr << "Too much! Please increase Max!" << endl;
    exit(-2);
  }

  cerr << "Before averaging" << endl
       << "Number of lines: " << n << endl
       << "Number of bytes: " << n*sizeof(Data) << endl;

  // stable_sort(d.v,d.v+n);

  Averages O;

  for (int i = 0; i < n;) {

    O.reset();

    while (d[i].time < minsweep)
      i++;

    O += d[i++];
    while (i < n and d[i-1] == d[i])
      O += d[i++];

    O.write();

    if (d[i-1].Lx != d[i].Lx or d[i-1].T != d[i].T)
      O.newset(nr++);		// Mark a new set?
  }

  //  system("xmgr -nxy resist");
}
