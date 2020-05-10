// plot-v.cc - Extract voltage differences from V.
// Usage: plot-v < V > VV; xmgr -nxy VV

#define BINARY

#include <math.h>
#include <fstream>
#include "util.h"

int main(int argc, char *argv[]) {

  istream& in = cin;

  cout.precision(8);

  int t = 0;
  double v[4];
  while (in.good() && in.Bread(v,4)) {
    cout << v[0] __ v[1] __ v[2] __ v[3] << endl;
  }
}
