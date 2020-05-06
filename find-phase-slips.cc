// find-phase-slips.cc - Extract voltage differences from Vx.
// $Id$

#define BINARY

#include <math.h>
#include <fstream>
#include "util.h"
#include "Table.h"
#include <stdlib.h>

int main(int argc, char *argv[]) {

  int Lx = 0, skip = 0;

  if (argc > 1) {
    Lx = atoi(argv[1]);
    if (argc > 2)
      skip = atoi(argv[2]);
  }
  else {
    cerr << "Usage: " << argv[0] << " Lx [skip] < Vx_...\n"
    exit(-1);
  }

  istream& in = cin;

  cout.precision(8);

  Table<double> V(Lx);

  int t = 0;
  while (in.good() && in.Bread(V.v,Lx)) {
    if (t++ < skip) continue;
    for (int x = 0; x < Lx-1; x++)
      cout _ V(x+1)-V(x);
    cout << endl;
  }

  //  system("xmgr -nxy resist");
}
