// rcsj.cc - RCSJ simulation of a 1D XY model with shunting capacitors.
// 1D version.

// $Id: rcsj.cc,v 1.32 2013/06/25 12:15:46 jack Exp $

// #define CHECK_OUT_OF_BOUNDS 2

// #define RANDOM_INITIAL_CONDITION
//  const double Ccomb = 10*C; // Comb capacitor 10 times junction capacitance. It is also possible to define multiple of Co
// x % 10 == 5 means statring from junction number 5 with periodicity of 10
#define FIND_PHASE_SLIPS
// #define PRINT_VOLTAGE

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include <cstdlib>
#include <cstdio>
#include <signal.h>

#include "mt19937.h"
#include "util.h"
#include "Table.h"
#include "Estimate.h"

static char rcsid[] =
"$Id: rcsj.cc,v 1.32 2013/06/25 12:15:46 jack Exp $";

static char usage[] =

  "Possible actions:\n"
  " #	...			-- Comment.\n"
  " !	...			-- Execute a shell command.\n"
  " H				-- Help.\n"
  " + lots of undocumented stuff\n"
  " Q				-- quit.\n";

class System* sys;
int please_leave = 0;

#ifndef NO_GRAPHICS
#include "graph.h"
VortexWin *win = 0;
#endif

int filenumber = 0;

#include "TriSolve.h"

typedef double angle;		// For now use a double to represent
				// the angle.

// Some variables defining the circuit and its properties:
enum layout{ordinary, midcap, comb, weak_link, all_different} circuit_layout = ordinary;

bool init_voltages = true;

layout read_layout(istream& in) {
  string s;
  in >> s;
  if (s == "ordinary") return ordinary;
  if (s == "midcap") return midcap;
  if (s == "comb") return comb;
  if (s == "weak_link") return weak_link;
  if (s == "all_different") return all_different;
  cerr << "What do you mean?" << endl;
  return ordinary;
}

class Tri : public TriSolve
{
public:
  Table<double> L, D, U;

  void init(int Lx, double C, double C0, double R, double Rterm, double Cterm, double dt, double Cm = 0) {
    N = Lx;
    L.init(N-1);
    D.init(N);
    U.init(N-1);

    // Fill in the capacitance and conductance matrices.
    // double M = C + dt/2/R;
    double M = C;		// XXX No R here! Euler.
    for (int x = 0; x < N; x++) {
      if (circuit_layout == comb && Cm > 0 && x % 10 == 5)
	// x % 10 == 5 means statring from junction number 5 with periodicity of 289  XXXXXXXXXXXXXXXXXXXXXXXXXXX
	D[x] = Cm;		// Capacitor comb every tenth.
      else
	D[x] = C0;		// ordinary circuit.
      if (x+1 < N) { // Right neighbour:
	D[x] += M;
	U[x] = - M;
      }
      if (x > 0) { // Left neighbour:
	D[x] += M;
	L[x-1] = - M;
      }
    }
    D[0]   += dt/2/Rterm; // First junction connected via Rterm to voltage source.
    D[N-1] += dt/2/Rterm; // Last junction terminated by Rterm to ground.

    // Terminal capacitance:
    if (Cterm > 0) { D[0] += Cterm - C0; D[N-1] += Cterm - C0; }

    // Middle capacitor if nonzero:
    if (circuit_layout == midcap && Cm > 0) D[N/2] += Cm - C0; // Replace C0 by Cm.

    // // Weak link at middle:
    // double Mw = Cw + dt/2/Rw;
    // D[N/2]   += Mw - M;
    // U[N/2]   -= Mw - M;
    // D[N/2-1] += Mw - M;
    // L[N/2-1] -= Mw - M;

    // L.write(cerr); cerr << '&' << endl;
    // D.write(cerr); cerr << '&' << endl;
    // U.write(cerr); cerr << '&' << endl;

    lu(L.v,D.v,U.v,N);

    // L.write(cerr); cerr << '&' << endl;
    // D.write(cerr); cerr << '&' << endl;
    // U.write(cerr); cerr << '&' << endl;
  }

};

class System
{
public:

  int Lx;			// System size.

  Table<double> V;		// Voltages;
  Table<angle> theta;		// Phases theta_r.
  Table<angle> new_theta;
  Table<angle> old_theta;

  Table<double> Is;		// Supercurrents (temporary)

  Table<double> Ic_array;	// Critical currents array of length Lx-1.
  Table<double> C0_array;	// Capaciatance to ground array of length Lx.
  Table<double> C_array;	// Capacitance of junctions array of length Lx-1.

  double T;			// Temperature.
  double J;			// Electric current (measured).
  double Ic, Icweak;		// Critical current of tunnel junctions.
  double U;			// Applied voltage.
  double C, C0, Cterm;		// Junction capacitance, and capacitance to ground.
  double R, Rterm;		// Junction resistance and terminal resistance.
  double Vgap;			// Gap voltage (or 0).
  double Vac, omega;		// Applied frequence dependent voltage.

  int Normaljunctionnumber;	// Junction which goes normal when a photon hits.
  double gapSupressionTime;	// Time scale for gap supression.
  double gapRecoveryTime;	// Time scale for gap recovery.
  double photonArrivalFreq;	// Frequency of photon arrivals.
  double photonArriavaltime;	// Time of next photon arrival.
  
  double Energy;		// Current energy (not used).

  int time;			// Time (number of md steps).
  int number, equ;		// Number of sweeps and sweeps to equilibrate.
  int sweeps;			// Time between measurements.
  double step;			// Step size for theta updates.

  Tri eqsys;

  ofstream rawfile;		// append everything to raw.
  obstream datafile;		// Binary data file.

  ofstream phout;		// Data file containing phase slip events (t, x).

  System() {}
  virtual ~System() {}

  bool setparam(istream&);
  void run_T_sweep(istream& cmd);
  void run_U_sweep(istream& cmd);
  int main(int argc, char *argv[]);

  void init(int Lx_) {

    Lx = Lx_;

    V.init(Lx);
    V = 0.0;
    theta.init(Lx);
    theta = 0.0;

#ifdef RANDOM_INITIAL_CONDITION
    for (int i = 0; i < theta.N; i++)
      theta[i] = twopi*rnd();
#endif

    new_theta.init(Lx);
    old_theta.init(Lx);
    Is.init(Lx+1);		// Currents, including leads.

    Ic_array.init(Lx-1);	// Critical currents.
    set_Ic();
    Energy = 0;

    init_voltages();
  }

  void set_defaults() {		// Default paramter values:
    step = 0.02;		// Time step in MD case.
    C = 0.1;
    C0 = C/100;
    R = 1;
    Rterm = R/100;		// 50 Ohms instead?
    Cterm = 0;			// 0 means same as C0.
    Vgap = 0;
    Ic = 1.0;
    Icweak = 0;			// 0 means 0.1*Ic.
    U = 0.0;
    Vac = 0.0;			// No AC.
    omega = 0.0;		// -"-
    // Ic_array = Ic;		// Same critical current by default.

    Normaljunctionnumber = Lx/2;
    gapSupressionTime = 0.0;
    gapRecoveryTime = 4.0; // XXX Units?
    photonArrivalFreq = 0;
  }

  void set_Ic() {
    if (Ic_array.v == 0) return;
    Ic_array = Ic;
    if (circuit_layout == weak_link) {
      double Icw = Icweak > 0 ? Icweak : Ic*0.1;
      Ic_array(Lx/2-1) = Icw;
    }
  }

  int read_Ic(istream& cmd) {
    string fname;
    cmd >> fname;
    if (fname == "random") {
      double stdev;
      cmd >> stdev >> fname;
      // setup random Ic array.
      for (int x = 0; x < Lx-1; x++) {
	Ic_array(x) = Ic + stdev*rnd.normal();
	if (Ic_array(x) <= 0) cerr << "Warning: Ic <= 0." << endl;
      }
      // save random Ic array.
      ofstream out(fname.data());
      if (!out) { cerr << "Cannot open " << fname.data() << "!" << endl; exit(1); }
      for (int x = 0; x < Lx-1; x++) {
	out << Ic_array[x] << endl;
      }
    }
    else {
      // read Ic array.
      ifstream in(fname.data());
      if (!in) { cerr << "Cannot open " << fname.data() << "!" << endl; exit(1); }
      for (int x = 0; x < Lx-1; x++) {
	in >> Ic_array[x];
      }
    }
    return !!cmd;
  }

  int read_C0(istream& cmd) {
    string fname;
    cmd >> fname;
    ifstream in(fname.data());
    if (!in) { cerr << "Cannot open " << fname.data() << "!" << endl; exit(1); }
    C0_array.init(Lx);
    for (int x = 0; x < Lx; x++)
      in >> C0_array[x];
    return !!cmd;
  }

  int read_C(istream& cmd) {
    string fname;
    cmd >> fname;
    ifstream in(fname.data());
    if (!in) { cerr << "Cannot open " << fname.data() << "!" << endl; exit(1); }
    C_array.init(Lx-1);
    for (int x = 0; x < Lx-1; x++)
      in >> C_array[x];
    return !!cmd;
  }

  int read_photon_param(istream& cmd) {
    cmd >> Normaljunctionnumber
	>> gapSupressionTime
	>> gapRecoveryTime
	>> photonArrivalFreq;
    photonArriavaltime = -100000.0;
    return !!cmd;
  }

  void init_voltages() {	// Initialize voltages to a linear profile.
    // const double Ic = 1.0;
    if (U < 2*Rterm*Ic) {
      for (int x = 0; x < Lx; x++)
	V(x) = U/2;
    }
    else {
      double Rtot = Rterm + (Lx-1)*R + Rterm;
      double Iarr = (U-2*Rterm*Ic)/Rtot;
      double Itot = Iarr + Ic;
      for (int x = 0; x < Lx; x++)
	V(x) = U - Itot*Rterm - Iarr*R*x;
    }
  }

  void setup(double T_, double U_) {
    T = T_;
    U = U_;

    const double Ccomb = 5*C; // Comb capacitor 10 times junction capacitance. It is also possible to define multiple of Co XXXXX
    
    // Setup capacitance matrix and LU factorize.
    // eqsys.init(Lx, C, C0, R, Rterm, step, 100*C); // Middle capacitance Cm = C >> C0.
    if (circuit_layout == comb || circuit_layout == midcap)
      eqsys.init(Lx, C, C0, R, Rterm, Cterm, step, Ccomb);
    else
      eqsys.init(Lx, C, C0, R, Rterm, Cterm, step);

    if (::init_voltages)
      init_voltages();
    System::reset();

  }

  virtual void reset() {
    J = 0;
    time = 0;
  }

  inline int xplus(const int x)  { return (x+1 == Lx) ? 0 : x+1; }
  inline int xminus(const int x) { return (x == 0) ? Lx-1 : x-1; }

  // Lx = number of sc islands. number of JJ = Lx-1

  void rcsj(int ant = 1) {

    for (int iii = 0; iii < ant; iii++) {

      if (photonArrivalFreq > 0 &&
	  photonArriavaltime + 1.0/photonArrivalFreq < time*step)
	photonArriavaltime = time*step + 1.0/photonArrivalFreq;

      // Voltage bias: Delta is twisted with time!

      // Calculate all supercurrents + noise, and the Delta update.

      double Ut = U + Vac*sin(omega*time*step); // Applied voltage, possibly time dependent.

      Is(0) = (Ut - V[0])/Rterm + sqrt(2*T/Rterm/step)*rnd.normal(); // Current through left lead.
      for (int x = 0; x < Lx-1; x++) { // Step through the Lx-1 junctions.
	double IR = 0, In = 0;
	if (abs(V(x) - V(x+1)) >= Vgap) {
	  IR = (V(x) - V(x+1))/R; // Current through resistor shunt.
	  In = sqrt(2*T/R/step)*rnd.normal();
	}
	Is(x+1) = Ic_array(x)*sin( theta(x) - theta(x+1) ) + IR + In;


	if (x == Normaljunctionnumber) {
	  double vgap = Vgap;
	  if (photonArriavaltime < time*step && time*step < photonArriavaltime + gapSupressionTime) {
	    double delta = time*step - photonArriavaltime;
	    vgap *= (1 - delta / gapSupressionTime);
	  } else if (photonArriavaltime + gapSupressionTime < time*step &&
		     time*step < photonArriavaltime + gapSupressionTime + gapRecoveryTime*10) {
	    double delta = time*step - photonArriavaltime - gapSupressionTime;
	    vgap *= (1 - exp(delta / gapRecoveryTime));
	    double IR = 0, In = 0;
	    if (abs(V(x) - V(x+1)) >= Vgap) {
	      IR = (V(x) - V(x+1))/R; // Current through resistor shunt.
	      In = sqrt(2*T/R/step)*rnd.normal();
	    }
	    Is(x+1) = Ic_array(x)*(vgap/Vgap)*sin( theta(x) - theta(x+1) ) + IR + In;
	  }
	}
      }
      Is(Lx) = V(Lx-1)/Rterm + sqrt(2*T/Rterm/step)*rnd.normal(); // Current through right lead.

      // Define the equation system:
      {
	double Cm = 5*C;
	const double dt = step;

	int N = Lx;
	double L[N-1], D[N], U[N-1]; // Allocate the matrix tridiag(L,D,U) on the stack!

	// Fill in the capacitance and conductance matrices.
	for (int x = 0; x < N; x++) {

	  if (circuit_layout == comb && Cm > 0 && x % 10 == 5)
	    // x % 10 == 5 means starting from junction number 5 with periodicity of 289  XXXXXXXXXXXXXXXXXXXXXXXXXXX
	    D[x] = Cm;		// Capacitor comb every tenth.
	  else if (circuit_layout == all_different)
	    D[x] = C0_array[x];	// Specify all capacitors individually.
	  else
	    D[x] = C0;		// ordinary circuit.

	  double M = C;  // (Euler)

	  if (x+1 < N) { // Right neighbour:
	    if (circuit_layout == all_different) M = C = C_array[x];
	    if (abs(V(x) - V(x+1)) >= Vgap) M = C + dt/2/R; // (Symmetric)
	    D[x] += M;
	    U[x] = - M;
	  }

	  M = C;	// (Euler)

	  if (x > 0) { // Left neighbour:
	    if (circuit_layout == all_different) M = C = C_array[x-1];
	    if (abs(V(x-1) - V(x)) >= Vgap) M = C + dt/2/R; // (Symmetric)
	    D[x] += M;
	    L[x-1] = - M;
	  }
	}
	D[0]   += dt/2/Rterm; // First junction connected via Rterm to voltage source.
	D[N-1] += dt/2/Rterm; // Last junction terminated by Rterm to ground.

	// Terminal capacitance:
	if (Cterm > 0) { D[0] += Cterm - C0; D[N-1] += Cterm - C0; }

	// Middle capacitor if nonzero:
	if (circuit_layout == midcap && Cm > 0) D[N/2] += Cm - C0; // Replace C0 by Cm.

	// Calculate phase update:
	for (int x = 0; x < Lx; x++) {
	  // The right hand side of the equation system = div (Is + Ir + In)
	  new_theta(x) = Is(x+1) - Is(x); // = div I
	}
	// new_theta(Lx-2) -= (RC + step/2)*dU; Yes, but dU = dU/dt *dt = 0 for constant applied voltage.

	// Solve the equation system using tridiagonal LU solver:
	tri_solve(L,D,U,N,&new_theta[0]);
	// new_theta now contains the solution.
      }

#ifdef FIND_PHASE_SLIPS
      // Leap-frog: This version tries to identify phase slips. But only in bulk.
      for (int i = 0; i < 1; i++) {
	V[i] += - new_theta[i]*step;
	new_theta[i] = theta[i] + V[i]*step;
      }
      for (int i = 1; i < Lx; i++) {
	V[i] += - new_theta[i]*step;
	new_theta[i] = theta[i] + V[i]*step;
	double dtheta = theta[i]-theta[i-1];
	int n = lrint(dtheta/twopi);
	dtheta = new_theta[i]-new_theta[i-1];
	if (n != lrint(dtheta/twopi))
	  phout << time _ i << endl;
      }
      for (int i = Lx; i < Lx; i++) {
	V[i] += - new_theta[i]*step;
	new_theta[i] = theta[i] + V[i]*step;
      }
#else
      // Leap-frog:
      for (int i = 0; i < Lx; i++) {
	V[i] += - new_theta[i]*step;
	new_theta[i] = theta[i] + V[i]*step;
      }
#endif
      // The current going into the array from the left, i.e., at the
      // point where the voltage is applied:
      J += (U - V[0])/Rterm*step; // No need to include the noise since it averages to zero.

#if 0
      // Alternatively: average over the array.
      J += Is(0)*step;
      for (int x = 1; x < Lx; x++) { // Step through the Lx-1 junctions.
	// Correct resistive current:
	Is(x) += ( (newtheta(x-1) - 2*theta(x-1) - old_theta(x-1)) -
		   (newtheta(x) - 2*theta(x) - old_theta(x)) )/(2*R*step);
	// Add capacitative current:
	Is(x) += C*( (newtheta(x-1) - 2*theta(x-1) - old_theta(x-1)) -
		     (newtheta(x) - 2*theta(x) - old_theta(x)) )/sqr(step);
	J += Is(x)*step;
      }
#endif

      // Do the update:
      swap(theta.v, new_theta.v);
      swap(old_theta.v, new_theta.v);

      time ++;
    }

  } // rcsj

  virtual void equil(int ant, char *filnamn = NULL) {
    if (win) win->ShowTitle("equil");
    System::reset();
    rcsj(ant);
  }

  void ramp_voltage_slowly(double time_, double U_end) {
    System::reset();
    int num = time_/step;
    if (num == 0) num = 1;
    double dU = (U_end - U)/num;
    double Vo = Vac, oo = omega;
    Vac = 0.0; // Keep AC component zero during ramp.
    for (int n = 0; n < num; n++) {
      U += dU;
      rcsj();
      if (win && win->timetoshow()) showsystem();
      if (please_leave) { please_leave = 0; break; }
    }
    Vac = Vo;
    omega = oo;
    U = U_end;
  }

  virtual void run(int ant, int _sweeps) {
    if (win) win->ShowTitle(c_to_str("T = %g, C = %g",T,C));
    reset();
    System::reset();
    System:: sweeps = _sweeps;
    ant /= sweeps;
    for (int i = 0; i < ant; i++) {
      rcsj(sweeps);
      samp();
      if (win && win->timetoshow()) showsystem();
      if (please_leave) { please_leave = 0; break; }
    }
    results();
  }

  virtual void samp() {}
  virtual void results() {}
  virtual void write_data() {}

  void showsystem() {
    sighold(SIGHUP);
    if (!win) { sigrelse(SIGHUP); return; }

    double size = 0.5;

    win-> setview(Lx);

#if 0
    for (int x = 0; x < Lx; x++) { // Plot the XY spin:

      double sx = cos(theta(x)), sy = sin(theta(x));
      sx *= size;
      sy *= size;
      double xx = x + sx*0.5, yy = 2 + 0*Lx/2 + sy*0.5;
      win->arrow(xx-sx,yy-sy,xx,yy);

    }
#endif

    // Plot the potential:
    for (int x = 0; x < Lx; x++)
      win-> line(x,0,x,V[x]/U * Lx/2);
      // win-> line(x,Lx*3/4.0,x,Lx*3/4 + Is[x] * Lx/100); // Current

    for (int x = 0; x < Lx-1; x++)
      win-> line(x,Lx*3/4.0,x,Lx*3/4 + (V[x]-V[x+1]) * Lx/8, win-> gcred); // potential difference

    if (win->info) {
      char s[200];
      sprintf(s,
	      "T = %g, U = %g\n"
	      "Lx = %d\n"
	      "time: %d\n"
	      "step: %g\n"
	      "C: %g\n",
	      T, U,
	      Lx,
	      time, step,
	      C);
      win->Info(5,s);
    }

    win -> Draw();

    while (win -> Events()) {
      if (win-> key == 'X')
	win-> please_close = 1;
      if (win -> key == 'F')
	please_leave = 1;
    }
    if (win->please_close) { delete win; win = 0; }
    sigrelse(SIGHUP);
  }

}; // System

class SampData : public System
{
public:

  // int time;
  Estimate curr;
  Estimate resistivity;

  SampData() {}

  void reset() {
    time = 0;

    curr.reset();
    resistivity.reset();
  }

  inline void samp() {
    curr += J;
    resistivity	+= sqr(J);

    J = 0;

#ifdef PRINT_VOLTAGE
    ofstream out(c_to_str("Vx_%6.4f",U),ios::app);
    out.Bwrite(V.v,Lx);
#endif
  }

  void results() {
    curr.results();
    resistivity.results();

    // Final results are obtaind as
    //     resistivity	/= sweeps*step*Lx*Ly*T;

  }

  void write_data() {
    rawfile <<
      T _ U _ Lx _
      curr _
      resistivity _
      0 _
      time _ step _ sweeps << endl << flush;

    datafile <<
      T << U << Lx <<
      curr <<
      resistivity <<
      0.0 <<
      time << step << sweeps << flush;
  }

};

static int repl_max=0, rr_L;

extern "C"
void handler(int s) {
  if (s==SIGUSR1) { please_leave = 1; return; }
  if (win) { delete win; win = 0; }
  else { win = new VortexWin(); win->setview(rr_L); win->replmax=repl_max; }
}

bool System :: setparam(istream& in) {
  string name;
  char c;
  double val;

  in >> name >> c;
  if (c != '=') return false;

  if (name == "U") in >> U;
  else if (name == "Vac") in >> Vac;
  else if (name == "omega") in >> omega;
  else if (name == "f") { in >> omega; omega *= twopi; }
  else if (name == "T") in >> T;
  else if (name == "dt") in >> step;
  else if (name == "C") in >> C;
  else if (name == "C0") in >> C0;
  else if (name == "R") in >> R;
  else if (name == "Rterm") in >> Rterm;
  else if (name == "Cterm") in >> Cterm;
  else if (name == "Ic") { in >> Ic; set_Ic(); }
  else if (name == "Icweak") { in >> Icweak; set_Ic(); }

  else if (name == "Vgap") in >> Vgap;
  else if (name == "circuit_layout") { circuit_layout = read_layout(in); set_Ic(); }
  else if (name == "init_voltages") in >> boolalpha >> ::init_voltages;
  // else if (name == "current_bias") in >> current_bias;

  else if (name == "L") { in >> Lx; init(Lx); }
  else if (name == "Ic_array") { read_Ic(in); }
  else if (name == "C0_array") { read_C0(in); }
  else if (name == "C_array") { read_C(in); }
  else if (name == "photon") { read_photon_param(in); }

  else cerr << "What is " << name << '?' << endl;

  return true;
}

void System :: run_T_sweep(istream& cmd) {
  rr_L = int(Lx);

  double T_step, T_end;

  cmd >> T >> T_step >> T_end >> newln;
  int nT = rint((T_end-T)/T_step);

  for (int n = 0; n < nT; n++) {
    setup(T, U);
    equil(equ);
    run(number,sweeps);
    write_data();
    T += T_step;
  }
}

void System :: run_U_sweep(istream& cmd) {
  rr_L = int(Lx);

  double U_start, U_step, U_end;
  cmd >> U_start >> U_step >> U_end;
  int nU = rint((U_end-U_start)/U_step);

  setup(T, U);
  for (int n = 0; n < nU; n++) {
#ifdef FIND_PHASE_SLIPS
    phout.close();
    phout.open(c_to_str("PS_%6.4f",U_start));
// CHANGE FILE NAME "PS_%8.3f",U_start
#endif
    ramp_voltage_slowly(abs(U_step)*1000, U_start);
    equil(equ);
    run(number,sweeps);
    write_data();
    U_start += U_step;
  }
}

int System :: main(int argc, char *argv[]) {
  // Do some initializations:
  set_defaults();

  // Open output files:
  rawfile.open("raw", ios::app); // append everything to raw.
  rawfile.precision(12);
  datafile.open("data", ios::app);	// Binary data file.
#ifdef FIND_PHASE_SLIPS
    phout.open("phase-slips");
#endif



  istream & cmd = cin;
  char c = 0;

  while ( cmd >> ws >> c ) {

    if (c == '#') {
      cmd >> newln; // a comment.
      continue;
    }

    if (c == '!') { // Execute shell command.
      string s, l;
      getline(cmd,s);
      // Allow an ending backslash to continue lines:
      while (s[s.size()-1] == '\\' && cmd >> c && c == '!') {
	s[s.size()-1] = '\n';
	getline(cmd,l);
	s.append(l);
      }
      s.append("\n");
      system(s.c_str());
      continue;
    }

    if (c == 'H' || c == 'h') {
      cerr << rcsid << endl << endl << usage;
      continue;
    }

    if (c == 'z') {
      unsigned int seed;
      cmd >> seed;
      rnd.setseed(seed);
      cerr << "Warning: Using random seed " << seed << endl;
      continue;
    }

    // else if (c == 'J') { // Set electric current.
    //   cmd >> J;
    //   continue;
    // }

    // else if (c == 'U') { // Set applied volatage.
    //   cmd >> U;
    //   continue;
    // }

    else if (c == ':') {
      setparam(cmd);
      continue;
    }

    else if (c == 'r') {
      int L;
      cmd >> L >> number >> sweeps >> equ;
      if (L != Lx) init(L);
      run_T_sweep(cmd);
      continue;
    }

    else if (c == 'u') {
      int L;
      cmd >> L >> number >> sweeps >> equ;
      if (L != Lx) init(L);
      run_U_sweep(cmd);
      continue;
    }

    else if (c == 'X')	{
      handler(SIGHUP);
      continue;
    }

    else if (c == 'Q')		// Quit
      break;

    else if (! cmd)
      break;

    else {
      cmd >> newln;
      cerr << "What " << c << '?' << endl;
      continue;
    }

  }
  return 0;
}

int main(int argc, char *argv[]) {
  Xcon _xcon_(argv[0],0);		// Create X-Window.
  _xcon_.count--;
  rr_L = 20;
  sigset(SIGHUP,handler); repl_max=0;
  sigset(SIGUSR1,handler);
  // handler(SIGHUP);

  SampData Sim;
  Sim.main(argc, argv);

  if (win) delete win; win = 0;
  return(0);
}
