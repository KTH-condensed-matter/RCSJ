// rcsj.cc - RCSJ simulation of a 1D XY model with shunting capacitors.
// 1D version.

// $Id: rcsj.cc,v 1.32 2013/06/25 12:15:46 jack Exp $

// #define CHECK_OUT_OF_BOUNDS 2

// #define RANDOM_INITIAL_CONDITION

#define FIND_PHASE_SLIPS
// #define PRINT_VOLTAGE

// Put a resistor Rterm also at the end of the array.
// #define END_RESISTOR
// Otherwise it is connected directly to ground and V(Lx-1) = theta(Lx-1) = 0.

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <array>

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
  " : param = ?			-- Set parameter.\n"
  " u Lx number sweeps equil Umin Ustep Umax		-- Run voltage sweep.\n"
  " r Lx number sweeps equil Tmin Tstep Tmax		-- Run temperature sweep.\n"
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

  double VL = 0;          // Voltage at the first node before bias tee.
  angle thetaL = 0;       // Phase at the first node before bias tee.
  double VA = 0;          // Voltage at the amplifier.

  Table<double> Is;		// Supercurrents (temporary)

  Table<double> Ic_array;	// Critical currents array of length Lx-1.
  Table<double> C0_array;	// Capaciatance to ground array of length Lx.
  Table<double> C_array;	// Capacitance of junctions array of length Lx-1.

  double T;			// Temperature.
  double J, Jtot, Jshunt;			// Electric currents (measured).
  double Ic, Icweak;		// Critical current of tunnel junctions.
  double U;			// Applied voltage.
  double C, C0, Cterm;		// Junction capacitance, and capacitance to ground.
  double R, Rterm;		// Junction resistance and terminal resistance.
  double Rqp;         // Quasiparticle resistance if > 0, otherwise oo.
  double Rshunt;      // Shunt resistance (for now between first and last SC island).
  double Lterm;       // Incuctance in bias tee. (The capacitance is Cterm.)
  double Ramp;        // Amplifier impedance. 50 Ohm in reduced units?
  double Camp;        // Capacator shunting Ramp to reduce noise.
  double Vgap;			// Gap voltage (or 0).
  double Vac, omega;		// Applied frequence dependent voltage.

  int normalJunctionNumber;	// Junction which goes normal when a photon hits.
  double gapSupressionTime;	// Time scale for gap supression.
  double gapRecoveryTime;	// Time scale for gap recovery.
  double photonArrivalTimeInterval;	// Time interval between photon arrivals.
  double photonArrivalTime;	// Time of next photon arrival.
  
  double Energy;		// Current energy (not used).

  int time;			// Time (number of md steps).
  int number, equ;		// Number of sweeps and sweeps to equilibrate.
  int sweeps;			// Time between measurements.
  double step;			// Step size for theta updates.

  Tri eqsys;

  ofstream rawfile;		// append everything to raw.
  obstream datafile;		// Binary data file.

  ofstream phout;		// Data file containing phase slip events (t, x).

  ofstream vout{"V"}; // Voltage at beginning of array V[0].

  double trigger_level = 0;
  int trigger_index = -1;
  array<double,8192> Vbuffer;  // Length needs to be a power of 2.
  int Vbuff_index = 0;

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
    Rqp = 0;        // 0 means oo.
    Rterm = R/100;	// 50 Ohms instead?
    Rshunt = 0;     // 0 means not present.
    Cterm = 0;			// 0 means same as C0.
    Lterm = 0;
    Ramp = 0;
    Camp = 0;
    Vgap = 0;
    Ic = 1.0;
    Icweak = 0;			// 0 means 0.1*Ic.
    U = 0.0;
    Vac = 0.0;			// No AC.
    omega = 0.0;		// -"-
    // Ic_array = Ic;		// Same critical current by default.

    normalJunctionNumber = Lx/2;
    gapSupressionTime = 0.0;
    gapRecoveryTime = 4.0; // XXX Units?
    photonArrivalTimeInterval = 0;
    photonArrivalTime = -100000.0;
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
    cmd >> normalJunctionNumber
	>> gapSupressionTime
	>> gapRecoveryTime
	>> photonArrivalTimeInterval;
    photonArrivalTime = -100000.0;
    return !!cmd;
  }

  double set_dt(istream& in) {
    double dt;
    string s;
    in >> s;
    if (s == "sqrt(LC)*") {
      in >> dt;
      double L_k = 1.0/Ic;     // Assumes hbar / 2e = 1...
      dt *= sqrt(L_k * C);
    }
    else
      dt = stod(s);
    return dt;
  }

#ifdef END_RESISTOR

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

#else

  void init_voltages() {	// Initialize voltages and phases to a reasonable state.
    if (U < Rterm*Ic) {   // Array is superconducting. V(x) = 0.
      angle dtheta = asin(U/(Rterm*Ic));
      theta(Lx-1) = 0.0;
      V(Lx-1) = 0;
      for (int x  = Lx-2; x >= 0; x--) {
        V(x) = 0;
        theta(x) = theta(x+1) + dtheta;
      }
    }
    else { // Might need some thinking...
      double I = U/(Rterm + Rshunt);      
      double Ua = I*Rshunt;
      double n = floor(Ua/R);
      double Rarr = n*R;
      double Rtot = Rterm + Rarr;
      if (Rshunt > 0)
        Rtot = Rterm + Rshunt * Rarr /(Rshunt + Rarr);
      double Iarr = (U-Rterm*I)/Rtot;
      double Itot = Iarr + Ic;
      for (int x = 0; x < Lx; x++)
        V(x) = U - Itot*Rterm - Iarr*R*x/n;
    }
  }
#endif

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
    J = Jtot = Jshunt = 0;
    time = 0;
  }

  inline int xplus(const int x)  { return (x+1 == Lx) ? 0 : x+1; }
  inline int xminus(const int x) { return (x == 0) ? Lx-1 : x-1; }

  // Lx = number of sc islands. number of JJ = Lx-1

  void rcsj(int ant = 1) {

    for (int iii = 0; iii < ant; iii++) {

      if (photonArrivalTimeInterval > 0 &&
          photonArrivalTime + 0.5*photonArrivalTimeInterval < time*step)
        photonArrivalTime = time*step + 0.5*photonArrivalTimeInterval;

      // Voltage bias: Delta is twisted with time!

      // Calculate all supercurrents + noise, and the Delta update.

      double Ut = U + Vac*sin(omega*time*step); // Applied voltage, possibly time dependent.

      // VL is the voltage at the node between the input resistors and the bias tee inductor.

      Is(0) = (Ut - VL)/Rterm + sqrt(2*T/Rterm/step)*rnd.normal(); // Current through left lead.
  
      // Maybe add a shunt to ground here:
      if (Rshunt > 0)
        Is(0) -= (VL - 0.0)/Rshunt + sqrt(2*T/Rshunt/step)*rnd.normal(); // Current through left lead.

      for (int x = 0; x < Lx-1; x++) { // Step through the Lx-1 junctions.
	
        double IR = 0, In = 0;
	      if (abs(V(x) - V(x+1)) >= Vgap) {
	        IR = (V(x) - V(x+1))/R; // Current through resistor shunt.
	        In = sqrt(2*T/R/step)*rnd.normal();
	      }
        else if (Rqp > 0) {
	        IR = (V(x) - V(x+1))/Rqp; // Quasiparticle current through resistor shunt.
	        In = sqrt(2*T/Rqp/step)*rnd.normal();
        }
      	Is(x+1) = Ic_array(x)*sin( theta(x) - theta(x+1) ) + IR + In;

      }

      double vgap = Vgap;
      if (photonArrivalTimeInterval > 0 && photonArrivalTime <= time*step) {
        int x = normalJunctionNumber;
        if (time*step < photonArrivalTime + gapSupressionTime) {
          double delta = time*step - photonArrivalTime;
          vgap *= (1 - delta / gapSupressionTime);
        } else if (photonArrivalTime + gapSupressionTime <= time*step &&
              time*step < photonArrivalTime + gapSupressionTime + gapRecoveryTime*10) {
          double delta = time*step - photonArrivalTime - gapSupressionTime;
          vgap *= (1 - exp(-delta / gapRecoveryTime));
        }
        double IR = 0, In = 0;
        // XXX Need to analyze the effect of reducing the gap...
        if (abs(V(x) - V(x+1)) >= vgap) {
          IR = (V(x) - V(x+1))/R; // Current through resistor shunt.
          In = sqrt(2*T/R/step)*rnd.normal();
        }
        else if (Rqp > 0) {
	        IR = (V(x) - V(x+1))/Rqp; // Quasiparticle current through resistor shunt.
	        In = sqrt(2*T/Rqp/step)*rnd.normal();
        }
        Is(x+1) = Ic_array(x)*(vgap/Vgap)*sin( theta(x) - theta(x+1) ) + IR + In;

        // Save Voltage at input, and the normal junction to file:
        save_voltage_vs_time(x);
      }

#ifdef END_RESISTOR
      Is(Lx) = V(Lx-1)/Rterm + sqrt(2*T/Rterm/step)*rnd.normal(); // Current through right lead.
#else
      Is(Lx) = 0.0; // Not used. Lx correspond to the ground. V(Lx-1) = 0 in this case.
#endif

      // Define the equation system:
      {
        double Cm = 5 * C;
        const double dt = step;

        int N = Lx;
        double L[N - 1], D[N], U[N - 1]; // Allocate the matrix tridiag(L,D,U) on the stack!

        // Fill in the capacitance and conductance matrices.
        for (int x = 0; x < N; x++) {

          if (circuit_layout == comb && Cm > 0 && x % 10 == 5)
            // x % 10 == 5 means starting from junction number 5 with periodicity of 289  XXXXXXXXXXXXXXXXXXXXXXXXXXX
            D[x] = Cm; // Capacitor comb every tenth.
          else if (circuit_layout == all_different)
            D[x] = C0_array[x]; // Specify all capacitors individually.
          else
            D[x] = C0; // ordinary circuit.

          double M = C; // (Euler)

          if (x + 1 < N) { // Right neighbour:
            if (circuit_layout == all_different)
              M = C = C_array[x];
            if (abs(V(x) - V(x + 1)) >= Vgap)
              M = C + dt / 2 / R; // (Symmetric)
            else if (Rqp > 0)
              M = C + dt / 2 / Rqp;
            D[x] += M;
            U[x] = -M;
          }

          M = C; // (Euler)

          if (x > 0) { // Left neighbour:
            if (circuit_layout == all_different)
              M = C = C_array[x - 1];
            if (abs(V(x - 1) - V(x)) >= Vgap)
              M = C + dt / 2 / R; // (Symmetric)
            else if (Rqp > 0)
              M = C + dt / 2 / Rqp;
            D[x] += M;
            L[x - 1] = -M;
          }
        }
        // In the bias tee configuration this must be disabled:
        if (Lterm > 0)
          D[0] += dt / 2 / Rterm;     // First junction connected via Rterm to voltage source.
#ifdef END_RESISTOR
        D[N - 1] += dt / 2 / Rterm; // Last junction terminated by Rterm to ground.
#endif

        // Handle photons here also.
        if (photonArrivalTimeInterval > 0 && vgap < Vgap) {
          int x = normalJunctionNumber;
          double dV = abs(V(x) - V(x + 1));
          if (dV >= vgap && dV < Vgap) {
            double M = dt / 2 / R;
            if (Rqp > 0) M -= dt / 2 / Rqp;
            D[x]   += 2*M;
            U[x]   -= M;
            L[x-1] -= M;
          }
        }

        // In the bias tee configuration the following two lines must be disabled.
        if (Lterm > 0 && Rshunt > 0) {
          D[0] += dt / 2 / Rshunt;     // First junction connected to ground.
        }

        // Terminal capacitance:
        if (Cterm > 0)
        {
          D[0] += Cterm - C0;
          D[N - 1] += Cterm - C0;
        }

        // Middle capacitor if nonzero:
        if (circuit_layout == midcap && Cm > 0)
          D[N / 2] += Cm - C0; // Replace C0 by Cm.

        // Calculate phase update:
        for (int x = 0; x < Lx; x++) {
          // The right hand side of the equation system = div (Is + Ir + In)
          new_theta(x) = Is(x + 1) - Is(x); // = div I
        }
#ifndef END_RESISTOR
        new_theta(Lx-1) = 0.0;
        N --; // Reduce the size of the equation system!!!
#endif

        // new_theta(Lx-2) -= (RC + step/2)*dU; Yes, but dU = dU/dt *dt = 0 for constant applied voltage.

        // Solve the equation system using tridiagonal LU solver:
        tri_solve(L, D, U, N, &new_theta[0]);
        // new_theta now contains the solution.
      }

      // Bias tee - Inductance:
      if (Lterm > 0) {
        double Gin = 1.0/Rterm + 1.0/Rshunt;
        VL += ( Is(0) - (thetaL - theta(0))/Lterm )/Gin;
        thetaL += step*VL;
      }

      // Bias tee - Capacitance:
      // The capacitance is Cterm, already accounted for.
      // Just add voltage over series resistor to ground:
      V[0] -= VA;
      if (Camp == 0)
        VA = Ramp * (Is(0) - Is(1)) + sqrt(2*T*Ramp/step)*rnd.normal(); // Voltage over amplifier.
      else
        VA += step/(Ramp * Camp + 0.5*step) * ( Ramp * (Is(0) - Is(1)) - VA + sqrt(2*T*Ramp/step)*rnd.normal() ); // Voltage over amplifier.
      V[0] += VA;

#ifdef FIND_PHASE_SLIPS
      // Leap-frog: This version tries to identify phase slips. But only in bulk.
      for (int i = 0; i < 1; i++) {
        V[i] += -new_theta[i] * step;
        new_theta[i] = theta[i] + V[i] * step;
      }
      for (int i = 1; i < Lx; i++) {
        V[i] += -new_theta[i] * step;
        new_theta[i] = theta[i] + V[i] * step;
        double dtheta = theta[i] - theta[i - 1];
        int n = lrint(dtheta / twopi);
        dtheta = new_theta[i] - new_theta[i - 1];
        if (n != lrint(dtheta / twopi))
          phout << time*step _ i << endl;
      }
      for (int i = Lx; i < Lx; i++) { // Never...
        V[i] += -new_theta[i] * step;
        new_theta[i] = theta[i] + V[i] * step;
      }
#else
      // Leap-frog:
      for (int i = 0; i < Lx; i++) {
        V[i] += -new_theta[i] * step;
        new_theta[i] = theta[i] + V[i] * step;
      }
#endif
      if (Lterm == 0)
        VL = V[0]; // Disable bias tee inductor.

      // The current going into the array from the left, i.e., at the
      // point where the voltage is applied:
      Jtot += (U - VL) / Rterm * step; // No need to include the noise since it averages to zero.
      Jshunt += VL / Rshunt * step; // No need to include the noise since it averages to zero.

      J += Is(1) * step; // Current through the first junction.

      trigger_oscilloscope(V[0]);

      // Do the update:
      swap(theta.v, new_theta.v);

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

  void save_voltage_vs_time(int x) {
    double t = time * step;
    vout.Bwrite(&t, 1);
    vout.Bwrite(&V[0], 1);
    vout.Bwrite(&V[x], 1);
    vout.Bwrite(&V[x + 1], 1);
  }

  void trigger_oscilloscope(double V) {
    if (trigger_level > 0) {
      Vbuffer[Vbuff_index++] = V;
      const int mask = Vbuffer.size()-1;
      Vbuff_index &= mask;
      if (trigger_index < 0 && V > trigger_level)
        trigger_index = (Vbuff_index -1 + Vbuffer.size() - 512) & mask;
      else if (trigger_index == Vbuff_index) {
        static ofstream out("Vtrig");
        for (int i = 0; i < Vbuffer.size(); i++)
          out << (i - 512) * step __ Vbuffer[(i + trigger_index) & mask] << endl;
        out << '&' << endl;
        trigger_index = -1;
      }
    }
  }

  virtual void samp() {}
  virtual void results() {}
  virtual void write_data() {}

  void showsystem() {
    sighold(SIGHUP);
    if (!win) { sigrelse(SIGHUP); return; }

    double size = 0.5;

    win-> setview(0,-2,Lx,Lx);

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

    static double Vscale = 1000;
    const int VN = Vbuffer.size(), mask = VN - 1;
    win-> line(0.0,0.0,Lx,0);
    for (int i = 0; i < VN-1; i++)
      win-> line(i*1.0*Lx/VN, Vbuffer[(i + Vbuff_index) % mask]*Vscale,(i+1.0)*Lx/VN, Vbuffer[(i + 1 + Vbuff_index) % mask]*Vscale, win-> gcgreen);

    static double ii = 0.0;
    ii += 0.001*(Is(0)/Ic-ii);

    if (win->info) {
      char s[200];
      sprintf(s,
	      "T = %g, U = %g\n"
	      "Lx = %d\n"
	      "time: %d\n"
	      "step: %g\n"
	      "C: %g\n"
        "I/Ic: %g\n",
	      T, U,
	      Lx,
	      time, step,
	      C,
        ii);
      win->Info(6,s);
    }

    win -> Draw();

    while (win -> Events()) {
      if (win-> key == 'X')
      	win-> please_close = 1;
      if (win-> key == 'F')
        please_leave = 1;
      if (win-> keysym == 0xff52) // up
        U += 1.0;
      if (win-> keysym == 0xff54) // down
        U -= 1.0;
      if (win-> keysym == 0xff55) // prior (fn-up)
        U += 0.1;
      if (win-> keysym == 0xff56) // next (fn-down)
        U -= 0.1;
      if (win-> mouse == 1) {
        normalJunctionNumber = int(win-> xc(win-> mouse_x));
        photonArrivalTime = time*step;
        photonArrivalTimeInterval = 1e12;  // Needs to be nonzero...
        clog << "=> Photon at x = " << normalJunctionNumber << ", time = " << photonArrivalTime << endl;
      }
      if (win-> mouse == 4)
        Vscale *= 1.2;
      if (win-> mouse == 5)
        Vscale /= 1.2;
    }
    if (win->please_close) { delete win; win = 0; }
    sigrelse(SIGHUP);
  }

}; // System

class SampData : public System
{
public:

  // int time;
  Estimate curr, jtot, jshunt;
  Estimate resistivity;
  Estimate voltage; // only over the array.

  SampData() {}

  void reset() {
    time = 0;

    curr.reset();
    jtot.reset();
    jshunt.reset();
    resistivity.reset();
    voltage.reset();
  }

  inline void samp() {
    curr += J;
    jtot += Jtot;
    jshunt += Jshunt;
    resistivity	+= sqr(J);
    voltage += V[0] - V[Lx-1];

    J = 0;
    Jtot = 0;
    Jshunt = 0;

#ifdef PRINT_VOLTAGE
    ofstream out(c_to_str("Vx_%6.4f",U),ios::app);
    out.Bwrite(V.v,Lx);
#endif
  }

  void results() {
    curr.results();
    jtot.results();
    jshunt.results();
    resistivity.results();
    voltage.results();

    // Final results are obtaind as
    //     resistivity	/= sweeps*step*Lx*Ly*T;

  }

  void write_data() {
    rawfile <<
      T _ U _ Lx _
      curr _
      jtot _
      jshunt _
      resistivity _
      voltage _
      time _ step _ sweeps << endl;

    datafile <<
      T << U << Lx <<
      curr <<
      jtot <<
      jshunt <<
      resistivity <<
      voltage <<
      time << step << sweeps;
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
  else if (name == "dt") step = set_dt(in);
  else if (name == "C") in >> C;
  else if (name == "C0") in >> C0;
  else if (name == "R") in >> R;
  else if (name == "Rqp") in >> Rqp;
  else if (name == "Rterm") in >> Rterm;
  else if (name == "Rshunt") in >> Rshunt;
  else if (name == "Cterm") in >> Cterm;
  else if (name == "Ic") { in >> Ic; set_Ic(); }
  else if (name == "Icweak") { in >> Icweak; set_Ic(); }
  else if (name == "Lterm") in >> Lterm;
  else if (name == "Ramp") in >> Ramp;
  else if (name == "Camp") in >> Camp;

  else if (name == "Vgap") in >> Vgap;
  else if (name == "circuit_layout") { circuit_layout = read_layout(in); set_Ic(); }
  else if (name == "init_voltages") in >> boolalpha >> ::init_voltages;
  // else if (name == "current_bias") in >> current_bias;

  else if (name == "L") { in >> Lx; init(Lx); }
  else if (name == "Ic_array") { read_Ic(in); }
  else if (name == "C0_array") { read_C0(in); }
  else if (name == "C_array") { read_C(in); }
  else if (name == "trigger_level") in >> trigger_level;
  else if (name == "photon") { read_photon_param(in); }

  else cerr << "What is " << name << '?' << endl;

  return true;
}

void System :: run_T_sweep(istream& cmd) {
  rr_L = int(Lx);

  double T_step, T_end;

  cmd >> T >> T_step >> T_end >> newln;
  int nT = 1 + rint((T_end-T)/T_step);

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
  int nU = 1 + rint((U_end-U_start)/U_step);

  setup(T, U);
  for (int n = 0; n < nU; n++) {
#ifdef FIND_PHASE_SLIPS
    phout.close();
    phout.open(c_to_str("PS_%6.4f",U_start));
// CHANGE FILE NAME "PS_%8.3f",U_start
    vout.close();
    vout.open(c_to_str("V_%6.4f",U_start)); // Set file name for voltage output.
#endif
    // ramp_voltage_slowly(abs(U_step)*1000, U_start);
    ramp_voltage_slowly(equ*step, U_start);
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
