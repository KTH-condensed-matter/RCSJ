//-*-C++-*-	$Id: util.h,v 1.2 1998/07/01 12:01:08 jack Exp $

#ifndef util_H
#define util_H

using namespace std;

#include <iostream>
#include <fstream>
#include <cmath>

const double pi = M_PI;
const double twopi = 2*pi;

#define Sqr(x) ((x)*(x))
#define _ << ' ' <<
#define __ << '\t' <<

// #ifndef hppa
// inline double sqr(const double x) { return x*x; }
// #endif

template <class T> inline T sqr(const T x) { return x*x; }

// template <class T> inline T min(const T x,const T y)
// { return (x < y) ? x : y; }
// template <class T> inline T max(const T x,const T y)
// { return (x > y) ? x : y; }

// template <class T> inline void swap(T& x,T& y)
// { T temp = x; x = y; y = temp; }

inline istream & newln(istream & s) // skip to new line.
{ s.ignore(400,'\n'); return s; }

inline istream & skipcomm(istream & s) { // Skip comments.
  while (s.peek() == '#')
    s.ignore(400,'\n');
  return s;
}

#define Bread(v,n)  read ((char *) v, (n)*sizeof(*v))
#define Bwrite(v,n) write((char *) v, (n)*sizeof(*v))

class obstream : public ofstream {
public:
  obstream() : ofstream() {}
  obstream(const char f[]) : ofstream(f) {}
  obstream(const char f[], const ios::openmode m) : ofstream(f,m) {}
  obstream& operator<<(const float x)  { Bwrite(&x,1); return *this; }
  obstream& operator<<(const double x) { Bwrite(&x,1); return *this; }
  obstream& operator<<(const short x)  { Bwrite(&x,1); return *this; }
  obstream& operator<<(const int x)    { Bwrite(&x,1); return *this; }
  obstream& operator<<(const long x)   { Bwrite(&x,1); return *this; }
  obstream& operator<<(ostream& (*f)(ostream&)) { f(*this); return *this; }
};

class ibstream : public ifstream {
public:
  ibstream() : ifstream() {}
  ibstream(const char f[]) : ifstream(f) {}
  ibstream& operator>>(float& x)  { Bread(&x,1); return *this; }
  ibstream& operator>>(double& x) { Bread(&x,1); return *this; }
  ibstream& operator>>(short& x)  { Bread(&x,1); return *this; }
  ibstream& operator>>(int& x)    { Bread(&x,1); return *this; }
  ibstream& operator>>(long& x)   { Bread(&x,1); return *this; }
};

#include <stdarg.h>
#include <stdio.h>
inline const char* c_to_str(const char* str,...) {
  va_list args;
  va_start(args,str);
  static char s[64];
  vsprintf(s,str,args);
  return s;
}

#endif
