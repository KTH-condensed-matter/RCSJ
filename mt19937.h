/* Mersenne Twister random number generator in -*-C++-*-

   MT19937: Real number version ([0,1)-interval, etc.)
   Adapted to C++ by Jack Lidmar.

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura.       
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
   FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
   OF THE POSSIBILITY OF SUCH DAMAGE.

   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp

*/

// REFERENCE:
//
// M. Matsumoto and T. Nishimura,
// "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
// Pseudo-Random Number Generator",
// ACM Transactions on Modeling and Computer Simulation,
// Vol. 8, No. 1, January 1998, pp 3--30.

// Normal distribution from
// George Marsaglia and Wai Wan Tsang,
// "The Monty Python Method for Generating Random Variables"
// ACM Transactions on Mathematical Software,
// Vol. 24, No. 3, September 1998, Pages 341-350.

#ifndef Random_H
#define Random_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>

/* Period parameters */  
#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

class mt19937
{
protected:

  uint32_t mt[MT_N+1];		/* the array for the state vector  */
  intptr_t mti;			/* mti==N+1 means mt[N] is not initialized */

public:

  mt19937() : mti(MT_N+1) {}

  void init() {
    using namespace std;
    ifstream in("/dev/urandom");
    in.read((char*) mt, MT_N*sizeof(uint32_t));
    if (!in) {
      cerr << "Cannot open /dev/urandom!" << endl;
      abort();
    }
    for (intptr_t i = 1; i <= MT_N; i++) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941U))
	- i; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
    }
    mt[0] |= 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 

    mti = MT_N;
    generate();
  }

  /* initializing the array with a seed */
  uint32_t setseed(uint32_t seed) {
    mt[0]= seed & 0xffffffffU;
    for (intptr_t i=1; i<MT_N; i++) {
      mt[i] = (1812433253U * (mt[i-1] ^ (mt[i-1] >> 30)) + i);
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      mt[i] &= 0xffffffffUL; /* for >32 bit machines */
    }
    mti = MT_N;
    return mt[MT_N-1];
  }

  /* initialize by an array with array-length */
  /* init_key is the array for initializing keys */
  /* key_length is its length */
  void init_by_array(uint32_t init_key[], uint32_t key_length) {
    setseed(19650218UL);
    intptr_t i=1, j=0;
    int k = (MT_N > key_length ? MT_N : key_length);
    for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525U))
	+ init_key[j] + j; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++; j++;
      if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
      if (j>=key_length) j=0;
    }
    for (k=MT_N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941U))
	- i; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
    }
    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
    mti = MT_N;
  }

  void generate() {
    uint32_t y;
    uint32_t mag01[2]={0x0, MATRIX_A};

    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    intptr_t kk;
    if (mti == MT_N+1)   // if setseed() has not been called.
      init();

    for (kk=0;kk<MT_N-MT_M;kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (;kk<MT_N-1;kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (mt[MT_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1];

    mti = 0;
  }

  inline double operator()() {	// Return real in [0,1)
    if (mti >= MT_N) generate();
    uint32_t y;
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    return ( (double)y * 2.3283064365386963e-10 ); /* reals: [0,1)-interval */
  }

  inline double closed() {	// Return real in [0,1]
    if (mti >= MT_N) generate();
    uint32_t y;
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    return ( (double)y * 2.3283064370807974e-10 ); /* reals */
  }

  inline double open() {	// Return real in (0,1)
    if (mti >= MT_N) generate();
    uint32_t y;
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    return ( ((double) y+1.0) * 2.3283064359965952e-10 ); /* reals */
  }

  inline uint32_t uint32() {
    if (mti >= MT_N) generate();
    uint32_t y;
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    return y;
  }

  /* generates a random number on [0,1) with 53-bit resolution*/
  double double53() {
    uint32_t a = uint32() >> 5, b = uint32() >> 6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
  }

  double normal() {
    double x,y,v;
    x = ((int32_t) uint32())*1.167239e-9;
    if (fabs(x)<1.17741) return(x);
    y = uint32()*2.328306e-10;
    if (log(y) < .6931472-.5*(x*x)) return(x);
    x = (x>0) ? .8857913*(2.506628-x) : -.8857913*(2.506628+x);
    if (log(1.8857913-y) < .5718733-.5*(x*x)) return(x);
    do {
      v = ((int32_t) uint32())*4.656613e-10;
      x = -log(fabs(v))*.3989423;
      y = -log(uint32()*2.328306e-10);
    } while(y+y<x*x);
    return( v>0? 2.506628+x : -2.506628-x );
  }

};

typedef mt19937 Random;

static Random rnd;

#endif
