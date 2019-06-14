#ifndef __ME_HH__
#define __ME_HH__

#include "constants.hh"

// gg -> gg
//---------
inline double Mgg2gg(double s, double t, double u){
  return 0.5 * 9.0/2.0*(3.0-(u*t)/(s*s)-(u*s)/(t*t)-(s*t)/(u*u));
}

// gg -> qqbar
//------------
inline double Mgg2qqbar(double s, double t, double u){
  return Nf*(1.0/6.0*(u/t+t/u) - 3.0/8.0*(u*u+t*t)/(s*s));
}

// qg -> qg
//---------
inline double Mqg2qg(double s, double t, double u){
  return (u*u+s*s)/(t*t) - 4.0/9.0*(u/s+s/u);
}

// qq' -> qq'
//-----------
inline double Mqqp2qqp(double s, double t, double u){
  return 4.0/9.0*(s*s+u*u)/(t*t);
}

// qqbar -> q'qbar'
//-----------------
inline double Mqqbar2qpqbarp(double s, double t, double u){
  return (Nf-1)*4.0/9.0*(t*t+u*u)/(s*s);
}

// qq -> qq
//---------
inline double Mqq2qq(double s, double t, double u){
  return 0.5*(4.0/9.0*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u)) - 8.0/27.0*(s*s)/(t*u));
}

// qqbar -> qqbar
//---------------
inline double Mqqbar2qqbar(double s, double t, double u){
  return 4.0/9.0*((t*t+u*u)/(s*s)+(s*s+u*u)/(t*t)) - 8.0/27.0*(u*u)/(s*t);
}

// qqbar -> gg
//------------
inline double Mqqbar2gg(double s, double t, double u){
  return 0.5*(32.0/27.0*(u/t+t/u) - 8.0/3.0*(t*t+u*u)/(s*s));
}

#endif  // __ME_HH__
