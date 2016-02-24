#pragma once

#include <vector>
#include <math.h>
#include <iostream>

#include "angle.hpp"
#include "global_config.hpp"
#include "error.hpp"
#include "blitz_typedefs.hpp"

using namespace mocc;

// Produce a vector of <theta,weight> pairs of size n_polar with Yamamoto
// quadrature within (0, PI/2). All weights sum to 1. Currently, only npol=3 
// is supported. 
std::vector<std::pair<real_t,real_t>> GenYamamoto( int n_polar ){

    std::vector<std::pair<real_t,real_t>> thetaWeightPairVec;
    
    if ( n_polar != 3 ) {
        throw EXCEPT("Only support Yamamoto quadrature when npol=3");
    }
    thetaWeightPairVec.emplace_back(0.167429147795000,4.623300000000000E-002);
    thetaWeightPairVec.emplace_back(0.567715121084000,0.283619000000000);
    thetaWeightPairVec.emplace_back(1.20253314678900,0.670148000000000);
    
    return thetaWeightPairVec;
}

// Produce a vector of <alpha,weight> pairs of size n_azimuthal with Chebyshev 
// quadrature within (0, PI/2). All weights sum to 1.
std::vector<std::pair<real_t,real_t>> GenChebyshev( int n_azimuthal ) {

    std::vector<std::pair<real_t,real_t>> alphaWeightPairVec;
    real_t weight = 1.0/n_azimuthal;
    real_t alpha;
    real_t delAlpha = 0.5*PI/(2*n_azimuthal);

    for( int i=0; i<n_azimuthal; i++ ) {
        alpha = delAlpha*(2*i+1);
        alphaWeightPairVec.emplace_back(alpha,weight);
    }
  
    return alphaWeightPairVec;
}

// Produce a vector of <theta,weight> pairs of size n_polar with Gaussian 
// quadrature within (0, PI/2). All weights sum to 1.
std::vector<std::pair<real_t,real_t>> GenGauss( int n_polar ){
    std::vector<std::pair<real_t,real_t>> thetaWeightPairVec;
    int N = n_polar*2 - 1;
    int N1 = n_polar*2;
    int N2 = n_polar*2 + 1;
    
    ArrayB1 xu(N1),y(N1),w(N1);
    real_t delxu=2.0/N;
    
    // Initial guess for y
    for( int i=0; i<N1; i++ ) {
        xu(i)=-1.0+i*delxu;
        y(i)=cos( (2*i+1)*PI/(2*N2) ) + 0.27/N1*sin(PI*xu(i)*N/N2);
    }

    // Legendre-Gauss Vandermonde Matrix
    ArrayB2 L(N1,N2);

    // Derivative of LGVM
    ArrayB1 Lg(N1);
    
    // Compute the zeros of the N+1 Legendre Polynomial using the recursion
    // relation and the Newton-Paphson method
    
    ArrayB1 y0(N1),Lp(N1);
    y0=2;

    // Iterate until new pioints are uniformly within epsilon of old points
    while( max(abs(y-y0))>FLOAT_EPS ) {
        L(blitz::Range::all(),0)=1;
        L(blitz::Range::all(),1)=y;

        for( int k=1; k<N2; k++ ) {
            L(blitz::Range::all(),k+1) = ( (2*k-1)*DotProduct(y,L(blitz::Range::all(),k))-
                    (k-1)*L(blitz::Range::all(),k-1) )/k;
        }

        for( int i=0; i<N1; i++ ) { 
            Lp(i) = N2*(L(i,N1)-DotProduct(y,L(blitz::Range::all(),N2)))/(1-y(i)*y(i));
        }

        y0=y;
        
        for( int i=0; i<N1; i++ ) {
            y(i) = y0(i)-L(i,N2)/Lp(i);
        }
    }
    
    for( int i=0; i<N1; i++ ) {
        w(i)=2.0/((1-y(i)*y(i))*Lp(i)*Lp(i))*N2*N2/(N1*N1);
    }

    for ( int i=n_polar; i<N1; i++ ) { 
        thetaWeightPairVec.emplace_back(y(i),w(i));
    }
    
    return thetaWeightPairVec;
}


// Produce a vector of angles with azi and pol pair vectors to represent a
// product quadrature set.
std::vector<Angle> GenProduct(const std::vector<std::pair<real_t,real_t>> 
           &azi, const std::vector<std::pair<real_t,real_t>> &pol){
       std::vector<Angle> angles;
       int n_azimuthal=azi.size();
       int n_polar=pol.size();
       real_t wsum=0.0;
       for( int i=0; i<n_azimuthal; i++ ) {
           for( int j=0; j<n_polar; j++ ) {
               Angle angle( azi[i].first, 
                            pol[i].first,
                            azi[i].second*pol[i].second );
               angles.push_back(angle);
               wsum += azi[i].second*pol[i].second;
           }
       }
       
       for( auto &a: angles ) {
           a.weight /= wsum;
       }

       return angles;
}
