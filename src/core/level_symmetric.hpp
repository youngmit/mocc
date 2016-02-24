#pragma once

#include <vector>
#include <math.h>
#include <iostream>

#include "angle.hpp"
#include "global_config.hpp"
#include "error.hpp"
#include "blitz_typedefs.hpp"

using namespace mocc;

// Hide all of these constants from the rest of the world
namespace {
    real_t mu_base[] = {0.577350269189626, 0.350021000000000,
        0.266636000000000, 0.218218218218218, 0.192450089729876,
        0.174077655955702, 0.161575000000000, 0.149071198499989};

    real_t w_unique[] = {
        1.0,
        1.0/3.0,
        0.1761262, 0.1572071,
        0.1209876, 0.0907408, 0.0925925,
        0.0893043, 0.0725281, 0.0450455, 0.0539274,
        0.0707734, 0.0558760, 0.0373436, 0.0502654, 0.0258553,
        0.0580031, 0.0488943, 0.0228095, 0.0393955, 0.0380920, 0.0258382, 0.0082759,
        0.0489967, 0.0413235, 0.0203158, 0.0265468, 0.0378883, 0.0135404, 0.0326129, 0.0103825
    };

    int w_offset[] = {0, 1, 2, 4, 7, 11, 16, 23};

    int w_map[] = {
        1,
        1, 1, 1,
        1, 2, 2, 1, 2, 1,
        1, 2, 2, 2, 3, 2, 1, 2, 2, 1,
        1, 2, 2, 3, 4, 3, 2, 4, 4, 2, 1, 2, 3, 2, 1,
        1, 2, 2, 3, 5, 3, 4, 6, 6, 4, 3, 6, 7, 6, 3, 2, 5, 6, 6, 5, 2, 1, 2, 3, 4, 3, 2, 1,
        1, 2, 2, 3, 5, 3, 4, 6, 6, 4, 4, 7, 8, 7, 4, 3, 6, 8, 8, 6, 3, 2, 5, 6, 7, 6, 5, 2, 1, 2, 3, 4, 4, 3, 2, 1
    };

    int w_map_offset[] = { 0, 1, 4, 10, 20, 35, 63 };
}



// Produce a vector of angles matching the level-symmetric quadrature of order
// 'order'
std::vector<Angle> GenSn( int order ){

    if( order%2 != 0 ) {
        throw EXCEPT("Sn quadrature order must be even.");
    }
    if( order > 16 ) {
        throw EXCEPT("Max supported Sn quadrature order is 16.");
    }

    // n is the number of base cosines to use
    const int n = order/2;

    // set up the list of base cosines
    VecF mu;
    mu.push_back(mu_base[n-1]);
    if( order > 2) {
        const real_t delta_mu = 2.0 * ( 1.0 -
                3.0*(mu[0]*mu[0]) ) / (real_t)(order-2);
        for (int i=1; i<n; i++) {
            mu.push_back( sqrt(mu[0]*mu[0] + i*delta_mu) );
        }
    }

    // Alias the w_unique array to get a slice for the order we are interested
    // in.
    real_t* weights = &w_unique[w_offset[n-1]];
    // Alias into the w_map to get our indices
    int* map = &w_map[w_map_offset[n-1]];

    // Apply the permutations of the base cosines to get actual angles. We will
    // do this once for the first octant, then reflect around.
    std::vector<Angle> angles;
    real_t wsum = 0.0;
    int k=0;
    for( int i=0; i<n; i++ ) {
        for( int j=0; j<=i; j++ ) {
            Angle angle( mu[i-j],
                         mu[j],
                         mu[n-i-1],
                         weights[map[k]-1]);

            angles.push_back(angle);
            wsum += angle.weight;
            k++;
        }
    }

    // One of these days, i should calculate the weights algorithmically, rather
    // than storing in a table. For now, make sure the angular integral comes
    // out to 4*PI. We defined the angles above for one octant, so make sure
    // that the weights sum to 1 to machine precision.
    for( auto &a: angles ) {
        a.weight /= wsum;
    }

    return angles;

}

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


// Produce a vector of angles matching the Chebyshev-Gaussian quadrature of
// order 'azi-order' for azimuthal angles and 'polar-order' for polar
// angles.
   std::vector<Angle> GenProduct(const std::vector<std::pair<real_t,real_t>> 
           &azi, const std::vector<std::pair<real_t,real_t>> &pol){
       std::vector<Angle> angles;

   return angles;
}


// Produce a vector of angles matching the Chebyshev-Yamamoto quadrature of
// order 'azimuthal-order' for azimuthal angles and 'polar-order' for polar
// angles. 
