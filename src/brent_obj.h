#ifndef BRENT_H
#define BRENT_H

//  Reference:
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//  Author:
//    Original FORTRAN77 version by Richard Brent.
//    C++ version by John Burkardt.
//
//  With only minor adjustment to avoid using static variables

# include <vector>
# include <cmath>
# include <cstdlib>
# include <ctime>
# include <iostream>
using namespace std;

class BrentObj {
public:
    double arg,eps,tol,tol1,tol2;
    double c,d,e;
    double fu,fv,fw,fx;
    double midpoint;
    double p,q,r,u,v,w,x;
    int32_t status;

BrentObj(double& a, double& b) {
    status = 0;
    if ( b <= a ) {
        printf("BrentObj::BrentObj - Fatal error: A < B is required\n");
        status = -1;
        exit ( 1 );
    }
    c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );
    eps = sqrt ( r8_epsilon ( ) );
    tol = r8_epsilon ( );
    v = a + c * ( b - a );
    w = v;
    x = v;
    e = 0.0;
    status = 1;
    arg = x;
}

double local_min_rc( double& a, double& b, double value ) {
    if ( status == 1 ) {
        //  STATUS = 1, initialize function value.
        fx = value;
        fv = fx;
        fw = fx;
    } else if ( status >= 2 ) {
        //  STATUS >= 2, update.
        fu = value;
        if ( fu <= fx ) {
            if ( x <= u ) {a = x;}
            else {b = x;}
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if ( u < x ) {a = u;}
            else{b = u;}
            if ( fu <= fw || w == x ) {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if ( fu <= fv || v == x || v == w ) {
                v = u;
                fv = fu;
            }
        }
    }
    //  Take the next step.
    midpoint = 0.5 * ( a + b );
    tol1 = eps * fabs ( x ) + tol / 3.0;
    tol2 = 2.0 * tol1;
    //  If the stopping criterion is satisfied.
    if ( fabs ( x - midpoint ) <= ( tol2 - 0.5 * ( b - a ) ) ) {
        status = 0;
        return arg;
    }
    //  Is golden-section necessary?
    if ( fabs ( e ) <= tol1 ) {
        if ( midpoint <= x ) {e = a - x;}
        else {e = b - x;}
        d = c * e;
    } else {
        //  Consider fitting a parabola.
        r = ( x - w ) * ( fx - fv );
        q = ( x - v ) * ( fx - fw );
        p = ( x - v ) * q - ( x - w ) * r;
        q = 2.0 * ( q - r );
        if ( 0.0 < q ) {p = - p;}
        q = fabs ( q );
        r = e;
        e = d;
        if ( ( fabs ( 0.5 * q * r ) <= fabs ( p ) ) || ( p <= q * ( a - x ) ) ||
             ( q * ( b - x ) <= p ) ) {
            //  Choose a golden-section step if the parabola is not advised.
            if ( midpoint <= x ) {e = a - x;}
            else {e = b - x;}
            d = c * e;
        } else {
            //  Choose a parabolic interpolation step.
            d = p / q;
            u = x + d;
            if ( ( u - a ) < tol2 ) {d = tol1 * r8_sign ( midpoint - x );}
            if ( ( b - u ) < tol2 ) {d = tol1 * r8_sign ( midpoint - x );}
        }
    }
    //  F must not be evaluated too close to X.
    if ( tol1 <= fabs ( d ) ) {u = x + d;}
    if ( fabs ( d ) < tol1 ) {u = x + tol1 * r8_sign ( d );}
    //  Request value of F(U).
    arg = u;
    status++;
    return arg;
    }

    double r8_epsilon () {
        const double value = 2.220446049250313E-016;
        return value;
    };
    double r8_max (double x, double y) {
        if ( y < x ) {return x;}
        else {return y;}
    };
    double r8_sign ( double x ) {
        if (x < 0.0) {return -1;}
        else {return 1.;}
    };
};

#endif
