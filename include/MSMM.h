#ifndef MSMM_H
#define MSMM_H

#include <iostream>
#include <cmath>
#include <boost/array.hpp>

typedef boost::array<double, 2> state_type;

class Model {
private:
protected:
    double _m; // hmotnost (jeji pocatecni hodnota)
    double _dm; // pritok
    double _v0; // rychlost pritoku
    double _g     = 1.0; // grav. zrychleni
    double _gamma = 0.05; // 'tuhost pruziny'
public:
    virtual void operator()(const state_type &y, state_type &dydx, double x) = 0;

    double get_k(double m) {
        if(m < 4.61)
            return -11.4 * m + 52.5;
        return 0.0;
    }

    // hmotnost
    void set_m(double m) {_m = m;}
    double get_m() {return _m;}

    // 'rychlost' zmeny hmotnosti (pritoku)
    void set_dm(double dm) {_dm = dm; }
    double get_dm() {return _dm;}

    // rychost pritoku
    void set_v0(double v0) {_v0 = v0;}
    double get_v0() {return _v0;}

    // grav. zrychleni
    void set_g(double g) {_g = g; }
    double get_g() {return _g;}
};

class MSMM : public Model {
private:
public:
    void operator()(const state_type &y, state_type &dydx, double x) {
        dydx[0] = y[1];
        dydx[1] = _g 
            - get_k(_m) * y[0] / _m 
            - _gamma * y[1] / _m 
            - (_dm * (y[1] - _v0)) / _m; 
    }
};

class MSMM2 : public Model {
private: 
public:
    void operator()(const state_type &y, state_type &dydx, double x) {
        dydx[0] = y[1];
        dydx[1] = _g 
            - get_k(_m) * y[0] / _m 
            - _gamma * y[1] / _m 
            - _dm * y[1] / _m;
    }
};

#endif

