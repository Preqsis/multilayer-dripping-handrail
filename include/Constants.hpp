#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace Constants {
    namespace mpi {
        // rank 
        const int MASTER    = 0;

        // compute tags
        const int COMPUTE   = 1;
        const int STOP      = 2;
        const int SKIP      = 3;
    }

    // Planck's const
    const double h      = 6.6260755e-27; // cm^2 * g * s^-1
    
    // Boltzmann's constant
    const double k      = 1.380658e-16; // erg * K^-1

    // Speed of light
    const double c      = 2.99792458e10; // cm * s^-1
    
    // Sun mass
    const double m_sun  = 1.9891e33; // g

    // Sun radius
    const double r_sun  = 69634000000; // cm
    
    // Gravitational constant
    const double G      = 6.67259e-8; // cm^3 * g^-1 * s^-2

    // Stefan-Boltzmann constant
    const double sigma  = 5.670374e-5; // erg * cm^-2 * s^-1 * K^-4

    // 
    const double R      = 8.314e7; // erg * mol^-1 * K^-1

    const double M_h    = 1.00794; // g * mol^-1
} 

#endif
