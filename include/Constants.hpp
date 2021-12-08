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
    
    // Gravitational constant
    const double G      = 6.67259e-8; // cm^3 * g^-1 * s^-2
} 

#endif
