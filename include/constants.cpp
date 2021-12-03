#ifndef CONSTANTS_CPP
#define CONSTANTS_CPP

namespace Constants {
    namespace mpi {
        // rank 
        const int MASTER    = 0;

        // compute tags
        const int COMPUTE   = 1;
        const int STOP      = 2;
        const int SKIP      = 3;
    }

    //const double h      = 6.62607004e-34; // m^2 * kg * s^-1
    const double h      = 6.62607004e-27; // cgs
    
    //const double k      = 5.670374419e-8; // W * m^-2 * K^-34
    const double k      = 5.670374419e-5; // cgs

    //const double c      = 2.99792458e8; // m * s^-1
    const double c      = 2.99792458e10; // cgs
    
    //const double m_sun  = 1.9891e30; // kg
    const double m_sun  = 1.9891e33; // g
    
    //const double G      = 6.6743e-11;
    const double G      = 6.67259e-12; // cgs
} 

#endif
