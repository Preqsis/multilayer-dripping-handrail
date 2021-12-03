#ifndef DISTRIBUTOR_HPP
#define DISTRIBUTOR_HPP

#include "BlobScheduler.hpp"

class Distributor {
private:
    double _q;
    double _zc;
    double*** _grid;
    std::vector<size_t> _dim;   // comms dimensions
    std::vector<double> _rProfile;      // rotation profile
    std::vector<double> _tProfile;      // temperature profile
    BlobScheduler* _scheduler;
    bool _hasScheduler;
public:
    Distributor();

    Distributor(double*** grid, std::vector<size_t> dim);

    Distributor(double*** grid, std::vector<size_t> dim, double q);

    Distributor(double*** grid, std::vector<size_t> dim, double q, BlobScheduler* scheduler);

    ~Distributor();

    void setBlobScheduler(BlobScheduler* scheduler);

    void setRotationProfile(std::vector<double> profile);

    void setTemperatureProfile(std::vector<double> profile);

    void setInflux(double q);

    double getRandW();

    void run (uint s);
};

#endif
