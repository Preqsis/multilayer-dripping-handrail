#ifndef DISTRIBUTOR_H
#define DISTRIBUTOR_H

#include <random>
#include "BlobScheduler.h"

class Distributor {
private:
    uint _idim;
    uint _jdim;
    uint _ic;
    uint _jc;
    uint _n;

    double _dx;
    double _q;
    double _zc;

    double** _drain;

    std::vector<size_t> _cdim;          // comms dimensions
    std::vector<size_t> _gdim;          // grid dimensions
    std::vector<double> _rProfile;      // rotation profile
    std::vector<double> _tProfile;      // temperature profile

    BlobScheduler* _scheduler;
    bool _hasScheduler = false;
public:
    Distributor(std::vector<size_t> gdim, std::vector<size_t> cdim, double** drain) {
        _gdim   = gdim;
        _cdim   = cdim;
        _drain  = drain;
        _zc     = 5.5;
        _q      = 0.5;
    }

    void setBlobScheduler(BlobScheduler* scheduler) {
        _scheduler      = scheduler;
        _hasScheduler   = true;
    }

    void setRotationProfile(std::vector<double> profile) {
        _rProfile.clear();
        _rProfile = profile;
    }

    void setTemperatureProfile(std::vector<double> profile) {
        _tProfile.clear();
        _tProfile = profile;
    }

    double getRandW() {
        std::random_device dev;
        std::default_random_engine engine(dev());
        std::uniform_real_distribution<double> dist(0.0, 0.3);
        return dist(engine);
    }

    void step (double** grid, uint n) {
        double mb;
        double w;
        double dp;
        uint k_out, k_in, k;
        uint k_lower_right, k_lower_left;
        double j_in, j_out;
        double part_lower_left, part_lower_right;

        // rotace vnejsiho prstence o "jeden pohyb"
        // zabranuje driftovani vliven nepresnosti datovych typu
        std::vector<double> tmp;
        for (uint j = 0; j < _gdim[1]; j++) {
            tmp.push_back(grid[j][9]);
        }
        std::rotate(tmp.begin(), tmp.begin() + 1, tmp.end());
        for (uint j = 0; j < _gdim[1]; j++) {
            grid[j][9] = tmp[j];
            grid[j][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi
        }

        // rotace ostatnich podle profilu
        for (uint i = 1; i < _gdim[0]; i++) {
            for (uint j = 0; j < _gdim[1]; j++) {
                k = i * _gdim[1] + j;
                grid[k][9] = std::fmod(grid[k][9] + _rProfile[i], 2.0 * M_PI);
                grid[k][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi
            }
        }

        // nalezeni pritokove bunky + pritok
        uint idx = 0;
        for (uint j = 0; j < _gdim[1]; j++) {
            if (grid[j][9] == 0.0) {
                idx = j;
                break;
            }
        }
        grid[idx][5] += _q;
        grid[idx][6] = _q; // u pritokove bunky presne odpovida _q, zbytek 0

        for (uint i=_gdim[0]-1; i < _gdim[0]; i--) { // od stredu ven, unsigned preskakuje na maximalni hodnotu!!!
            if (i < _gdim[0] - 1) {
                k_out   = i * _gdim[1] + 0;
                k_in    = (i + 1) * _gdim[1] + 0;
                dp      = (double)_gdim[1] * (grid[k_out][9] - grid[k_in][9]) / (2.0 * M_PI);

                for (uint j = 0; j < _gdim[1]; j++) {
                    k = i * _gdim[1] + j;

                    if (grid[k][3] < _zc) { // k odtreni nedochazi
                        continue;
                    }

                    // urceni vnitrnich prilehajicih bunek
                    j_in                    = (double)j + dp;
                    j_in                    = (j_in > 0.0) ? j_in : (double)_gdim[1] + j_in;

                    k_lower_left            = (i + 1) * _gdim[1] + (uint)std::floor(j_in) % _gdim[1];
                    k_lower_right           = (i + 1) * _gdim[1] + (uint)std::floor(j_in + 1.0) % _gdim[1];
                    part_lower_right        = std::fmod(j_in, 1.0);
                    part_lower_left         = 1.0 - part_lower_right;
                    
                    // odtrzena hmotnost
                    w                       = getRandW();
                    mb                      = grid[k][5] * (0.8 - w);

                    // odebrat od zdrojove bunky
                    grid[k][5]              -= mb;

                    // 'vynulovat' rychlost zdrojove bunky
                    grid[k][4]              = 0.0;
                    
                    // reset vychyleni zdrojove bunky
                    grid[k][3]              = 2.0;

                    grid[k_lower_left][5]   += part_lower_left * mb; // m
                    grid[k_lower_left][6]   = part_lower_left * mb; // dm
                    grid[k_lower_right][5]  += part_lower_right* mb; // m
                    grid[k_lower_right][6]  = part_lower_right* mb; // dm

                    // zapsat do drain datasetu
                    _drain[_n][i]           += mb;
                }
            } else {
                for (uint j = 0; j < _gdim[1]; j++) {
                    k = i * _gdim[1] + j;

                    if (grid[k][3] < _zc) { // k odtreni nedochazi
                        continue;
                    }

                    // odtrzena hmotnost
                    w               = getRandW();
                    mb              = grid[k][5] * (0.8 - w);

                    // odebrat od zdrojove bunky
                    grid[k][5]      -= mb;

                    // 'vynulovat' rychlost zdrojove bunky
                    grid[k][4]      = 0.0;
                    
                    // reset vychyleni zdrojove bunky
                    grid[k][3]      = 2.0;
                    
                    // hmotu na drain
                    _drain[_n][i]   += mb;
                }
            }
        }

        // blob?
        if (_hasScheduler) {
            _scheduler->run(n);
        }
    }
};

#endif
