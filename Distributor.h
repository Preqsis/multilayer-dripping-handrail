#ifndef DISTRIBUTOR_H
#define DISTRIBUTOR_H

#include <random>
#include "BlobScheduler.h"

class Distributor {
private:
    double _dx;
    double _q;
    double _zc;

    double*** _grid;
    double** _drain;

    std::vector<size_t> _dim;          // comms dimensions
    std::vector<size_t> _drain_dim;          // grid dimensions
    std::vector<double> _rProfile;      // rotation profile
    std::vector<double> _tProfile;      // temperature profile

    BlobScheduler* _scheduler;
    bool _hasScheduler = false;
public:
    Distributor(double*** grid, double** drain, std::vector<size_t> dim, std::vector<size_t> drain_dim) {
        _grid       = grid;
        _drain      = drain;
        _dim        = dim;
        _drain_dim  = drain_dim;

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

    void step (uint n) {
        double mb;
        double w;
        double dp;
        uint j_right, j_left;
        double j_in, j_out;
        double part_left, part_right;

        // rotace vnejsiho prstence o "jeden pohyb"
        // zabranuje driftovani vliven nepresnosti datovych typu
        std::vector<double> tmp;
        for (uint j = 0; j < _dim[1]; j++) {
            tmp.push_back(_grid[0][j][9]);
        }
        std::rotate(tmp.begin(), tmp.begin() + 1, tmp.end());
        for (uint j = 0; j < _dim[1]; j++) {
            _grid[0][j][9] = tmp[j];
            _grid[0][j][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi
        }

        // rotace ostatnich podle profilu
        for (uint i = 1; i < _dim[0]; i++) {
            for (uint j = 0; j < _dim[1]; j++) {
                _grid[i][j][9] = std::fmod(_grid[i][j][9] + _rProfile[i], 2.0 * M_PI);
                _grid[i][j][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi
            }
        }

        // nalezeni pritokove bunky + pritok
        uint idx = 0;
        for (uint j = 0; j < _dim[1]; j++) {
            if (_grid[0][j][9] == 0.0) {
                idx = j;
                break;
            }
        }
        _grid[0][idx][5] += _q;
        _grid[0][idx][6] = _q; // u pritokove bunky presne odpovida _q, zbytek 0

        for (uint i=_dim[0]-1; i < _dim[0]; i--) { // od stredu ven, unsigned preskakuje na maximalni hodnotu!!!
            if (i < _dim[0] - 1) {
                // posun uhlu mezi prstenci
                dp      = (double)_dim[1] * (_grid[i][0][9] - _grid[i+1][0][9]) / (2.0 * M_PI);

                for (uint j = 0; j < _dim[1]; j++) {
                    //k = i * _gim[1] + j;

                    if (_grid[i][j][3] < _zc) { // k odtreni nedochazi
                        _drain[i][j] = 0.0;
                        continue;
                    }


                    // urceni vnitrnich prilehajicih bunek
                    j_in        = (double)j + dp;
                    j_in        = (j_in > 0.0) ? j_in : (double)_dim[1] + j_in;

                    j_left      = (uint)std::floor(j_in) % _dim[1];
                    j_right     = (uint)std::floor(j_in + 1.0) % _dim[1];

                    part_right  = std::fmod(j_in, 1.0);
                    part_left   = 1.0 - part_right;

                    // odtrzena hmotnost
                    w                       = getRandW();
                    mb                      = _grid[i][j][5] * (0.8 - w);
                    
                    // odebrat od zdrojove bunky
                    _grid[i][j][5]          -= mb;

                    // 'vynulovat' rychlost zdrojove bunky
                    _grid[i][j][4]          = 0.0;
                    
                    // reset vychyleni zdrojove bunky
                    _grid[i][j][3]          = 2.0;

                    _grid[i+1][j_left][5] += part_left * mb; // m
                    _grid[i+1][j_left][6] = part_left * mb; // dm

                    _grid[i+1][j_right][5] += part_right* mb; // m
                    _grid[i+1][j_right][6] = part_right* mb; // dm

                    // hmotu na drain
                    _drain[i][j]            = mb;
                }
            } else {
                for (uint j = 0; j < _dim[1]; j++) {
                    //k = i * _gdim[1] + j;

                    if (_grid[i][j][3] < _zc) { // k odtreni nedochazi
                        _drain[i][j] = 0.0;
                        continue;
                    }

                    // odtrzena hmotnost
                    w               = getRandW();
                    mb              = _grid[i][j][5] * (0.8 - w);

                    // odebrat od zdrojove bunky
                    _grid[i][j][5]  -= mb;

                    // 'vynulovat' rychlost zdrojove bunky
                    _grid[i][j][4]  = 0.0;
                    
                    // reset vychyleni zdrojove bunky
                    _grid[i][j][3]      = 2.0;
                    
                    // hmotu na drain
                    _drain[i][j]            = mb;
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
