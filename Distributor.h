#ifndef DISTRIBUTOR_H
#define DISTRIBUTOR_H

#include <random>

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

    std::vector<uint> _probs;
    std::vector<double> RingTemp;
    std::vector<double> _RotationProfile;
public:
    Distributor(uint idim, uint jdim, uint ic, uint jc, double** drain, double dx) {
        _idim = idim;
        _jdim = jdim;
        _ic = ic;
        _jc = jc;
        _dx = dx;
        _q = 0.5;
        _zc = 5.5;
        _drain = drain;
        _n = 0;
    }

    void setRotationProfile(std::vector<double> profile) {
        _RotationProfile.clear();
        _RotationProfile = profile;
    }

    double getRandW() {
        std::random_device dev;
        std::default_random_engine engine(dev());
        std::uniform_real_distribution<double> dist(0.0, 0.3);
        return dist(engine);
    }

    void step(double** grid) {
        // rotace vnejsiho prstence o "jeden pohyb"
        // zabranuje driftovani vliven nepresnosti datovych typu
        uint k;
        for (uint j = 0; j < _jdim; j++) {
            k = j;
            RingTemp.push_back(grid[k][9]);
        }
        std::rotate(RingTemp.begin(), RingTemp.begin() + 1, RingTemp.end());
        for (uint j = 0; j < _jdim; j++) {
            k = j;
            grid[k][9] = RingTemp[j];
            grid[k][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi
        }
        RingTemp.clear();

        // rotace ostatnich podle profilu
        for (uint i = 1; i < _idim; i++) {
            for (uint j = 0; j < _jdim; j++) {
                k = i * _jdim + j;
                grid[k][9] = std::fmod(grid[k][9] + _RotationProfile[i], 2.0 * M_PI);
                grid[k][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi
            }
        }

        // nalezeni pritokove bunky + pritok
        uint idx = 0;
        for (uint j = 0; j < _jdim; j++) {
            k = j;
            if (grid[k][9] == 0.0) {
                idx = k;
                break;
            }
        }
        grid[idx][5] += _q;
        grid[idx][6] = _q; // u pritokove bunky presne odpovida _q, zbytek 0

        double mb;
        double w;
        double dp;
        uint k_out, k_in;
        uint k_lower_right, k_lower_left;
        double j_in, j_out;
        double part_lower_left, part_lower_right;
        for (uint i=_idim-1; i < _idim; i--) { // od stredu ven, unsigned preskakuje na maximalni hodnotu!!!
            if (i < _idim - 1) {
                k_out = i * _jdim + 0;
                k_in = (i + 1) * _jdim + 0;
                dp = (double)_jdim * (grid[k_out][9] - grid[k_in][9]) / (2.0 * M_PI);

                for (uint j = 0; j < _jdim; j++) {
                    k = i * _jdim + j;

                    if (grid[k][3] < _zc) { // k odtreni nedochazi
                        continue;
                    }

                    // urceni vnitrnich prilehajicih bunek
                    j_in = (double)j + dp;
                    j_in = (j_in > 0.0) ? j_in : (double)_jdim + j_in;

                    k_lower_left = (i + 1) * _jdim + (uint)std::floor(j_in) % _jdim;
                    k_lower_right = (i + 1) * _jdim + (uint)std::floor(j_in + 1.0) % _jdim;

                    part_lower_right = std::fmod(j_in, 1.0);
                    part_lower_left = 1.0 - part_lower_right;
                    
                    // odtrzena hmotnost
                    //mb = grid[k][5] - (0.2 * grid[k][5] + 0.3);
                    //mb = grid[k][5] * 0.5; // hausnumero
                    w = getRandW();
                    mb = grid[k][5] * (0.8 - w);

                    // odebrat od zdrojove bunky
                    grid[k][5] -= mb;

                    // 'vynulovat' rychlost zdrojove bunky
                    grid[k][4] = 0.0;
                    
                    // reset vychyleni zdrojove bunky
                    grid[k][3] = 2.0;

                    grid[k_lower_left][5] += part_lower_left * mb; // m
                    grid[k_lower_left][6] = part_lower_left * mb; // dm
                    grid[k_lower_right][5] += part_lower_right* mb; // m
                    grid[k_lower_right][6] = part_lower_right* mb; // dm

                    // zapsat do drain datasetu
                    _drain[_n][i] += mb;
                }
            } else {
                for (uint j = 0; j < _jdim; j++) {
                    k = i * _jdim + j;

                    if (grid[k][3] < _zc) { // k odtreni nedochazi
                        continue;
                    }

                    // odtrzena hmotnost
                    //mb = grid[k][5] - (0.2 * grid[k][5] + 0.3);
                    //mb = grid[k][5] * 0.5; // hausnumero
                    w = getRandW();
                    mb = grid[k][5] * (0.8 - w);

                    // odebrat od zdrojove bunky
                    grid[k][5] -= mb;

                    // 'vynulovat' rychlost zdrojove bunky
                    grid[k][4] = 0.0;
                    
                    // reset vychyleni zdrojove bunky
                    grid[k][3] = 2.0;
                    
                    // hmotu na drain
                    _drain[_n][i] += mb;
                }
            }
        }

        _n += 1;
        
    }
};

#endif
