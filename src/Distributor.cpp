#include <random>
#include "BlobScheduler.hpp"
#include "Distributor.hpp"

Distributor::Distributor() {
    _q              = 0.5;
    _zc             = 5.5;
    _hasScheduler   = false;
}

Distributor::Distributor(double*** grid, std::vector<size_t> dim) : Distributor() {
    _grid           = grid;
    _dim            = dim;
}

Distributor::Distributor(double*** grid, std::vector<size_t> dim, double q) :  Distributor(grid, dim) {
    setInflux(q);
}

Distributor::Distributor(double*** grid, std::vector<size_t> dim, double q, BlobScheduler* scheduler) : Distributor(grid, dim, q) {
    setBlobScheduler(scheduler);
}

Distributor::~Distributor() {
    if (_hasScheduler)
        delete _scheduler;
}

void Distributor::setBlobScheduler(BlobScheduler* scheduler) {
    _scheduler      = scheduler;
    _hasScheduler   = true;
}

void Distributor::setRotationProfile(std::vector<double> profile) {
    _rProfile.clear();
    _rProfile = profile;
}

void Distributor::setTemperatureProfile(std::vector<double> profile) {
    _tProfile.clear();
    _tProfile = profile;
}

void Distributor::setInflux(double q) {
    _q = q;
}

double Distributor::getRandW() {
    std::random_device dev;
    std::default_random_engine engine(dev());
    std::uniform_real_distribution<double> dist(0.0, 0.3);
    return dist(engine);
}

void Distributor::run (uint s) {
    double mb, w, dp, j_in, part_left, part_right;
    uint i, j, j_right, j_left, idx;

    // rotace vnejsiho prstence o "jeden pohyb"
    // zabranuje driftovani vliven nepresnosti datovych typu
    std::vector<double> tmp;
    for (j = 0; j < _dim[1]; j++) {
        tmp.push_back(_grid[0][j][8]);
    }
    std::rotate(tmp.begin(), tmp.begin() + 1, tmp.end());
    for (j = 0; j < _dim[1]; j++) {
        _grid[0][j][8] = tmp[j];
        _grid[0][j][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi
    }

    // rotace ostatnich podle profilu
    for (i = 1; i < _dim[0]; i++) {
        for (j = 0; j < _dim[1]; j++) {
            _grid[i][j][8] = std::fmod(_grid[i][j][8] + _rProfile[i], 2.0 * M_PI);
            _grid[i][j][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi
        }
    }

    // nalezeni pritokove bunky + pritok
    idx = 0;
    for (j = 0; j < _dim[1]; j++) {
        if (_grid[0][j][8] == 0.0) {
            idx = j;
            break;
        }
    }
    _grid[0][idx][5] += _q;
    _grid[0][idx][6] = _q; // u pritokove bunky presne odpovida _q, zbytek 0

    for (i = _dim[0]-1; i < _dim[0]; i--) { // od stredu ven, unsigned preskakuje na maximalni hodnotu!!!
        if (i < _dim[0] - 1) {
            // posun uhlu mezi prstenci
            dp      = (double)_dim[1] * (_grid[i][0][8] - _grid[i+1][0][8]) / (2.0 * M_PI);

            for (j = 0; j < _dim[1]; j++) {
                if (_grid[i][j][3] < _zc) { // k odtreni nedochazi
                    _grid[i][j][9] = 0.0; 
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
                _grid[i][j][5]          -= mb; // m
                _grid[i][j][6]          -= mb; // dm

                // 'vynulovat' rychlost zdrojove bunky
                _grid[i][j][4]          = 0.0;
                
                // reset vychyleni zdrojove bunky
                _grid[i][j][3]          = 2.0;

                _grid[i+1][j_left][5] += part_left * mb; // m
                _grid[i+1][j_left][6] = part_left * mb; // dm

                _grid[i+1][j_right][5] += part_right* mb; // m
                _grid[i+1][j_right][6] = part_right* mb; // dm

                // hmotu na drain
                _grid[i][j][9] = mb;
            }
        } else {
            for (j = 0; j < _dim[1]; j++) {
                if (_grid[i][j][3] < _zc) { // k odtreni nedochazi
                    _grid[i][j][9] = 0.0; 
                    continue;
                }

                // odtrzena hmotnost
                w               = getRandW();
                mb              = _grid[i][j][5] * (0.8 - w);

                // odebrat od zdrojove bunky
                _grid[i][j][5]  -= mb;
                _grid[i][j][6]  -= mb;

                // 'vynulovat' rychlost zdrojove bunky
                _grid[i][j][4]  = 0.0;
                
                // reset vychyleni zdrojove bunky
                _grid[i][j][3]  = 2.0;
                
                // hmotu na drain
                _grid[i][j][9] = mb;
            }
        }
    }

    // blob?
    if (_hasScheduler) {
        _scheduler->run(s);
    }
}
