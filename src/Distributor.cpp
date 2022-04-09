#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <math.h>

#include "rapidjson/document.h"
namespace json = rapidjson;

#include "Constants.hpp"
namespace cs = Constants;

#include "Distributor.hpp"

Distributor::Distributor() {
    _zc         = 5.5;
    _hasBlobs   = false;
}

Distributor::Distributor(double*** grid, std::vector<size_t> dim) : Distributor() {
    _grid           = grid;
    _dim            = dim;
}

Distributor::Distributor(double*** grid, std::vector<size_t> dim, std::string blob_file) : Distributor(grid, dim) {
    setBlobSchedule(blob_file);
    _hasBlobs   = true;
}

Distributor::~Distributor() {}

void Distributor::setRotationProfile(std::vector<double> profile) {
    _rProfile.clear();
    _rProfile = profile;
}

void Distributor::setParams(double M, double r_in, double r_out, double Q, double q) {
    _M      = M;
    _r_in   = r_in;
    _r_out  = r_out;
    _Q      = Q;
    _q      = q;
    _dt     = 2.0 * M_PI * std::sqrt(std::pow(_r_out, 3.0) / cs::G / _M / cs::m_sun);
    _qs     = _Q * _dt / _q;
    _T_in   = std::pow((3.0 * cs::G * cs::m_sun * _M * _Q) / (8.0 * M_PI * cs::sigma * std::pow(_r_in, 3.0)), 0.25);
}

double Distributor::get_dt() {
    return _dt;
}

void Distributor::setBlobSchedule(std::string blob_file) {
    std::ifstream infile(blob_file);
    std::string content((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
    json::Document d;
    d.Parse(content.c_str());
    const json::Value& blobs = d["blobs"];
    uint step;
    BlobData* b;
    for (json::Value::ConstValueIterator it0 = blobs.Begin(); it0 != blobs.End(); ++it0) {
        step                        = (*it0)["step"].GetInt();
        const json::Value& data     = (*it0)["data"].GetArray();
        std::vector<BlobData*> tmp;
        for (json::Value::ConstValueIterator it1 = data.Begin(); it1 != data.End(); ++it1) {
            b = new BlobData();
            b->i = (*it1)["i"].GetInt();
            b->r = (*it1)["r"].GetInt();
            b->a = (*it1)["a"].GetDouble();
            b->m = (*it1)["m"].GetDouble();
            tmp.push_back(b);
        }
        _blobSchedule[step] = tmp;
    }
}

void Distributor::addBlob(BlobData* b) {
    // kontrola podminek vsech bunek
    uint ri, rc, i, j;
    double lx, ly, l;
    for (i = 0; i < _dim[0]; i++) {
        for (j = 0; j < _dim[1]; j++) {
            ri  = _dim[0] - i;
            rc  = _dim[0] - b->i;
            lx  = ri * cos(_grid[i][j][8]) - rc * cos(b->a);
            ly  = ri * sin(_grid[i][j][8]) - rc * sin(b->a);
            l   = sqrt(pow(lx, 2.0) + pow(ly, 2.0));
            if (l > b->r) {
                continue;
            }
            _grid[i][j][5] += (b->m / (0.8 * l));
            _grid[i][j][6] += (b->m / (0.8 * l));
        }
    }
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

void Distributor::runBlobs(uint step) {
    stype::iterator it = _blobSchedule.find(step);
    if (it != _blobSchedule.end()) {
        for (auto blob : it->second) {
            this->addBlob(blob);
        }
    }
}

void Distributor::run(uint s) {
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

                // DORESIT !!! --> co se zapornym dm???
                //_grid[i][j][6]          -= mb; // dm 

                // 'vynulovat' rychlost zdrojove bunky
                _grid[i][j][4]          = 0.0;
                
                // reset vychyleni zdrojove bunky
                _grid[i][j][3]          = 2.0;

                _grid[i+1][j_left][5] += part_left * mb; // m
                _grid[i+1][j_left][6] = part_left * mb; // dm

                _grid[i+1][j_right][5] += part_right * mb; // m
                _grid[i+1][j_right][6] = part_right * mb; // dm

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
                //_grid[i][j][6]  -= mb;

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
    if (_hasBlobs) {
        runBlobs(s);
    }

    // temperature
    /*
    double f = 1.0; // inner boundary effect??
    double ef = 1.0; // energy conversion efficiency aka. A'
    double A = ef * cs::G * _M * cs::m_sun / (8.0 * M_PI *std::pow(_r_in, 3.0) * cs::sigma * std::pow(_T_in, 4.0) * _dt);
    double temp;
    for (i = 0; i < _dim[0]; i++) {
        for (j = 0; j < _dim[1]; j++) {
            
            //temp = _T_in * std::pow((1.0 + A * _grid[i][j][6] * _qs) / std::pow(_grid[i][j][7] / _r_in, 3.0), 0.25) * std::pow(f, 0.25);
            //if (!isnan(temp) && _grid[i][j][6] < 0.0) {
            //    std::cout << i << ", " << j << ", " << _grid[i][j][6] << std::endl;
            //}
            
            _grid[i][j][10] = _T_in * std::pow((1.0 + A * _grid[i][j][6] * _qs) / std::pow(_grid[i][j][7] / _r_in, 3.0), 0.25) * std::pow(f, 0.25);

            //if (_grid[i][j][6] > 0.0) {
            //    std::cout << " ---->> ";
            //}
            //std::cout << s << ", " << i << ", " << j << ", " << _grid[i][j][10] << ", " << _T_in << std::endl;
        }
    }
    */
}
