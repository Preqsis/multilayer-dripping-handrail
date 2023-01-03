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

void Distributor::setParams(double M, double r_in, double r_out, double Q, double q, double T_flow, double psi) {
    _M      = M;
    _r_in   = r_in;
    _r_out  = r_out;
    _Q      = Q;
    _q      = q;
    _dt     = 2.0 * M_PI * std::sqrt(std::pow(_r_out, 3.0) / cs::G / _M) / _dim[1];
    // TODO: do proper qs computation
    _qs     = _Q * _dt / _q / 0.1; // divide by dx = 0.1
    _dr     = (r_out - r_in) / (double) _dim[1];
    _T_flow = T_flow;
    _psi    = psi;
}

double Distributor::get_dt() {
    return _dt;
}

double Distributor::get_qs() {
    return _qs;
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

            _grid[i][j][5] += (b->m / (0.8 * (l < 2 ? 2 : l)));
            _grid[i][j][6] += (b->m / (0.8 * (l < 2 ? 2 : l)));
            _grid[i][j][10] = getMixedTemp(_grid[i][j][10], _T_flow, _grid[i][j][5], b->m);
        }
    }
}

void Distributor::setInflux(double q) {

}

double Distributor::getRandW() {
    std::random_device dev;
    std::default_random_engine engine(dev());
    std::uniform_real_distribution<double> dist(0.1, 0.3);
    return dist(engine);
}

double Distributor::getMixedTemp(double T1, double T2, double m1, double m2) {
    return (m1 * T1 + m2 * T2) / (m1 + m2);
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

    // chalzeni absolutne cerneho telesa
    double S, n;
    for (i = 0; i < _dim[0]; i++) {
        S = 2.0 * M_PI * _grid[i][0][7] * _dr / _dim[1];  
        for (j = 0; j < _dim[1]; j++) {
            if (_grid[i][j][5] <= 0.0) {
                continue;
            }
            n = _grid[i][j][5] * _qs / cs::M_h;
            _grid[i][j][10] = 1.0 / std::pow(((2.0 * cs::sigma * S * _dt) / (n * cs::R) + 1.0 / std::pow(_grid[i][j][10], 3.0) ), 1.0/3.0);
        }
    }

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
        _grid[0][j][9] = 0.0; 
    }

    // rotace ostatnich podle profilu
    for (i = 1; i < _dim[0]; i++) {
        for (j = 0; j < _dim[1]; j++) {
            _grid[i][j][8] = std::fmod(_grid[i][j][8] + _rProfile[i], 2.0 * M_PI);
            _grid[i][j][6] = 0; // dm u vsech vynulovat, uvazujeme pohyb casti nikoliv tok mezi nimi

            _grid[i][j][9] = 0.0; 
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
    
    // teplota v pritokove bunce
    _grid[0][idx][10] = (_grid[0][idx][5] * _grid[0][idx][10] + _q * _T_flow) / (_grid[0][idx][5] + _q);

    // hmota v pritokove bunce
    _grid[0][idx][5] += _q;
    _grid[0][idx][6] = _q; // u pritokove bunky presne odpovida _q, zbytek 0

    for (i = _dim[0]-1; i < _dim[0]; i--) { // od stredu ven, unsigned preskakuje na maximalni hodnotu!!!
        if (i < _dim[0] - 1) {
            // posun uhlu mezi prstenci
            dp      = (double)_dim[1] * (_grid[i][0][8] - _grid[i+1][0][8]) / (2.0 * M_PI);

            for (j = 0; j < _dim[1]; j++) {
                if (_grid[i][j][3] < _zc) { // k odtreni nedochazi
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
                // mb                      = _grid[i][j][5] * (0.8 - w);
                mb                      = _grid[i][j][5] * _psi;
                if (mb <= 1e-8) { // jaky limit presne? !!!!!!!!!
                    mb                  = _grid[i][j][5];
                    _grid[i][j][5]      = 0.0;
                    _grid[i][j][10]     = 0.0;
                } else {
                    _grid[i][j][5]      -= mb; // m
                }
                
                // odebrat od zdrojove bunky

                // DORESIT !!! --> co se zapornym dm???
                //_grid[i][j][6]          -= mb; // dm 

                // 'vynulovat' rychlost zdrojove bunky
                _grid[i][j][4]          = 0.0;
                
                // reset vychyleni zdrojove bunky
                _grid[i][j][3]          = 2.0;
                
                // temperature change by mixing
                //_grid[i+1][j_left][10] = (_grid[i+1][j_left][5] * _grid[i+1][j_left][10] + part_left * mb * _grid[i][j][10]) / (_grid[i+1][j_left][5] + part_left * mb);
                //_grid[i+1][j_right][10] = (_grid[i+1][j_right][5] * _grid[i+1][j_right][10] + part_right * mb * _grid[i][j][10]) / (_grid[i+1][j_right][5] + part_right * mb);

                _grid[i+1][j_left][10] = getMixedTemp(_grid[i][j][10], _grid[i+1][j_left][10], part_left * mb, _grid[i+1][j_left][5]);
                _grid[i+1][j_right][10] = getMixedTemp(_grid[i][j][10], _grid[i+1][j_right][10], part_right * mb, _grid[i+1][j_right][5]);


                // total cell mass
                _grid[i+1][j_left][5] += part_left * mb; // m
                _grid[i+1][j_right][5] += part_right * mb; // m

                // cell mass change
                _grid[i+1][j_left][6] += part_left * mb; // dm
                _grid[i+1][j_right][6] += part_right * mb; // dm

                // hmotu na drain
                _grid[i][j][9] = mb;
                //_grid[i+1][j_left][9] += part_left * mb; // drain
                //_grid[i+1][j_right][9] += part_right * mb; // drain
            }
        } else {
            for (j = 0; j < _dim[1]; j++) {
                if (_grid[i][j][3] < _zc) { // k odtreni nedochazi
                    continue;
                }

                // odtrzena hmotnost
                w               = getRandW();
                // mb              = _grid[i][j][5] * (0.8 - w);
                mb              = _grid[i][j][5] * _psi;
                if (mb <= 1e-8) {
                    mb                  = _grid[i][j][5];
                    _grid[i][j][5]      = 0.0;
                    _grid[i][j][10]     = 0.0;
                } else {
                    _grid[i][j][5]      -= mb; // m
                }

                // odebrat od zdrojove bunky
                //_grid[i][j][5]  -= mb;
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

    // free-free emission heating
    double A = 1.0; // heating efficiency
    double K = A * cs::G * _M * _dr * cs::M_h / (3.0 * cs::R); // constant part of the equation
    for (i = 0; i < _dim[0]; i++) {
        for (j = 0; j < _dim[1]; j++) {
            if (_grid[i][j][5] > 0.0) {
                _grid[i][j][10] +=  K * _grid[i][j][9] / (std::pow(_grid[i][j][7], 2.0) * _grid[i][j][5]);
            }
        }
    }
}
