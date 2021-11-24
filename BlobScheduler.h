#ifndef BLOBSCHEDULER_H
#define BLOBSCHEDULER_H

#include <iostream>
#include <fstream>
#include <sstream>

#include "rapidjson/document.h"
namespace json = rapidjson;

class BlobData{
private:
    uint        _i;
    uint        _r;
    double      _a;
    double      _m;
public:
    BlobData(uint i, uint r, double a, double m) {
        _i = i;
        _r = r;
        _a = a;
        _m = m;
    }

    uint get_i() {return _i;}
    uint get_r() {return _r;}
    uint get_a() {return _a;}
    uint get_m() {return _m;}
};

typedef std::vector<BlobData*> btype;
typedef std::map<uint, btype> stype;

class BlobScheduler {
private:
    stype                   _schedule;
    double**                _grid;
    std::vector<size_t>     _dim;
    bool                    _hasSchedule = false;
public:
    BlobScheduler(std::vector<size_t> dim) {
        _dim = dim;
    }

    BlobScheduler(std::vector<size_t> dim, double** grid) : BlobScheduler(dim) {
        setGrid(grid);
    }

    BlobScheduler(std::vector<size_t> dim, double** grid, std::string jsonSchedule) : BlobScheduler(dim, grid) {
        setSchedule(jsonSchedule);
        _hasSchedule = true;
    }

    void setSchedule(std::string jsonSchedule) {
        std::ifstream infile(jsonSchedule);
        std::string content((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
        json::Document d;
        d.Parse(content.c_str());
        const json::Value& blobs = d["blobs"];
        uint step, ic, r;
        double az, mc;
        for (json::Value::ConstValueIterator it0 = blobs.Begin(); it0 != blobs.End(); ++it0) {
            step                        = (*it0)["step"].GetInt();
            const json::Value& data     = (*it0)["data"].GetArray();
            std::vector<BlobData*> tmp;
            for (json::Value::ConstValueIterator it1 = data.Begin(); it1 != data.End(); ++it1) {
                ic  = (*it1)["i"].GetInt();
                r   = (*it1)["r"].GetInt();
                az  = (*it1)["a"].GetDouble();
                mc  = (*it1)["m"].GetDouble();
                tmp.push_back(new BlobData(ic, r, az, mc));
            }
            _schedule[step] = tmp;
        }
    }

    bool hasShedule() {
        return _hasSchedule;
    }

    void setGrid(double** grid) {
        _grid = grid;
    }

    void run (uint step) {
        stype::iterator it = _schedule.find(step);
        if (it != _schedule.end()) {
            for (auto blob : it->second) {
                addBlob(blob);
            }
        }
    }

    void addBlob(BlobData* b) {
        // nalezeni stredove bunky
        double da;
        std::vector<double> das; // pro rozdily mezi zvolenym azm. a azm. bunky
        for (uint k = b->get_i() * _dim[1]; k < b->get_i() * _dim[1] + _dim[1]; k++) {
            da = _grid[k][9] - b->get_a();
            das.push_back((da < 0.0) ? da * -1.0 : da);
        }
        uint j  = std::distance(das.begin(), std::min_element(das.begin(), das.end()));
        uint kc = b->get_i() * _dim[1] + j;

        // kontrola podminek vsech bunek
        uint ri, rc;
        double lx, ly, l;
        for (uint k = 0; k < _dim[0] * _dim[1]; k++) {
            ri  = _dim[0] - _grid[k][0];
            rc  = _dim[0] - _grid[kc][0];
            lx  = ri * cos(_grid[k][9]) - rc * cos(b->get_a());
            ly  = ri * sin(_grid[k][9]) - rc * sin(b->get_a());
            l   = sqrt(pow(lx, 2.0) + pow(ly, 2.0));
            if (l > b->get_r()) {
                continue;
            }
            _grid[k][5] += (b->get_m() / (0.8 * l));
            _grid[k][6] += (b->get_m() / (0.8 * l));
        }
    }
};

#endif
