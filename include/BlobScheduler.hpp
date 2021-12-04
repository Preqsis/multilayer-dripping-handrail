#ifndef BLOBSCHEDULER_HPP
#define BLOBSCHEDULER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <math.h>

#include "rapidjson/document.h"
namespace json = rapidjson;

#include "BlobData.h"

typedef std::vector<BlobData*> btype;
typedef std::map<uint, btype> stype;

class BlobScheduler {
private:
    stype                   _schedule;
    double***               _grid;
    std::vector<size_t>     _dim;
    bool                    _hasSchedule;
public:
    BlobScheduler() {
        _hasSchedule = false;
    }

    BlobScheduler(double*** grid) : BlobScheduler() {
        this->setGrid(grid);
    }

    BlobScheduler(double*** grid, std::vector<size_t> dim) : BlobScheduler(grid) {
        _dim = dim;
    }

    BlobScheduler(double*** grid, std::vector<size_t> dim, std::string schedule_file) : BlobScheduler(grid, dim) {
        this->setSchedule(schedule_file);
        _hasSchedule = true;
    }

    void setSchedule(std::string schedule_file) {
        std::ifstream infile(schedule_file);
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

    void setGrid(double*** grid) {
        _grid = grid;
    }

    void run (uint step) {
        stype::iterator it = _schedule.find(step);
        if (it != _schedule.end()) {
            for (auto blob : it->second) {
                this->addBlob(blob);
            }
        }
    }

    void addBlob(BlobData* b) {
        // nalezeni stredove bunky
        uint ic = b->get_i();

        // kontrola podminek vsech bunek
        uint ri, rc, i, j;
        double lx, ly, l;
        for (i = 0; i < _dim[0]; i++) {
            for (j = 0; j < _dim[1]; j++) {
                ri  = _dim[0] - i;
                rc  = _dim[0] - ic;
                lx  = ri * cos(_grid[i][j][8]) - rc * cos(b->get_a());
                ly  = ri * sin(_grid[i][j][8]) - rc * sin(b->get_a());
                l   = sqrt(pow(lx, 2.0) + pow(ly, 2.0));
                if (l > b->get_r()) {
                    continue;
                }
                _grid[i][j][5] += (b->get_m() / (0.8 * l));
                _grid[i][j][6] += (b->get_m() / (0.8 * l));
            }
        }
    }
};

#endif
