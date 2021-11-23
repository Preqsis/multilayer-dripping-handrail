#ifndef BLOBSCHEDULER_H
#define BLOBSCHEDULER_H

#include <iostream>
#include <fstream>
#include <sstream>

#include "rapidjson/document.h"
namespace json = rapidjson;

typedef std::map<uint, double> btype;
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
        uint step, i, j, k;
        double m;
        for (json::Value::ConstValueIterator it0 = blobs.Begin(); it0 != blobs.End(); ++it0) {
            step                        = (*it0)["step"].GetInt();
            const json::Value& data     = (*it0)["data"].GetArray();
            std::map<uint, double> tmp;
            for (json::Value::ConstValueIterator it1 = data.Begin(); it1 != data.End(); ++it1) {
                tmp[(*it1)["i"].GetInt() * _dim[1] + (*it1)["j"].GetInt()] = (*it1)["m"].GetDouble();
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
            for (btype::iterator bt = it->second.begin(); bt != it->second.end(); ++bt) {
                _grid[bt->first][5] += bt->second; // m
                _grid[bt->first][6] += bt->second; // dm

                std::cout << step << ", " << it->first << ", " << bt->first << ", " << bt->second << std::endl;
            }
        }
    }
};

#endif
