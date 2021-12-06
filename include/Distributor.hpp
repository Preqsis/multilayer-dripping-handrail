#ifndef DISTRIBUTOR_HPP
#define DISTRIBUTOR_HPP

struct BlobData {
    uint i;
    uint r;
    double a;
    double m;
};

typedef std::vector<BlobData*> btype;
typedef std::map<uint, btype> stype;

class Distributor {
private:
    double _q;
    double _zc;
    double*** _grid;
    std::vector<size_t> _dim;           // comms dimensions
    std::vector<double> _rProfile;      // rotation profile
    std::vector<double> _tProfile;      // temperature profile
    stype _blobSchedule;
    bool _hasBlobs;
public:
    Distributor();

    Distributor(double*** grid, std::vector<size_t> dim);

    Distributor(double*** grid, std::vector<size_t> dim, double q);

    Distributor(double*** grid, std::vector<size_t> dim, double q, std::string blob_file);

    ~Distributor();

    void setRotationProfile(std::vector<double> profile);

    void setTemperatureProfile(std::vector<double> profile);

    void setBlobSchedule(std::string blob_file);

    void addBlob(BlobData* b);

    void setInflux(double q);

    double getRandW();

    void runBlobs(uint steps);

    void run (uint s);
};

#endif
