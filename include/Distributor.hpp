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
    double _M;
    double _r_in;
    double _r_out;
    double _Q;
    double _dt;
    double _T_in;
        
    double _q;
    double _zc;

    double _qs;                         // "inner / outer" mass scaling

    double*** _grid;
    std::vector<size_t> _dim;           // comms dimensions
    std::vector<double> _rProfile;      // rotation profile
    stype _blobSchedule;
    bool _hasBlobs;
public:
    Distributor();

    Distributor(double*** grid, std::vector<size_t> dim);

    Distributor(double*** grid, std::vector<size_t> dim, std::string blob_file);

    ~Distributor();

    void setRotationProfile(std::vector<double> profile);

    void setParams(double M, double r_in, double r_out, double Q, double q);

    double get_dt();

    void setBlobSchedule(std::string blob_file);

    void addBlob(BlobData* b);

    void setInflux(double q);

    double getRandW();

    void runBlobs(uint steps);

    void run (uint s);
};

#endif
