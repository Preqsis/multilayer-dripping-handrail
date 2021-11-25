#ifndef BLOBDATA_H
#define BLOBDATA_H

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

#endif
