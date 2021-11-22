#ifndef ARGUMENTPARSER_H
#define ARGUMENTPARSER_H

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "Argument.h"

class ArgumentParser {
private:
    std::map<std::string, Argument<bool>*> _bools;
    std::map<std::string, Argument<int>*> _ints;
    std::map<std::string, Argument<double>*> _doubles;
    std::map<std::string, Argument<std::string>*> _strings;
public:
    ArgumentParser() {}

    /*
    ~ArgumentParser() {
        // delete bools
        for(std::map<std::string, Argument<bool>*>::iterator it = _bools.begin(); it != _bools.end(); it++) {
            delete it->second;
            _bools.erase(it);
        }
        // delete bools
        for(std::map<std::string, Argument<int>*>::iterator it = _ints.begin(); it != _ints.end(); it++) {
            delete it->second;
            _ints.erase(it);
        }
        // delete bools
        for(std::map<std::string, Argument<double>*>::iterator it = _doubles.begin(); it != _doubles.end(); it++) {
            delete it->second;
            _doubles.erase(it);
        }
        // delete bools
        for(std::map<std::string, Argument<std::string>*>::iterator it = _strings.begin(); it != _strings.end(); it++) {
            delete it->second;
            _strings.erase(it);
        }
    }*/

    void addArgument(Argument<bool>* arg) {
        addBool(arg);
    }

    void addArgument(Argument<int>* arg) {
        addInt(arg);
    }

    void addArgument(Argument<double>* arg) {
        addDouble(arg);
    }

    void addArgument(Argument<std::string>* arg) {
        addString(arg);
    }

    void addBool(Argument<bool>* arg) {
        _bools[arg->getName()] = arg;
    }

    void addBool(std::initializer_list<Argument<bool>*> arg) {
        for (auto a: arg) 
            addBool(a);
    }

    void addInt(Argument<int>* arg) {
        _ints[arg->getName()] = arg;
    }
    
    void addInt(std::initializer_list<Argument<int>*> arg) {
        for (auto a: arg) 
            addInt(a);
    }

    void addDouble(Argument<double>* arg) {
        _doubles[arg->getName()] = arg;
    }

    void addDouble(std::initializer_list<Argument<double>*> arg) {
        for (auto a: arg) 
            addDouble(a);
    }

    void addString(Argument<std::string>* arg) {
        _strings[arg->getName()] = arg;
    }

    void addString(std::initializer_list<Argument<std::string>*> arg) {
        for (auto a: arg) 
            addString(a);
    }

    bool inBools(const std::string& key) {
        std::map<std::string, Argument<bool>*>::iterator it = this->_bools.find(key);
        return it != this->_bools.end();
    }

    bool inInts(const std::string& key) {
        std::map<std::string, Argument<int>*>::iterator it = this->_ints.find(key);
        return it != this->_ints.end();
    }

    bool inDoubles(const std::string& key) {
        std::map<std::string, Argument<double>*>::iterator it = this->_doubles.find(key);
        return it != this->_doubles.end();
    }

    bool inStrings(const std::string& key) {
        std::map<std::string, Argument<std::string>*>::iterator it = this->_strings.find(key);
        return it != this->_strings.end();
    }

    void setArgument(std::string key, std::string value) {
        key = stripKey(key);
        if (this->inBools(key)) {
            bool x = (value == "1" || value == "true" || value == "True" || value == "TRUE");
            this->_bools[key]->setValue(x);
        } else if (this->inInts(key)) {
            std::stringstream ss(value);
            int x = 0;
            ss >> x;
            this->_ints[key]->setValue(x);
        } else if (this->inDoubles(key)) {
            std::stringstream ss(value);
            double x = 0.0;
            ss >> x;
            this->_doubles[key]->setValue(x);
        } else if (this->inStrings(key)) {
            this->_strings[key]->setValue(value);
        }
    }

    void setArgument(std::string key, bool value) {
        this->_bools[stripKey(key)]->setValue(value);
    }

    bool isKey(const std::string& s) {
        return s.find("--") != std::string::npos || s.find("-") != std::string::npos;    
    }

    std::string stripKey(std::string key) {
        while (key[0] == '-')
            key.erase(0, 1);
        return key;
    }

    bool isValue(const std::string& s)  {
        return !this->isKey(s);
    }

    bool hasValue(const std::string& s) {
        return s.find("=") != std::string::npos;
    } 

    bool parse(int &argc, char **argv) {
        std::string lastKey = "";
        bool ex             = false;
        for (uint i=1; i<argc; ++i) {
            std::string tmp = std::string(argv[i]);

            if (this->isKey(tmp)) { // je argument
                if (this->hasValue(tmp)) { // prirazene pres =
                    ex = false;
                    std::string key = tmp.substr(0, tmp.find("="));
                    tmp.erase(0, tmp.find("=") + 1);
                    this->setArgument(key, tmp);
                } else if (inBools(stripKey(tmp))) {
                    ex = false;
                    this->setArgument(tmp, true);
                } else { // pristi je hodnota tohohle
                    lastKey = tmp; 
                    ex      = true; 

                    // napoveda?
                    if (lastKey == "--help" || lastKey == "-h" || lastKey == "-help") {
                        return false;
                    }
                } 
            } else { // je hodnota
                if (!ex) continue; // osirela hodnota
                this->setArgument(lastKey, tmp);
                ex = false;
            }
        }

        return true;
    }

    bool isSet(const std::string& key) {
        if (inBools(key))
            return _bools[key]->hasValue();
        else if (inInts(key))
            return _ints[key]->hasValue();
        else if (inDoubles(key))
            return _doubles[key]->hasValue();
        else if (inStrings(key))
            return _strings[key]->hasValue();
        return false;
    }

    friend std::ostream& operator<<(std::ostream& os, const ArgumentParser& p) {
        os << "arguments: " << std::endl;
        os << std::left << std::setw(24) << "  -h, --help" << "Show this help msg. and exit." << std::endl;

        for (auto const& [key, a] : p._bools) 
            os << "  --" << std::left << std::setw(20)<< a->getName() << a->getHelp() << std::endl;

        for (auto const& [key, a] : p._ints) 
            os << "  --" << std::left << std::setw(20) << a->getName() << a->getHelp() << std::endl;

        for (auto const& [key, a] : p._doubles) 
            os << "  --" << std::left << std::setw(20) << a->getName() << a->getHelp() << std::endl;

        for (auto const& [key, a] : p._strings)
            os << "  --" << std::left << std::setw(20) << a->getName() << a->getHelp() << std::endl;

        return os;
    } 
    
    bool b(const std::string& key) {
        return _bools[key]->getValue();
    }

    int i(const std::string& key) {
        return _ints[key]->getValue();
    }

    double d(const std::string& key) {
        return _doubles[key]->getValue();
    }

    std::string s(const std::string& key) {
        return _strings[key]->getValue();
    }
};

#endif
