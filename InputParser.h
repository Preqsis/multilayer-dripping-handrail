#ifndef INPUTPARSER_H
#define INPUTPARSER_H

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "Argument.h"

typedef long unsigned int lint;

class InputParser {
private:
    std::map<std::string, Argument<bool>*> _bools;
    std::map<std::string, Argument<int>*> _ints;
    std::map<std::string, Argument<double>*> _doubles;
    std::map<std::string, Argument<std::string>*> _strings;
public:
    InputParser(
        std::initializer_list<Argument<bool>*> bools, 
        std::initializer_list<Argument<int>*> ints, 
        std::initializer_list<Argument<double>*> doubles, 
        std::initializer_list<Argument<std::string>*> strings
    ) {
        for (auto a: bools) {this->addBool(a);}
        for (auto a: ints) {this->addInt(a);}
        for (auto a: doubles) {this->addDouble(a);}
        for (auto a: strings) {this->addString(a);}
    }

    InputParser(
        std::initializer_list<Argument<bool>*> bools, 
        std::initializer_list<Argument<int>*> ints, 
        std::initializer_list<Argument<double>*> doubles, 
        std::initializer_list<Argument<std::string>*> strings,
        int &argc,
        char **argv
    ) : InputParser(bools, ints, doubles, strings) {
        this->run(argc, argv);
    }

    void addBool(Argument<bool>* a) {this->_bools[a->getName()] = a;}
    void addInt(Argument<int>* a) {this->_ints[a->getName()] = a;}
    void addDouble(Argument<double>* a) {this->_doubles[a->getName()] = a;}
    void addString(Argument<std::string>* a) {this->_strings[a->getName()] = a;}

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
        // zbavit pocatecnich zbytecnych pomlcek
        while (key[0] == '-') {
            key.erase(0, 1);
        }

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

    bool isKey(const std::string& s) {
        return s.find("--") != std::string::npos || s.find("-") != std::string::npos;    
    }

    bool isValue(const std::string& s)  {
        return !this->isKey(s);
    }

    bool hasValue(const std::string& s) {
        return s.find("=") != std::string::npos;
    } 

    bool run(int &argc, char **argv) {
        int i=1;
        std::string lastKey = "";
        bool ex = false;

        for (i; i < argc; ++i) {
            std::string tmp = std::string(argv[i]);

            if (this->isKey(tmp)) { // je argument
                if (this->hasValue(tmp)) { // prirazene pres =
                    ex = false;
                    std::string key = tmp.substr(0, tmp.find("="));
                    tmp.erase(0, tmp.find("=") + 1);
                    this->setArgument(key, tmp);
                } else { // pristi je hodnota tohohle
                    lastKey = tmp; 
                    ex = true; 

                    // napoveda
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

    friend std::ostream& operator<<(std::ostream& os, const InputParser& p) {
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
