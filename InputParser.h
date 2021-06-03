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
    std::map<std::string, Argument<lint>*> _ints;
    std::map<std::string, Argument<double>*> _doubles;
    std::map<std::string, Argument<std::string>*> _strings;
public:
    InputParser(
        std::initializer_list<Argument<lint>*> ints, 
        std::initializer_list<Argument<double>*> doubles, 
        std::initializer_list<Argument<std::string>*> strings
    ) {
        for (auto a: ints) {this->addInt(a);}
        for (auto a: doubles) {this->addDouble(a);}
        for (auto a: strings) {this->addString(a);}
    }

    InputParser(
        std::initializer_list<Argument<lint>*> ints, 
        std::initializer_list<Argument<double>*> doubles, 
        std::initializer_list<Argument<std::string>*> strings,
        int &argc,
        char **argv
    ) : InputParser(ints, doubles, strings) {
        this->run(argc, argv);
    }

    void addInt(Argument<lint>* a) {this->_ints[a->getName()] = a;}
    void addDouble(Argument<double>* a) {this->_doubles[a->getName()] = a;}
    void addString(Argument<std::string>* a) {this->_strings[a->getName()] = a;}

    bool inInts(std::string key) {
        std::map<std::string, Argument<lint>*>::iterator it = this->_ints.find(key);
        return it != this->_ints.end();
    }

    bool inDoubles(std::string key) {
        std::map<std::string, Argument<double>*>::iterator it = this->_doubles.find(key);
        return it != this->_doubles.end();
    }

    bool inStrings(std::string key) {
        std::map<std::string, Argument<std::string>*>::iterator it = this->_strings.find(key);
        return it != this->_strings.end();
    }

    void setArgument(std::string key, std::string value) {
        // zbavit pocatecnich zbytecnych pomlcek
        while (key[0] == '-') {
            key.erase(0, 1);
        }

        if (this->inInts(key)) {
            std::stringstream ss(value);
            lint x = 0;
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

    bool isKey(std::string s) {
        return s.find("--") != std::string::npos || s.find("-") != std::string::npos;    
    }

    bool isValue(std::string s)  {
        return !this->isKey(s);
    }

    bool hasValue(std::string s) {
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

    friend std::ostream& operator<<(std::ostream& os, const InputParser& p) {
        os << "arguments: " << std::endl;
        os << std::left << std::setw(24) << "  -h, --help" << "Show this help msg. and exit." << std::endl;

        for (auto const& [key, a] : p._ints) {
            os << "  --" << std::left << std::setw(20)  << a->getName()  << a->getHelp() << std::endl;
        }

        for (auto const& [key, a] : p._doubles) {
            os << "  --" << std::left << std::setw(20)  << a->getName()  << a->getHelp() << std::endl;
        }

        for (auto const& [key, a] : p._strings) {
            os << "  --" << std::left << std::setw(20)  << a->getName()  << a->getHelp() << std::endl;
        }

        return os;
    } 
};

#endif
