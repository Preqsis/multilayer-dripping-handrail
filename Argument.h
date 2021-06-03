#ifndef ARGUMENT_H
#define ARGUMENT_H

template<class T>
class Argument {
private:
    std::string _name;
    std::string _help;
    T _value;
public:
    Argument(std::string name) {
        this->_name = name;
    }

    Argument(std::string name, T value) : Argument(name) {
        this->setValue(value);
    }

    Argument(std::string name, T value, std::string help) : Argument(name, value) {
        this->_help = help;
    }

    std::string getName() {return this->_name;}

    T getValue() {return this->_value;}

    void setValue(T value) {this->_value = value;} 

    std::string getHelp() {return this->_help;}
};

#endif
