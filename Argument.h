#ifndef ARGUMENT_H
#define ARGUMENT_H

template<class T>
class Argument {
private:
    std::string _name;
    std::string _help;
    T _value;
    bool _has_value = false;
public:
    Argument(std::string name) {
        this->setName(name);
        this->setHelp("Not specified");
    }

    Argument(std::string name, std::string help) : Argument(name) {
        this->setHelp(help);
    }

    Argument(std::string name, std::string help, T value) : Argument(name, help) {
        this->setValue(value);
    }

    void setName(std::string name) {this->_name = name;}

    std::string getName() {return this->_name;}

    void setValue(T value) {
        this->_value = value;
        this->_has_value = true;
    } 

    T getValue() {return this->_value;}

    void setHelp(std::string help) {this->_help = help;} 

    std::string getHelp() {return this->_help;}

    bool hasValue() {return this->_has_value;}
};

#endif
