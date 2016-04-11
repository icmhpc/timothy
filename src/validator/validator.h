//
// Created by Grzegorz Bokota on 11.02.16.
//

#ifndef TIMOTHY_VALIDATOR_H
#define TIMOTHY_VALIDATOR_H



#include <string>
#include "../ini_manipulator/ini_manipulator.h"

namespace timothy {
    namespace validator{
        std::pair<bool,std::string> validate_config(timothy::ini_manipulator::iniConfiguration & c);
    }
}



#endif //TIMOTHY_VALIDATOR_H
