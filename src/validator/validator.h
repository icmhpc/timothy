//
// Created by Grzegorz Bokota on 11.02.16.
//

#ifndef TIMOTHY_VALIDATOR_H
#define TIMOTHY_VALIDATOR_H

extern "C"{
#include "../ini_parser/ini_parser.h"
};

#include <string>

std::pair<bool,std::string> validate_config(parsed_config * c);


#endif //TIMOTHY_VALIDATOR_H
