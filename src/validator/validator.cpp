//
// Created by Grzegorz Bokota on 11.02.16.
//

#include <python2.7/Python.h>
#include "validator.h"
#include <sstream>
#include <map>

static const std::string error = "[Error] ";
static const std::string warning = "[Warning] ";


bool validateGlobal(parsed_config * c, std::iostream * s);

std::pair<bool,std::string> validate_config(parsed_config * c){
  std::stringstream errors;
  bool res = true;
  if (!sectionExist("Global", c)){
    errors << "No GLOBAL section in config file" << std::endl;
    res = false;
  } else {
    res &= validateGlobal(c, &errors);

  }

  return std::make_pair<int,std::string>(res, errors.str());
};



bool validateGlobal(parsed_config * c, std::iostream * s){
  bool res = true;
  std::string section_name = "GLOBAL";
  std::string prefix = " section do not contain ";
  std::map<std::string, int (*) (const char *, const char *, parsed_config *) > to_check = {
          {"number_of_steps", isNumberField}, {"size_x", isNumberField}, {"size_y", isNumberField},
          {"size_z", isNumberField}
  };
  for(auto it : to_check){
    if (fieldExist(section_name.c_str(), it.first.c_str(), c)){
      res &= it.second(section_name.c_str(), it.first.c_str(), c)
    } else {
      res = false;
      *s << error << section_name << prefix << "\"" << it.first << "\"" << std::endl;

    }
  }
  if (fieldExist(section_name.c_str(), "python_file", c)){
    if (isStringField(section_name.c_str(), "python_file", c)){
      //TODO dopisać wczytanie pliku pythonowego
    } else {
      res = false;
      *s << error << section_name << " python_file have to be string field" << std::endl;
    }

  }
  return res;
}