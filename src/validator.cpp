//
// Created by Grzegorz Bokota on 04.02.16.
//
#include "validator/validator.h"
#include "ini_manipulator/ini_manipulator.h"
#include <iostream>
#include <cstdlib>
#include <fstream>

namespace mp = timothy::ini_manipulator;

int main(int argc, char** argv){
  timothy::ini_manipulator::iniConfiguration c;
  try {
    if (argc > 1) {
      std::fstream file;
      file.open(argv[1]);
      c = mp::iniConfiguration(file);
      file.close();
    } else {
      c = mp::iniConfiguration(std::cin);
    }
  } catch (mp::fileFormatError e) {
    std::cerr << "[PARSING ERROR]" << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << c << std::endl;
  std::pair<int, std::string>validate_result = timothy::validator::validate_config(c);
  if (!validate_result.first){
    std::cout << validate_result.second;
  }
  return EXIT_SUCCESS;
}

