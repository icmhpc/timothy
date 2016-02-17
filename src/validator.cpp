//
// Created by Grzegorz Bokota on 04.02.16.
//
#include "validator/validator.h"
#include <iostream>
#include <cstdlib>


int main(int argc, char** argv){
  struct parsed_config config;
  int res;
  if (argc > 1){
    res = readFromPath(argv[1], &config);
  } else {
    res = readFromFile(stdin, &config);
  }
  if (res != (int) OK){
    fprintf(stderr, "Parser error %d in line %zu\n", res, config.number_of_sections);
    return EXIT_FAILURE;
  }
  printf("%d\n", sectionExist( (char *) "GLOBAL", &config));
  prettyPrint(stdout, &config);
  std::pair<int, std::string>validate_result = validate_config(&config);
  if (!validate_result.first){
    std::cout << validate_result.second;
  }
  deleteParsedConfig(&config);
  return EXIT_SUCCESS;
}

