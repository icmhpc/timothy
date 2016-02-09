//
// Created by czaki on 04.02.16.
//
#include <stdio.h>
#include <stdlib.h>
#include "ini_parser/ini_parser.h"

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
  printf("%d\n", sectionExist("GLOBAL", &config));
  prettyPrint(stdout, &config);
  deleteParsedConfig(&config);
  return EXIT_SUCCESS;
}

