//
// Created by czaki on 04.02.16.
//
#include <stdio.h>
#include <stdlib.h>
#include "ini_parser/ini_parser.h"

int main(int argc, char** argv){
  struct parsed_config config;
  if (argc > 1){
    readFromPath(argv[1], &config);
  } else {
    readFromFile(stdin, &config);
  }
  printf("%d\n", sectionExist("GLOBAL", &config));
  prettyPrint(stdout, &config);
  deleteParsedConfig(&config);
  return EXIT_SUCCESS;
}

