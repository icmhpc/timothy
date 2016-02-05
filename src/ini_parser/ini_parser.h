//
// Created by Grzegorz Bokota CENT on 02.02.16.
//

#ifndef TIMOTHY_INI_PARSER_H
#define TIMOTHY_INI_PARSER_H
#include<stdio.h>

enum parse_error {
    OK,
    NO_FILE,
    ERROR_IN_SECTION_NAME,
    TO_LONG_SECTION_NAME,
    TO_LONG_FIELD_NAME,
    MEMORY_ERROR,
    BAD_FILE_FORMAT,
    TO_LONG_FIELD_VALUE
};

enum ini_fields_type{
    NUMBER_FIELD,
    FLOAT_FIELD,
    BOOLEAN_FIELD,
    STRING_FIELD
};

union ini_data{
    int i;
    float f;
    double d;
    char * str;
};

struct fields_of_ini{
    char * name;
    union ini_data data;
    enum ini_fields_type type;
};

struct section_of_ini{
    char * name;
    size_t number_of_fields;
    struct fields_of_ini * fields;
};

struct parsed_config{
    size_t number_of_sections;
    struct section_of_ini * sections_list;
};

int readFromFile(FILE * f, struct parsed_config * res);

int readFromPath(char * path, struct parsed_config * res);

void deleteParsedConfig(struct parsed_config * c);

char ** getSectionsNamesByPrefix(char * prefix, struct parsed_config *c);

int sectionExist(char * name, struct parsed_config * c);
int fieldExist(char * section_name, char * field_name,  struct parsed_config * c);

int isNumberField(char * section_name, char * field_name,  struct parsed_config * c);
int isFloatField(char * section_name, char * field_name,  struct parsed_config * c);
int isBooleanField(char * section_name, char * field_name,  struct parsed_config * c);
int isStringField(char * section_name, char * field_name,  struct parsed_config * c);

void prettyPrint(FILE *f, struct parsed_config *c);

#endif //TIMOTHY_INI_PARSER_H

