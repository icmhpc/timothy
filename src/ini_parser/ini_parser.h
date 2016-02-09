//
// Created by Grzegorz Bokota CENT on 02.02.16.
//

#ifndef TIMOTHY_INI_PARSER_H
#define TIMOTHY_INI_PARSER_H
#include<stdio.h>
#include <stdbool.h>

/*!
 * enum type to identify result of parsing ini file
 */
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
    bool b;
    //float f;
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

/*!
 * Read ini configuration from given file pointer. Result is in res argument.
 * On error different than MEMORY_ERROR put line number in  parsed_config::number_of_sections for res variable.
 * \return state code from "enum parse_error"
 */
int readFromFile(FILE * f, struct parsed_config * res);

/*!
 * Open file and call readFromFile().
 */
int readFromPath(char * path, struct parsed_config * res);

/*!
 * Clean inside structures of parsed_config
 */
void deleteParsedConfig(struct parsed_config * c);

/*!
 * Return list of sections names which begins which given prefix;
 */
char ** getSectionsNamesByPrefix(char * prefix, struct parsed_config *c);

/*!
 * Check that section with given name exists
 */
int sectionExist(char * name, struct parsed_config * c);

/*!
 * Check that field in given section exists
 */
int fieldExist(char * section_name, char * field_name,  struct parsed_config * c);

/*!
 * Check that given field contains integral number
 * \return true if field exist and have correct type
 */
int isNumberField(char * section_name, char * field_name,  struct parsed_config * c);
/*!
 * Check that given field contains float
 * \return true if field exist and have correct type
 */
int isFloatField(char * section_name, char * field_name,  struct parsed_config * c);
/*!
 * Check that given field contains bool
 * \return true if field exist and have correct type
 */
int isBooleanField(char * section_name, char * field_name,  struct parsed_config * c);
/*!
 * Check that given field contains string
 * \return true if field exist and have correct type
 */
int isStringField(char * section_name, char * field_name,  struct parsed_config * c);
/*!
 * Check that given field is numeric
 * \return Field type if field exist, -1 otherwise
 */
int getFieldType(char * section_name, char * field_name,  struct parsed_config * c);

/*!
 * Get boolean value from given field. Behaviour undefined if fields is other type or do not exist.
 * In second case can cause Segmentation Fault
 */
bool getBoolValue(char * section_name, char * field_name,  struct parsed_config * c);
/*!
 * Get double value from given field. Behaviour undefined if fields is other type or do not exist.
 * In second case can cause Segmentation Fault
 */
double getFloatValue(char * section_name, char * field_name,  struct parsed_config * c);
/*!
 * Get integer value from given field. Behaviour undefined if fields is other type or do not exist.
 * In second case can cause Segmentation Fault
 */
int getNumericValue(char * section_name, char * field_name,  struct parsed_config * c);
/*!
 * Get string value from given field. Behaviour undefined if fields is other type or do not exist.
 * In second case can cause Segmentation Fault
 *
 * Return reference to string inside structure. Pleas do not modify.
 */
char const * getStringValue(char * section_name, char * field_name,  struct parsed_config * c);

//TODO think about interface to single section

/*!
 * Give access to RAW  section structure. Be careful.
 */
struct section_of_ini * getSection(char * section_name, struct parsed_config * c);

/*!
 * Give access to RAW  field structure. Be careful.
 */
struct fields_of_ini * getFieldFromSection(char * field_name, struct section_of_ini* section);

/*!
 * Print ini config in lexicographical order.
 */
void prettyPrint(FILE *f, struct parsed_config *c);

#endif //TIMOTHY_INI_PARSER_H

