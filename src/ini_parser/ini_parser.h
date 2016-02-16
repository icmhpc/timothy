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
int readFromPath(const char *path, struct parsed_config *res);

/*!
 * Function used to free memory after end of parser use.
 */
void deleteParsedConfig(struct parsed_config * c);

/*!
 * Return list of sections names which begins which given prefix;
 */
char ** getSectionsNamesByPrefix(const char *prefix, const struct parsed_config *c);

/*!
 * Check if section with given name exists
 */
int sectionExist(const char *name, const struct parsed_config *c);

/*!
 * Check if field in given section exists
 */
int fieldExist(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Check if given field contains integral number
 * \return true if field exist and have correct type
 */
int isNumberField(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Check if given field contains float
 * \return true if field exist and have correct type
 */
int isFloatField(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Check if given field contains bool
 * \return true if field exist and have correct type
 */
int isBooleanField(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Check if given field contains string
 * \return true if field exist and have correct type
 */
int isStringField(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Check if given field is numeric
 * \return Field type (enum ini_fields_type) if field exist, -1 otherwise
 */
int getFieldType(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Get boolean value from given field. Behaviour undefined if fields is other type or do not exist.
 * In second case can cause Segmentation Fault
 */
bool getBoolValue(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Get double value from given field. Behaviour undefined if fields is other type or do not exist.
 * In second case can cause Segmentation Fault
 */
double getFloatValue(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Get integer value from given field. Behaviour undefined if fields is other type or do not exist.
 * In second case can cause Segmentation Fault
 */
int getNumericValue(const char *section_name, const char *field_name, const struct parsed_config *c);

/*!
 * Get string value from given field. Behaviour undefined if fields is other type or do not exist.
 * In second case can cause Segmentation Fault
 *
 * Return reference to string inside structure. Please do not modify.
 */
char const * getStringValue(const char *section_name, const char *field_name, const struct parsed_config *c);

//TODO think about interface to single section

/*!
 * Give access to RAW section structure. Be careful.
 */
struct section_of_ini * getSection(const char *section_name, struct parsed_config *c);

/*!
 * Give access to RAW field structure. Be careful.
 */
struct fields_of_ini * getFieldFromSection(const char *field_name, struct section_of_ini *section);

/*!
 * Print ini config in lexicographical order.
 */
void prettyPrint(FILE *f, const struct parsed_config *c);

#endif //TIMOTHY_INI_PARSER_H

