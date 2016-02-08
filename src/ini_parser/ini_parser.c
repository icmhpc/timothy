//
// Created by czaki on 02.02.16.
//

#include "ini_parser.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

static const size_t max_buff_size = 400;
static const size_t max_buff_value_size = 4000;

void flushComment(FILE * f){
  int character;
  do {
    character = getc(f);
  } while (character != EOF && character != '\n');
  if (character == '\n')
    ungetc('\n', f);
}

void deleteFieldsList(struct fields_of_ini * fields, size_t size){
  for(size_t i =0; i < size; i++){
    free(fields[i].name);
    if (STRING_FIELD == fields[i].type){
      free(fields[i].data.str);
    }
  }
}

void deleteSectionList(struct section_of_ini * section_list, size_t size){
  for(size_t i = 0; i < size; i++){
    deleteFieldsList(section_list[i].fields, section_list[i].number_of_fields);
    free(section_list[i].fields);
    free((section_list[i].name));
  };
}

void deleteParsedConfig(struct parsed_config * c){
  if (NULL == c->sections_list)
    return;
  deleteSectionList(c->sections_list, c->number_of_sections);
  free(c->sections_list);
}


char * readSectionName(FILE *f, enum parse_error * err){
  char buff[max_buff_size];
  int c;
  size_t pos = 0;
  do {
    buff[pos] = (char) (c = getc(f));
    pos += 1;
  } while (c != ']' && c != '\n' && c != EOF && pos < max_buff_size);
  if (c == EOF || c == '\n'){
    *err = ERROR_IN_SECTION_NAME;
    return NULL;
  }
  if (pos > max_buff_size || (pos == max_buff_size && buff[max_buff_size - 1] != ']')) {
    *err = TO_LONG_SECTION_NAME;
    return NULL;
  }
  // remove
  buff[pos-1] = 0;
  char * res = calloc(pos, sizeof(char));
  strcpy(res, buff+1);
  return res;
};

enum parse_error readField(FILE * f, struct fields_of_ini * field, int *line_num){
  char name_buffer[max_buff_size];
  char value_buffer[max_buff_value_size];
  size_t pos = 0;
  while (isalnum(name_buffer[pos] = getc(f))){
    pos += 1;
    if (pos >= max_buff_size){
      return TO_LONG_FIELD_NAME;
    }
  }
  name_buffer[pos] = '\0';

  while (isblank(value_buffer[0] = getc(f))){};
  if (value_buffer[0] != '='){
    return BAD_FILE_FORMAT;
  }
  while (isblank(value_buffer[0] = getc(f))){};
  pos = 0;
  do {
    if (value_buffer[pos] == '\\'){
      // read every blank characters
      while (isblank(value_buffer[pos] = getc(f))) { };

      if ('\r' == value_buffer[pos])
        value_buffer[pos] = (char) getc(f);
      if ('\n' == value_buffer[pos]){
        *line_num += 1;
        // read every blank characters
        while (isblank(value_buffer[pos] = getc(f))) { };
      } else {
        return BAD_FILE_FORMAT;
      }
      continue;
    }
    if (isblank(value_buffer[pos]) && isblank(value_buffer[(pos > 0 ? pos - 1 : 0)])){
      continue;
    }
    if ((value_buffer[pos] == ';' || value_buffer[pos] == '#' ) && value_buffer[0] != '"'){
      ungetc(value_buffer[pos], f);
      flushComment(f);
      pos -=1;
    }
    pos++;
    if (pos >= max_buff_value_size){
      return TO_LONG_FIELD_VALUE;
    }
    value_buffer[pos] = (char) getc(f);
  } while (value_buffer[pos] != '\n' && value_buffer[pos] != EOF);
  *line_num += 1;
  value_buffer[pos] = '\0';
  pos -= 1;
  while (pos >0 && isblank(value_buffer[pos])) {value_buffer[pos] = '\0'; pos -= 1;};
  if (value_buffer[0] == '"'){
    field->type = STRING_FIELD;
    if (pos !=0 && value_buffer[pos] == '"'){
      value_buffer[pos] = '\0';
      field->data.str = malloc(strlen(value_buffer + 1) + 1);
      strcpy(field->data.str, value_buffer + 1);
    } else {
      //TODO Maybe error?
      field->data.str = malloc(strlen(value_buffer) + 1);
      strcpy(field->data.str, value_buffer);
    }
    goto name_copy;
  }


  if (isalpha(value_buffer[0])){
    if ('T' == value_buffer[0] && (strcmp(value_buffer, "TRUE") || strcmp(value_buffer, "True"))){
      field->type = BOOLEAN_FIELD;
      field->data.i = true;
      goto name_copy;
    }
    if ('F' == value_buffer[0] && (strcmp(value_buffer, "FALSE") || strcmp(value_buffer, "False"))){
      field->type = BOOLEAN_FIELD;
      field->data.i = false;
      goto name_copy;
    }
    if ('t' == value_buffer[0] && strcmp(value_buffer, "true")){
      field->type = BOOLEAN_FIELD;
      field->data.i = true;
      goto name_copy;
    }
    if ('f' == value_buffer[0] && strcmp(value_buffer, "false")){
      field->type = BOOLEAN_FIELD;
      field->data.i = false;
      goto name_copy;
    }
  }
  {
    int int_res;
    double float_res;
    int len, len2;
    sscanf(value_buffer, "%lf%n", &float_res, &len);
    sscanf(value_buffer, "%d%n", &int_res, &len2);
    if (strlen(value_buffer) == (size_t) len){
      if (strlen(value_buffer) == (size_t) len2){
        field->type = NUMBER_FIELD;
        field->data.i = int_res;
      } else {
        field->type = FLOAT_FIELD;
        field->data.d = float_res;
      }
      goto name_copy;
    }
  }

  field->type = STRING_FIELD;
  field->data.str = malloc(strlen(value_buffer) + 1);
  strcpy(field->data.str, value_buffer);
name_copy:
  field->name = malloc(strlen(name_buffer) + 1);
  strcpy(field->name, name_buffer);

  return OK;
};

int fieldCompare(const void * v1, const void * v2){
  struct fields_of_ini * f1 = (struct fields_of_ini *) v1;
  struct fields_of_ini * f2 = (struct fields_of_ini *) v2;
  return strcmp(f1->name, f2->name);
};

int sectionCompare(const void * v1, const void * v2){
  struct section_of_ini * s1 = (struct section_of_ini *) v1;
  struct section_of_ini * s2 = (struct section_of_ini *) v2;
  return strcmp(s1->name, s2->name);
}

int sectionNameCompare(const void * v1, const void * v2){
  char * name  = (char *) v1;
  struct section_of_ini * s = (struct section_of_ini *) v2;
  return strcmp(name, s->name);
}

int fieldNameCompare(const void * v1, const void * v2){
  char * name  = (char *) v1;
  struct fields_of_ini * s = (struct fields_of_ini *) v2;
  return strcmp(name, s->name);
}

enum parse_error readSection(FILE * f, struct section_of_ini * res, int *line_num){
  enum parse_error err = OK;
  res->number_of_fields = 0;
  res->fields = NULL;
  res->name = readSectionName(f, &err);
  if (res->name == NULL){
    free(res);
    return err;
  }
  size_t buffer_size = 10;
  //size_t pos = 0;
  struct fields_of_ini * fields;
  struct fields_of_ini * realloc_fields;
  fields = calloc(buffer_size, sizeof(struct fields_of_ini));
  size_t fields_num = 0;
  int c;
  do {
    switch (c = getc(f)) {
      case ';' : flushComment(f); break;
      case '#' : flushComment(f); break;
      case '\r': break;
      case '\n': line_num += 1; break;
      case '[' : ungetc(c, f); break;
      default:
        if (isalnum(c)) {
          ungetc(c,f);
          if (fields_num >= buffer_size){
#ifdef INI_LONG_SECTIONS
            //TODO Add to documentation
            buffer_size *= 2;
#else
            buffer_size += 10;
#endif /* INI_LONG_SECTIONS */
            realloc_fields = realloc(fields, buffer_size * sizeof(struct fields_of_ini));
            if (realloc_fields != NULL){
              fields = realloc_fields;
            } else {
              free(res->name);
              deleteFieldsList(fields, fields_num-1);
              free(fields);
              return MEMORY_ERROR;
            }
          }
          if ((err = readField(f, &fields[fields_num], line_num)) != OK){
            free(res->name);
            deleteFieldsList(fields, fields_num);
            return err;
          };
          fields_num += 1;

        }
    }
  } while (c != EOF && c != '[');
  qsort(fields, fields_num, sizeof(struct fields_of_ini), fieldCompare);
  res->fields = realloc(fields, fields_num * sizeof(struct fields_of_ini));
  res->number_of_fields = fields_num;
  return OK;
};


int readFromFile(FILE * f, struct parsed_config * res){
  res->sections_list = NULL;
  res->number_of_sections = 0;
  int character;
  size_t number_of_section = 0;
  size_t buffer_size = 10;
  struct section_of_ini * section_list;
  struct section_of_ini * realloc_section_list;
  int err;
  int line_number = 0;

  section_list = calloc(buffer_size, sizeof(struct section_of_ini));

  while ((character = getc(f)) != EOF){
    switch (character) {
      case ' ' : break;
      case ';' : flushComment(f); break;
      case '#' : flushComment(f); break;
      case '\n' : line_number += 1;
      case '[' :
        ungetc(character, f);
        if (number_of_section >= buffer_size) {
#ifdef INI_MANY_SECTIONS
          buffer_size *= 2
#else
          buffer_size += 10;
#endif
          realloc_section_list = realloc(section_list, buffer_size);
          if (realloc_section_list != NULL){
            section_list = realloc_section_list;
          } else {
            deleteSectionList(section_list, number_of_section);
            free(section_list);
            return MEMORY_ERROR;
          }
        }

        err = readSection(f, &section_list[number_of_section], &line_number);
        if (err != OK) {
          deleteSectionList(section_list, number_of_section);
          free(section_list);
          return err;
        }
        number_of_section += 1;

        break;
      default:break;
    }
  }
  qsort(section_list,number_of_section, sizeof(struct section_of_ini), sectionCompare);
  realloc_section_list = realloc(section_list, number_of_section * sizeof(struct section_of_ini));
  if (realloc_section_list == NULL) {
    deleteSectionList(section_list, number_of_section);
    return MEMORY_ERROR;
  }
  res->number_of_sections = number_of_section;
  res->sections_list = realloc_section_list;
  return OK;
}

int readFromPath(char * path, struct parsed_config * res)
{
  FILE * f = fopen(path, "r");
  if (f == NULL)
    return NO_FILE;
  int err_code = readFromFile(f, res);
  fclose(f);
  return err_code;
};

int sectionExist(char * name, struct parsed_config * c) {
  return NULL != bsearch(name, c->sections_list, c->number_of_sections, sizeof(struct section_of_ini), sectionNameCompare);
}

struct fields_of_ini * getField(char * section_name, char * field_name,  struct parsed_config * c){
  struct section_of_ini * sect = (struct section_of_ini *)
          bsearch(section_name, c->sections_list, c->number_of_sections, sizeof(struct section_of_ini), sectionNameCompare);
  if (sect == NULL){
    return NULL;
  }
  return bsearch(field_name, sect->fields, sect->number_of_fields, sizeof(struct fields_of_ini), fieldNameCompare);

}

int fieldExist(char * section_name, char * field_name,  struct parsed_config * c){
  return  NULL != getField(section_name, field_name, c);
}

void prettyPrint(FILE *f, struct parsed_config *c){
  for (size_t i = 0; i < c->number_of_sections; i++){
    fprintf(f,"[%s]\n", c->sections_list[i].name);
    for (size_t j = 0; j < c->sections_list[i].number_of_fields; j++){
      struct fields_of_ini * field = &c->sections_list[i].fields[j];
      fprintf(f, "    %s = ", field->name);
      switch (field->type){
        case NUMBER_FIELD : fprintf(f, "%d\n", field->data.i); break;
        case FLOAT_FIELD : fprintf(f, "%e\n", field->data.d); break;
        case BOOLEAN_FIELD : fprintf(f, "%s\n", (field->data.i ? "TRUE" : "FALSE")); break;
        case STRING_FIELD : fprintf(f, "\"%s\"\n", field->data.str); break; //TODO break long lines
        default:
          fprintf(f, "UNSUPPORTED FORMAT\n");
          break;
      }
    }
  }
};




