#ifndef TIMOTHY_DICTS
#define TIMOTHY_DICTS

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

struct str_uint16_pair{
  char * first;
  uint16_t second;
};
typedef struct str_uint16_pair str_uint16_pair;

struct str_uint16_dict{
  str_uint16_pair * members;
  size_t size;
};

typedef struct str_uint16_dict str_uint16_dict;

int compare_str_uint16_pair(const void * value1, const void * value2){
  str_uint16_pair * s1 = (str_uint16_pair *) value1;
  str_uint16_pair * s2 = (str_uint16_pair *) value2;
  return strcmp(s1->first, s2->first);
}

int compare_str_and_str_uint16_pair(const void * key, const void * value){
  char * s1 = (char *) key;
  str_uint16_pair * s2 = (str_uint16_pair *) value;
  return strcmp(s1, s2->first);
}

/*!
 * Function used to create constant dictionary
 */
str_uint16_dict * create_dict(str_uint16_pair * in_data, size_t size){
  size_t strings_size =0 ;
  size_t i;
  str_uint16_dict * result;
  str_uint16_pair * data = malloc(sizeof(str_uint16_pair) * size);
  data = memcpy(data, in_data, sizeof(str_uint16_pair) * size);
  qsort(data, size, sizeof(str_uint16_pair), compare_str_uint16_pair);
  char * strings;
  for(i=0; i < size; i++){
    strings_size += 1 + strlen(data[i].first);
  }
  result = (str_uint16_dict *) malloc(sizeof(str_uint16_dict) + size * sizeof(str_uint16_pair) + strings_size * sizeof(char));
  result->size = size;
  result->members = (str_uint16_pair *) (result + 1);
  strings = (char *) (result->members + size); 
  for(i=0; i < size; i++){
    result->members[i].second = data[i].second;
    strcpy(strings, data[i].first);
    result->members[i].first = strings;
    strings += 1 + strlen(data[i].first);
  }
  free(data);
  return result;
};

/*!
 * function used to get value from str_uint_dict return -1 on fail
 */
int32_t get_value(str_uint16_dict * dict, char* key){
  str_uint16_pair * res = bsearch(key, dict->members, dict->size, sizeof(str_uint16_pair), compare_str_and_str_uint16_pair);
  if (res == NULL)
    return  (int32_t) -1;
  return (int32_t) res->second;
}

/*!
 * Function to calculate total size of dictionary. Useful to send wia MPI or network
 */
size_t sizeof_dict(str_uint16_dict * dict){
  size_t result = sizeof(str_uint16_dict);
  result += sizeof(str_uint16_pair) * sizeof(str_uint16_pair);
  for(size_t i = 0; i < dict->size; i++){
    result += strlen(dict->members[i].first);
  }
  return result;
}

/*!
 * Function to set properly pointers inside struct
 */
void regenerate_dict(str_uint16_dict * dict){
  dict->members = (str_uint16_pair *) (dict + 1);
  char * strings = (char*) (dict->members + dict->size);
  for (size_t i = 0; i < dict->size; i++){
    dict->members[i].first = strings;
    strings += strlen(strings) + 1;
  }
}

#endif /* ifndef TIMOTHY_DICTS */
