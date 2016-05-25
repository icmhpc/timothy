//
// Created by Grzegorz Bokota on 18.02.16.
//

#include "ini_manipulator.h"
#include <algorithm>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string.h>
#include "../endian.h"



namespace timothy {
    namespace ini_manipulator {

        namespace {
            const std::map<std::string, size_t> type_size = {
                    {"int", sizeof(int)}, {"uint", sizeof(unsigned)}, {"float", sizeof(float)},
                    {"double", sizeof(double)}, {"char", sizeof(char)},
                    {"uint8", sizeof(uint8_t)}, {"int8", sizeof(int8_t)},
                    {"uint16", sizeof(uint16_t)}, {"int16", sizeof(int16_t)},
                    {"uint32", sizeof(uint32_t)}, {"int32", sizeof(int32_t)}
            };
            const std::string indent  = "    ";
        }

        std::string to_string(ini_fields_type type){
          switch (type){
            case NUMBER_FIELD: return "Integer";
            case FLOAT_FIELD: return "Double";
            case BOOLEAN_FIELD: return "Boolean";
            case STRING_FIELD: return "String";
            case UNSPECIFIED_FIELD: return "Unspecified";
            case STRUCT_FIELD: return "Struct";
          }
          return "Error [timothy::ini_manipulator::to_string]";
        }

        iniField::iniField(const std::string name, const int val) {
          this->name = name;
          this->value.i = val;
          this->type = NUMBER_FIELD;
        }

        iniField::iniField(const std::string name, const bool val) {
          this->name = name;
          this->value.b = val;
          this->type = BOOLEAN_FIELD;
        }

        iniField::iniField(const std::string name, const double val) {
          this->name = name;
          this->value.d = val;
          this->type = FLOAT_FIELD;
        }

        iniField::iniField(const std::string name, const char *val) {
          this->name = name;
          this->value.str = new std::string(val);
          this->type = STRING_FIELD;
        }

        iniField::iniField(const std::string name, const std::string val) {
          this->name = name;
          this->value.str = new std::string(val);
          this->type = STRING_FIELD;
        }

        iniField::~iniField() {
          if (this->type == STRING_FIELD){
            delete this->value.str;
          }
          if (this->type == STRUCT_FIELD){
            delete this->value.str;
            free(this->ptr);
          }
        }

        iniField::iniField(const iniField & f) : name(f.name), type(f.type) {
          switch (f.type){
            case NUMBER_FIELD: this->value.i = f.value.i; break;
            case FLOAT_FIELD: this->value.d = f.value.d; break;
            case BOOLEAN_FIELD: this->value.b = f.value.b; break;
            case STRING_FIELD: this->value.str = new std::string(*f.value.str); break;
            case UNSPECIFIED_FIELD:break;
            case STRUCT_FIELD:
              this->value.str = new std::string(*f.value.str);
              this->ptr = malloc(f.size);
              this->size = f.size;
              break;
          }
        }

        void flushComment(std::istream &s){
          if (s.peek() == ';' || s.peek() == '#')
            s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        bool isComment(int c){
          return c == ';' || c == '#';
        }

        bool isComment(std::istream &s){
          return isComment(s.peek());
        }

        iniField::iniField(std::istream &s) {
          std::string tmp_name;
          std::string tmp_value;
          ini_fields_type tmp_type = UNSPECIFIED_FIELD;
          bool escaped = false;
          while (isspace(s.peek())) s.get();
          while (isalnum(s.peek()) || s.peek() == '_') tmp_name += s.get();
          while (isspace(s.peek())) s.get();
          if (s.peek() != '=')
            throw fileFormatError("Field name have to contain only alphanumeric characters: " + tmp_name);
          while (isspace(tmp_name.back())) tmp_name.pop_back();
          if (s.eof() || tmp_name.length() == 0 )
            throw fileFormatError("invalid field format ");
          s.get(); // read '='
          while (isspace(s.peek())) s.get();

          if (s.peek() == '"' ){
            s.get();
            tmp_type = STRING_FIELD;
            while (true){
              if (escaped){
                std::string error("Unsupported escape characters ");
                switch (s.peek()){
                  case 'n' : tmp_value += "\n"; break;
                  case 't' : tmp_value += "\t"; break;
                  case 'r' : tmp_value += "\r"; break;
                  case '\"' : tmp_value += "\""; break;
                  case '\'' : tmp_value += "\'"; break;
                  case '\\' : tmp_value += "\\"; break;
                  case ' ' :
                    while (isblank(s.peek())) s.get();
                    if (s.peek() == '\r' || s.peek() == '\n')
                      s.get();
                    else
                      throw fileFormatError("After break line sequence \"\\ \" can bee only blank characters");
                    if (s.peek() == '\n') s.get();
                    while (isblank(s.peek())) s.get();
                    break;
                  default: error.push_back((char)s.peek()); throw fileFormatError(error);
                }
                s.get();
                escaped = false;
              } else {
                if (s.peek() == '\"'){
                  s.get();
                  break;
                }
                if (s.peek() ==  '\\'){
                  escaped = true;
                  s.get();
                  continue;
                }
                tmp_value.push_back((char) s.get());
              }
            }
          } else {
            if (s.peek() == '{'){ //TODO check pragma pack.
              s.get();
              tmp_type = STRUCT_FIELD;
              std::vector<std::pair<std::string, std::string>> struct_fields;
              std::string field_type;
              std::string field_value;
              while (s.peek() != '}' && s.peek() != EOF) {
                while (isspace(s.peek())) s.get();
                while (isalnum(s.peek()) || s.peek() == '_') { field_type += s.get(); }
                while (isblank(s.peek())) s.get();
                if (s.peek() != ':') {
                  std::cerr << field_type << " " << s.peek() << std::endl;
                  throw fileFormatError("struct field must have format \"type : value\"");
                }
                s.get();
                while (isblank(s.peek())) s.get();
                while (isprint(s.peek())) field_value += s.get();
                while (isspace(field_value.back())) field_value.pop_back();
                struct_fields.push_back(std::make_pair(field_type, field_value));
                field_type = "";
                field_value = "";
                while (isspace(s.peek())) s.get();
              }
              if (s.peek() != '}')
                throw  fileFormatError("Bad struct format");
              s.get();
              size_t struct_size = 0;
              for (auto it : struct_fields){
                try {
                  struct_size += type_size.at(it.first);
                } catch (std::out_of_range e){
                  throw fileFormatError("Unknown type " + it.first);
                }
              }
              this->size = struct_size;
              this->ptr = malloc(struct_size);
              tmp_value = "{\n";
              char * tmp_pointer = (char *) this->ptr;
              std::string parsed;
              for (auto it : struct_fields) {
                size_t pos = 0;
                if (it.first == "float"){
                  *(float *) tmp_pointer = (float) std::stod(it.second, &pos);
                  if (pos != it.second.length())
                    throw fileFormatError("Can not convert \"" + it.second + "\" to float");
                  parsed = std::to_string(*(float *) tmp_pointer);
                  tmp_pointer += sizeof(float)/sizeof(char);
                  tmp_value += indent + indent + it.first + " : " + parsed + "\n";
                  continue;
                }
                if (it.first == "double"){
                  *(double *) tmp_pointer = std::stod(it.second, &pos);
                  if (pos != it.second.length())
                    throw fileFormatError("Can not convert \"" + it.second + "\" to double");
                  parsed = std::to_string(*(double *) tmp_pointer);
                  tmp_pointer += sizeof(double)/sizeof(char);
                  tmp_value += indent + indent + it.first + " : " + parsed + "\n";
                  continue;
                }
                uint64_t swap;
                if (it.first[0] == 'u'){
                  unsigned long long val = std::stoull(it.second, &pos);
                  if (pos != it.second.length())
                    throw fileFormatError("Can not convert \"" + it.second + "\" to unsigned integer type");
                  parsed = std::to_string(val); //TODO casting on smaller numbers
                  swap = htobe64(val);
                } else {
                  long long val = std::stoll(it.second, &pos);
                  if (pos != it.second.length())
                    throw fileFormatError("Can not convert \"" + it.second + "\" to signed integer type");
                  parsed = std::to_string(val);
                  swap = htobe64(val);
                }
                char * array_swap = (char *) &swap;
                pos = 8 - type_size.at(it.first);
                switch (type_size.at(it.first)){
                  case 2 : be16toh( *(uint16_t *) &array_swap[6]); break;
                  case 4 : be32toh( *(uint32_t *) &array_swap[4]); break;
                  case 8 : be64toh( *(uint64_t *) &array_swap[0]); break;
                  default: break;
                }
                memcpy(tmp_pointer, array_swap+pos, type_size.at(it.first));
                tmp_pointer += type_size.at(it.first);
                tmp_value += indent + indent + it.first + " : " + parsed + "\n";

              }
              tmp_value += indent + "}";
              this->value.str = new std::string(tmp_value);
            } else {
              while (true) {
                if (s.peek() == '\n') {
                  s.get();
                  break;
                }
                tmp_value.push_back((char) s.get());
              }
            }
          }
          while (isspace(tmp_value.back())) tmp_value.pop_back();
          this->name = tmp_name;
          if (tmp_type != UNSPECIFIED_FIELD){
            this->type = tmp_type;
            switch (tmp_type){
              case NUMBER_FIELD:break;
              case FLOAT_FIELD:break;
              case BOOLEAN_FIELD:break;
              case STRING_FIELD: this->value.str = new std::string(tmp_value); break;
              case STRUCT_FIELD:break;
              case UNSPECIFIED_FIELD:break;
            }
            return;
          }
          if (tmp_value == "true" || tmp_value == "True" || tmp_value == "TRUE" ){
            this->type = BOOLEAN_FIELD;
            this->value.b = true;
            return;
          }
          if (tmp_value == "false" || tmp_value == "False" || tmp_value == "FALSE" ){
            this->type = BOOLEAN_FIELD;
            this->value.b = false;
            return;
          }
          if (isdigit(tmp_value[0]) || tmp_value[0] == '.' || tmp_value[0] == '-' || tmp_value[0] == '+'){
            bool has_dot = false;
            bool has_e = false;
            tmp_type = NUMBER_FIELD;
            auto ch = tmp_value.begin();
            if (*ch == '+' || *ch == '-')
              ch++;
            for (; ch != tmp_value.end(); ch++){
              if (isdigit(*ch))
                continue;
              if (*ch == '.'){
                if (!has_dot) {
                  has_dot = true;
                  tmp_type = FLOAT_FIELD;
                } else {
                  tmp_type = UNSPECIFIED_FIELD;
                  break;
                }
                continue;
              }
              if (*ch == 'e'){
                if (has_e){
                  tmp_type = UNSPECIFIED_FIELD;
                  break;
                } else {
                  has_e = true;
                  has_dot = true;
                  tmp_type = FLOAT_FIELD;
                  if (*(ch+1) == '+' || (*(ch+1) == '-')){
                    ch++;
                  }
                }
                continue;
              }
              tmp_type = UNSPECIFIED_FIELD;
              break;
            }
            if (tmp_type == NUMBER_FIELD){
              std::istringstream ss(tmp_value);
              ss >> this->value.i;
              this->type = NUMBER_FIELD;
              return;
            }
            if (tmp_type == FLOAT_FIELD){
              std::istringstream ss(tmp_value);
              ss >> this->value.d;
              this->type = FLOAT_FIELD;
              return;
            }
          }

          this->value.str = new std::string(tmp_value);
          this->type = STRING_FIELD;
        }


        bool iniField::isRestrictBoolean() const { return this->type == BOOLEAN_FIELD; }

        bool iniField::isRestrictFloat() const { return this->type == FLOAT_FIELD; }

        bool iniField::isRestrictString() const { return this->type == STRING_FIELD; }

        bool iniField::isRestrictNumeric() const { return this->type == NUMBER_FIELD; }

        bool iniField::isBoolean() const { return this->isRestrictBoolean() || isRestrictNumeric(); }

        bool iniField::isNumeric() const { return this->isRestrictNumeric(); }

        bool iniField::isString() const { return this->isRestrictString(); }

        bool iniField::isFloat() const { return this->isRestrictFloat() || this->isRestrictNumeric(); }

        bool iniField::isStruct() const {
          return this->type == STRUCT_FIELD;
        }


        bool iniField::getBoolValue() const {
          switch (this->type) {
            case NUMBER_FIELD :
              return (bool) this->value.i;
            case BOOLEAN_FIELD:
              return this->value.b;
            default:
              throw type_error("getBoolValue on non boolean compatybile field");
          }
        }

        long iniField::getIntValue() const {
          switch (this->type) {
            case NUMBER_FIELD :
              return this->value.i;
            default:
              throw type_error("getIntValue on non numeric field");
          }
        }

        double iniField::getFloatValue() const {
          switch (this->type) {
            case NUMBER_FIELD :
              return (double) this->value.i;
            case FLOAT_FIELD :
              return this->value.d;
            default:
              throw type_error("getFloatValue on non numeric field");
          }
        }

        std::string iniField::getStringValue() const {
          switch (this->type) {
            case STRING_FIELD :
              return *this->value.str;
            default:
              throw type_error("getStringValue on non string field");
          }
        }

        void iniField::setValue(const int i) {
          this->value.i = i;
          this->type = NUMBER_FIELD;
        }

        void iniField::setValue(const bool b) {
          this->value.b = b;
          this->type = BOOLEAN_FIELD;
        }

        void iniField::setValue(const double d) {
          this->value.d = d;
          this->type = FLOAT_FIELD;
        }

        void iniField::setValue(const std::string &s) {
          this->value.str = new std::string(s);
          this->type = STRING_FIELD;
        }

        /*bool iniField::operator<(const std::string &s) const {
          return  s > this->name;
        }*/

        /*bool iniField::operator<(const iniField &field) {
          return false;
        }*/

        type_error::type_error(std::string &what_arg) : logic_error(what_arg) { ; }

        type_error::type_error(const char *what_arg) : logic_error(what_arg) { ; }

        fileFormatError::fileFormatError(std::string s) : msg(s) { ; }

        const char *fileFormatError::what() const noexcept {
          return msg.c_str();
        }


        iniSection::iniSection(const std::string &name) : name(name) { ; }

        iniSection::iniSection(const std::string &name, const std::vector<iniField> fields) : name(name) {
          for (auto it : fields){
            this->fields[it.getName()] = it;
            this->fields_order.push_back(it.getName());
          }
          this->fields_order.shrink_to_fit();
        }

        void removeSpaces(std::istream &s){
          while (isspace(s.peek())) {s.get();};
        }

        iniSection::iniSection(std::istream &s) {
          if (s.get() != '['){
            throw fileFormatError("Section name should start from \'[\'");
          }
          int c;
          while ((c = s.get()) != ']' && (isalnum(c) || '_' == c)){
            this->name += c;
          }
          if (c != ']' && c != '_'){
            throw fileFormatError("Section name should contain only alphanumeric characters or \'_\'");
          }

          while (s.peek() != '[' && s.peek() != EOF){
            do {
              flushComment(s);
              removeSpaces(s);
            } while (isComment(s));
            if (s.peek() == '[' || s.peek() == EOF)
              break;
            auto field = iniField(s);
            this->fields_order.push_back(field.name);
            this->fields.insert(std::make_pair(field.name, field));
            removeSpaces(s);
          }
        }

        bool iniSection::hasField(const std::string name) const{
          return fields.count(name) > 0;
        }


        std::ostream & operator<<(std::ostream &o, const iniField &f) {
          o << f.name << " = "; //{" << to_string(f.type) << "} ";
          switch (f.type){
            case BOOLEAN_FIELD : o << f.value.b; break;
            case FLOAT_FIELD : o << std::showpoint << f.value.d; break;
            case NUMBER_FIELD : o << f.value.i; break;
            case STRING_FIELD : o << *f.value.str; break;
            case UNSPECIFIED_FIELD : o << "<UNSPECIFIED FIELD>"; break;
            case STRUCT_FIELD: o << *f.value.str; break;
          }
          return o;
        }



        bool iniSection::isBoolean(const std::string name) const {
          return fields.at(name).isBoolean();
        }

        bool iniSection::isFloat(const std::string name) const {
          return fields.at(name).isFloat();
        }

        bool iniSection::isString(const std::string name) const {
          return fields.at(name).isString();
        }

        bool iniSection::isNumeric(const std::string name) const {
          return fields.at(name).isNumeric();
        }
        bool iniSection::isStruct(const std::string name) const {
          return fields.at(name).isStruct();
        }

        bool iniSection::isRestrictBoolean(const std::string name) const {
          return fields.at(name).isRestrictBoolean();
        }

        bool iniSection::isRestrictFloat(const std::string name) const {
          return fields.at(name).isRestrictFloat();
        }

        bool iniSection::isRestrictString(const std::string name) const {
          return fields.at(name).isRestrictString();
        }

        bool iniSection::isRestrictNumeric(const std::string name) const {
          return fields.at(name).isRestrictNumeric();
        }

        long iniSection::getIntValue(const std::string name) const {
          return fields.at(name).getIntValue();
        }

        double iniSection::getFloatValue(const std::string name) const {
          return fields.at(name).getFloatValue();
        }

        std::string iniSection::getStringValue(const std::string name) const {
          return fields.at(name).getStringValue();
        }

        bool iniSection::getBoolValue(const std::string name) const {
          return fields.at(name).getBoolValue();
        }

        std::ostream & operator<<(std::ostream &o, const iniSection &f) {
          o << "[" << f.name << "]" << std::endl;
          for( auto & it : f.fields_order ){
            o << indent << f.fields.at(it) << std::endl;
          }
          return o;
        }

        iniConfiguration::iniConfiguration(std::istream &s) {
          while (s.peek() != EOF){
            if (isComment(s.peek())){
              flushComment(s);
              continue;
            }
            if (isspace(s.peek())){
              removeSpaces(s);
              continue;
            }
            if (s.peek() == '['){
              //std::cerr << "BUKA \"[\"" << std::endl;
              auto section = iniSection(s);
              this->sections_order.push_back(section.name);
              this->sections.insert(std::make_pair(section.name, section));
              continue;
            }
            if (s.peek() == EOF)
              return;
            std::string message("wrong value on section level ");
            message.push_back((char)s.get());
            throw fileFormatError(message);
          }
        }


        bool iniConfiguration::hasField(const std::string section_name, const std::string field_name) const{
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).hasField(field_name);
          return false;
        }

        bool iniConfiguration::isBoolean(const std::string section_name, const std::string field_name) const {
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).isBoolean(field_name);
          return false;
        }

        bool iniConfiguration::isFloat(const std::string section_name, const std::string field_name) const {
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).isFloat(field_name);
          return false;
        }

        bool iniConfiguration::isString(const std::string section_name, const std::string field_name) const {
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).isString(field_name);
          return false;
        }

        bool iniConfiguration::isNumeric(const std::string section_name, const std::string field_name) const {
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).isNumeric(field_name);
          return false;
        }

        bool iniConfiguration::isRestrictBoolean(const std::string section_name, const std::string field_name) const {
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).isRestrictBoolean(field_name);
          return false;
        }

        bool iniConfiguration::isRestrictFloat(const std::string section_name, const std::string field_name) const {
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).isRestrictFloat(field_name);
          return false;
        }

        bool iniConfiguration::isRestrictString(const std::string section_name, const std::string field_name) const {
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).isRestrictString(field_name);
          return false;
        }

        bool iniConfiguration::isRestrictNumeric(const std::string section_name, const std::string field_name) const {
          if (sections.count(section_name) > 0 )
            return sections.at(section_name).isRestrictNumeric(field_name);
          return false;
        }

        long iniConfiguration::getIntValue(const std::string section_name, const std::string field_name) const {
          return sections.at(section_name).getIntValue(field_name);
        }

        double iniConfiguration::getFloatValue(const std::string section_name, const std::string field_name) const {
          return sections.at(section_name).getFloatValue(field_name);
        }

        std::string iniConfiguration::getStringValue(const std::string section_name,
                                                     const std::string field_name) const {
          return sections.at(section_name).getStringValue(field_name);
        }

        bool iniConfiguration::getBoolValue(const std::string section_name, const std::string field_name) const {
          return sections.at(section_name).getBoolValue(field_name);
        }

        bool iniConfiguration::hasFieldsWithType(const std::vector<section_info> &data, std::stringstream &errors) const {
          bool any_error = false;
          for (auto & section_requirement : data){
            if (this->sections.count(section_requirement.first) == 0){
              any_error = true;
              errors << "Section "  << section_requirement.first << " do not exists" << std::endl;
              continue;
            }
            auto & section = this->sections.at(section_requirement.first);
            for (auto field_requirement = section_requirement.second.begin();
                    field_requirement != section_requirement.second.end(); field_requirement++){
              if (section.hasField(field_requirement->name)){
                bool local_any_error = false;
                switch (field_requirement->type){
                  case NUMBER_FIELD: local_any_error =  section.isNumeric(field_requirement->name); break;
                  case FLOAT_FIELD: local_any_error =  section.isFloat(field_requirement->name); break;
                  case BOOLEAN_FIELD: local_any_error = section.isBoolean(field_requirement->name); break;
                  case STRING_FIELD: local_any_error =  section.isString(field_requirement->name); break;
                  case STRUCT_FIELD: local_any_error =  section.isStruct(field_requirement->name); break;
                  case UNSPECIFIED_FIELD: local_any_error = true; break;
                }
                if (!local_any_error){
                  any_error = true;
                  errors << "Field " << section_requirement.first << "::" << field_requirement->name;
                  errors << " has wrong type " << to_string(section.fields.at(field_requirement->name).type);
                  errors << " instead of " << to_string(field_requirement->type) << std::endl;
                } else {
                  if (field_requirement->test != nullptr){
                    if (!field_requirement->test->check(section.fields.at(field_requirement->name))){
                      errors << field_requirement->test->error_message(
                              section_requirement.first + "::" + field_requirement->name) << std::endl;
                    }
                  }
                }
              } else if(field_requirement->required){
                any_error = true;
                errors << "Field " << section_requirement.first << "::" << field_requirement->name;
                errors << " does not exists" << std::endl;
              }
            }
          }
          return !any_error;
        }

        std::vector<std::string> iniConfiguration::getListOfSections() {
          return this->sections_order;
        }

        std::ostream & operator<<(std::ostream &o, const iniConfiguration &f) {
          for (auto & it : f.sections_order){
            o << f.sections.at(it);
            o << std::endl;
          }
          return o;
        }


    }
}



