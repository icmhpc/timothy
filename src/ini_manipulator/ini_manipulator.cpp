//
// Created by Grzegorz Bokota on 18.02.16.
//

#include "ini_manipulator.h"
#include <algorithm>
#include <sstream>

const static std::string indent = "    ";


namespace timothy {
    namespace ini_manipulator {
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
        }

        iniField::iniField(std::istream &s) {
          std::stringbuf st;
          std::string tmp_name;
          std::string tmp_value;
          while (isspace(s.peek())) s.get();
          s.get(st, '=');
          tmp_name = st.str();
          while (isspace(tmp_name.back())) tmp_name.pop_back();
          if (s.eof() || tmp_name.length() == 0 )
            throw fileFormatError("invalid field format ")
          s.get(); // read '='
          while (isspace(s.peek())) s.get();

          if (s.peek() == '"' ){

          } else {
            s.get(st, '/');
            tmp_value = st.str();
            while (s.peek() == '/') {
              s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
              tmp_value += "\n";
              s.get(st, '/');
              tmp_value = st.str();
            }
          }

        }


        bool iniField::isRestrictBoolean() const { return this->type == BOOLEAN_FIELD; }

        bool iniField::isRestrictFloat() const { return this->type == FLOAT_FIELD; }

        bool iniField::isRestrictString() const { return this->type == STRING_FIELD; }

        bool iniField::isRestrictNumeric() const { return this->type == NUMBER_FIELD; }

        bool iniField::isBoolean() const { return this->isRestrictBoolean() || isRestrictNumeric(); }

        bool iniField::isNumeric() const { return this->isRestrictNumeric(); }

        bool iniField::isString() const { return this->isRestrictString(); }

        bool iniField::isFloat() const { return this->isRestrictFloat() || this->isRestrictNumeric(); }

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

        int iniField::getIntValue() const {
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

        const char *fileFormatError::what() const {
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

        bool iniSection::hasField(const std::string name) const{
          return fields.count(name) > 0;
        }


        std::ostream &timothy::ini_manipulator::operator<<(std::ostream &o, const iniField &f) {
          o << f.name << " = ";
          switch (f.type){
            case BOOLEAN_FIELD : o << f.value.b; break;
            case FLOAT_FIELD : o << f.value.d; break;
            case NUMBER_FIELD : o << f.value.i; break;
            case STRING_FIELD : o << *f.value.str; break;
            case UNSPECIFIED_FIELD : o << "<UNSPECIFIED FIELD>"; break;
          }
          o << std::endl;
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

        int iniSection::getIntValue(const std::string name) const {
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

        int iniConfiguration::getIntValue(const std::string section_name, const std::string field_name) const {
          return sections.at(section_name).getIntValue(field_name);
        }

        double iniConfiguration::getFloatValue(const std::string section_name, const std::string field_name) const {
          return sections.at(section_name).getFloatValue(field_name);
        }

        std::string iniConfiguration::getStringValue(const std::string section_name,
                                                     const std::string field_name) const {
          return return sections.at(section_name).getStringValue(field_name
        }

        bool iniConfiguration::getBoolValue(const std::string section_name, const std::string field_name) const {
          return return sections.at(section_name).getBoolValue(field_name
        }


    }
}



