//
// Created by czaki on 18.02.16.
//

#ifndef TIMOTHY_INI_MANIPULATOR_H
#define TIMOTHY_INI_MANIPULATOR_H

#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <istream>

namespace timothy {
    namespace ini_manipulator {
        class type_error : public std::logic_error{
        public:
            explicit type_error (std::string & what_arg);
            explicit type_error (const char * what_arg);
        };
        class fileFormatError : public std::exception {
            std::string msg;
        public:
            fileFormatError(std::string);
            virtual const char * what() const noexcept;

        };
        enum ini_fields_type {
            NUMBER_FIELD,
            FLOAT_FIELD,
            BOOLEAN_FIELD,
            STRING_FIELD,
            UNSPECIFIED_FIELD
        };

        union ini_data {
            int i;
            bool b;
            double d;
            std::string * str;
        };

        class iniField {
        private:
            std::string name;
            ini_data value;
            ini_fields_type type;
        public:
            iniField(const std::string name) : name(name),  type(UNSPECIFIED_FIELD) {};
            iniField(const std::string name, const int val);
            iniField(const std::string name, const bool val);
            iniField(const std::string name, const double val);
            iniField(const std::string name, const char* val);
            iniField(const std::string name, const std::string val);
            iniField(std::istream & s);
            virtual ~iniField();
            std::string getName() const { return name ;}
            bool isBoolean() const;
            bool isFloat()const;
            bool isString() const;
            bool isNumeric() const;
            bool isRestrictBoolean() const;
            bool isRestrictFloat() const;
            bool isRestrictString() const;
            bool isRestrictNumeric() const ;
            int getIntValue() const;
            double getFloatValue() const;
            std::string getStringValue() const;
            bool getBoolValue() const;
            void setValue(const int);
            void setValue(const bool);
            void setValue(const double);
            void setValue(const std::string &);
            friend std::ostream & operator<< (std::ostream & o, const iniField &f);
            //bool operator< (const std::string & ) const;
            //bool operator< (const iniField & ) const;
        };
        class iniSection {
            std::string name;
            std::map<std::string, iniField> fields;
            std::vector<std::string> fields_order;
        public:
            iniSection(const std::string & name);
            iniSection(const std::string & name, std::vector<iniField> fields);
            bool hasField(const std::string name) const;
            bool isBoolean(const std::string name) const;
            bool isFloat(const std::string name)const;
            bool isString(const std::string name) const;
            bool isNumeric(const std::string name) const;
            bool isRestrictBoolean(const std::string name) const;
            bool isRestrictFloat(const std::string name) const;
            bool isRestrictString(const std::string name) const;
            bool isRestrictNumeric(const std::string name) const ;
            int getIntValue(const std::string name) const;
            double getFloatValue(const std::string name) const;
            std::string getStringValue(const std::string name) const;
            bool getBoolValue(const std::string name) const;
            template<typename T>
            void setValue(const std::string name, const T & t ) {
              return fields.at(name).setValue(t);
            }


        };

        class iniConfiguration {
            std::map<std::string, iniSection> sections;
            std::vector<std::string> sections_order;
        public:
            bool hasField(const std::string section_name, const std::string field_name) const;
            bool isBoolean(const std::string section_name, const std::string field_name) const;
            bool isFloat(const std::string section_name, const std::string field_name) const;
            bool isString(const std::string section_name, const std::string field_name) const;
            bool isNumeric(const std::string section_name, const std::string field_name) const;
            bool isRestrictBoolean(const std::string section_name, const std::string field_name) const;
            bool isRestrictFloat(const std::string section_name, const std::string field_name) const;
            bool isRestrictString(const std::string section_name, const std::string field_name) const;
            bool isRestrictNumeric(const std::string section_name, const std::string field_name) const;
            int getIntValue(const std::string section_name, const std::string field_name) const;
            double getFloatValue(const std::string section_name, const std::string field_name) const;
            std::string getStringValue(const std::string section_name, const std::string field_name) const;
            bool getBoolValue(const std::string section_name, const std::string field_name) const;
            template<typename T>
            void setValue(const std::string section_name, const std::string field_name, const T & t ) {
              return sections.at(section_name).setValue(field_name,t);
            }
        };
    }
}
#endif //TIMOTHY_INI_MANIPULATOR_H
