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
            STRUCT_FIELD,
            UNSPECIFIED_FIELD
        };

        union ini_data {
            long i;
            bool b;
            double d;
            std::string * str;
        };

        class iniField {
        private:
            std::string name;
            ini_data value;
            ini_fields_type type;
            void  * ptr;
            size_t size;
        public:
            iniField(const std::string name) : name(name),  type(UNSPECIFIED_FIELD) {};
            iniField() : name("empty field"),  type(UNSPECIFIED_FIELD) {};
            iniField(const std::string name, const int val);
            iniField(const std::string name, const bool val);
            iniField(const std::string name, const double val);
            iniField(const std::string name, const char* val);
            iniField(const std::string name, const std::string val);
            iniField(const iniField & f);
            iniField(std::istream & s);
            virtual ~iniField();
            std::string getName() const { return name ;}
            bool isBoolean() const;
            bool isFloat()const;
            bool isString() const;
            bool isNumeric() const;
            bool isStruct() const;
            bool isRestrictBoolean() const;
            bool isRestrictFloat() const;
            bool isRestrictString() const;
            bool isRestrictNumeric() const ;
            long getIntValue() const;
            double getFloatValue() const;
            std::string getStringValue() const;
            bool getBoolValue() const;
            void setValue(const int);
            void setValue(const bool);
            void setValue(const double);
            void setValue(const std::string &);
            ini_fields_type getType() { return type;};
            friend std::ostream & operator<< (std::ostream & o, const iniField &f);
            //bool operator< (const std::string & ) const;
            //bool operator< (const iniField & ) const;
            friend class iniSection;
            friend class iniConfiguration;
        };
        class iniSection {
            std::string name;
            std::map<std::string, iniField> fields;
            std::vector<std::string> fields_order;
        public:
            iniSection(const std::string & name);
            iniSection(const std::string & name, std::vector<iniField> fields);
            iniSection(std::istream & s);
            bool hasField(const std::string name) const;
            bool isBoolean(const std::string name) const;
            bool isFloat(const std::string name)const;
            bool isString(const std::string name) const;
            bool isNumeric(const std::string name) const;
            bool isStruct(const std::string name) const;
            bool isRestrictBoolean(const std::string name) const;
            bool isRestrictFloat(const std::string name) const;
            bool isRestrictString(const std::string name) const;
            bool isRestrictNumeric(const std::string name) const;
            long getIntValue(const std::string name) const;
            double getFloatValue(const std::string name) const;
            std::string getStringValue(const std::string name) const;
            bool getBoolValue(const std::string name) const;
            template<typename T>
            void setValue(const std::string name, const T & t ) {
              return fields.at(name).setValue(t);
            }
            friend class iniConfiguration;
            friend std::ostream & operator<< (std::ostream & o, const iniSection &f);
        };

        class iniConfiguration {
            std::map<std::string, iniSection> sections;
            std::vector<std::string> sections_order;
        public:
            class tester {
            public:
                virtual bool check(const iniField &) = 0;
                virtual std::string error_message(std::string name) = 0;
            };
            class field_info{
            public:
                std::string name;
                ini_fields_type type;
                bool required;
                tester * test;
                field_info(std::string _n, ini_fields_type _t, bool _r, tester * t)
                        : name(_n), type(_t), required(_r), test(t) {};
            };

            typedef std::pair<const std::string, const std::vector<field_info > > section_info;
            iniConfiguration(std::istream & s);
            iniConfiguration() {};
            bool hasField(const std::string section_name, const std::string field_name) const;
            bool isBoolean(const std::string section_name, const std::string field_name) const;
            bool isFloat(const std::string section_name, const std::string field_name) const;
            bool isString(const std::string section_name, const std::string field_name) const;
            bool isNumeric(const std::string section_name, const std::string field_name) const;
            bool isRestrictBoolean(const std::string section_name, const std::string field_name) const;
            bool isRestrictFloat(const std::string section_name, const std::string field_name) const;
            bool isRestrictString(const std::string section_name, const std::string field_name) const;
            bool isRestrictNumeric(const std::string section_name, const std::string field_name) const;
            long getIntValue(const std::string section_name, const std::string field_name) const;
            double getFloatValue(const std::string section_name, const std::string field_name) const;
            std::string getStringValue(const std::string section_name, const std::string field_name) const;
            bool getBoolValue(const std::string section_name, const std::string field_name) const;
            std::vector<std::string> getListOfSections();
            template<typename T>
            void setValue(const std::string section_name, const std::string field_name, const T & t ) {
              return sections.at(section_name).setValue(field_name,t);
            }
            bool hasFieldsWithType(const std::vector<section_info> &data, std::stringstream &errors) const;
            friend std::ostream & operator<< (std::ostream & o, const iniConfiguration &f);

        };


    }
}
#endif //TIMOTHY_INI_MANIPULATOR_H
