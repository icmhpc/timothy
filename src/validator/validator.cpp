//
// Created by Grzegorz Bokota on 11.02.16.
//


#include "validator.h"
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>

static const std::string error = "[Error] ";
static const std::string warning = "[Warning] ";

namespace timothy {
    namespace validator {
        namespace mp= timothy::ini_manipulator;
        namespace {
            typedef mp::iniConfiguration::tester tester;
            class isPowerOfTwo : public tester{
            public:
                virtual bool check (const mp::iniField & f){
                  auto num = f.getIntValue();
                  return num != 0 && !(num & (num - 1));
                }
                virtual std::string error_message(std::string name){
                  return  "Value of field " + name + " has to be power of two";
                }
            };
            class is2or3 : public tester{
            public:
                virtual bool check (const mp::iniField & f){
                  auto num = f.getIntValue();
                  return num == 2 || num == 3;
                }
                virtual std::string error_message(std::string name){
                  return  "Value of field " + name + " has to be 2 or 3";
                }
            };
            class activeOrPassive : public tester{
            public:
                virtual bool check (const mp::iniField & f){
                  auto str = f.getStringValue();
                  return str == "active" || str == "passive";
                }
                virtual std::string error_message(std::string name){
                  return  "Value of field " + name + " has to be \"active\" or \"passive\"";
                }
            };
            class nonNegativeFloat : public tester{
            public:
                virtual bool check (const mp::iniField & f){
                  auto fl = f.getFloatValue();
                  return fl >= 0;
                }
                virtual std::string error_message(std::string name){
                  return  "Value of field " + name + " has to be non negative number";
                }
            };
            class nonNegativeNumber : public tester{
            public:
                virtual bool check (const mp::iniField & f){
                  auto fl = f.getIntValue();
                  return fl >= 0;
                }
                virtual std::string error_message(std::string name){
                  return  "Value of field " + name + " has to be non negative number";
                }
            };
            isPowerOfTwo isPowerOfTwo;
            is2or3 is2or3;
            activeOrPassive activeOrPassive;
            nonNegativeFloat nonNegativeFloat;
            nonNegativeNumber nonNegativeNumber;

            std::vector<mp::iniConfiguration::field_info> global_check = {
                    {"size_x",              mp::NUMBER_FIELD, true, &isPowerOfTwo},
                    {"size_y",              mp::NUMBER_FIELD, true, &isPowerOfTwo},
                    {"size_z",              mp::NUMBER_FIELD, true, &isPowerOfTwo},
                    {"max_number_of_cells", mp::NUMBER_FIELD, true, &nonNegativeNumber},
                    {"dimension", mp::NUMBER_FIELD, true, &is2or3},
                    {"seconds_per_step", mp::NUMBER_FIELD, true, nullptr}

            };
            std::vector<mp::iniConfiguration::field_info> active_cell_check = {
                    {"type",     mp::STRING_FIELD, true, &activeOrPassive},
                    {"g1_phase", mp::FLOAT_FIELD,  true, &nonNegativeFloat},
                    {"s_phase",  mp::FLOAT_FIELD,  true, &nonNegativeFloat},
                    {"g2_phase", mp::FLOAT_FIELD,  true, &nonNegativeFloat},
                    {"m_phase",  mp::FLOAT_FIELD,  true, &nonNegativeFloat},
                    {"pre_function", mp::STRING_FIELD, false, nullptr},
                    {"pre_function_data", mp::STRUCT_FIELD, false, nullptr},
                    {"post_function", mp::STRING_FIELD, false, nullptr},
                    {"post_function_data", mp::STRUCT_FIELD, false, nullptr}
            };

            std::vector<mp::iniConfiguration::field_info> passive_cell_check;
            std::vector<mp::iniConfiguration::field_info> active_env_check = {
                    {"dc",               mp::FLOAT_FIELD, true, &nonNegativeFloat},
                    {"bc",               mp::FLOAT_FIELD, true, &nonNegativeFloat},
                    {"ic_mean",          mp::FLOAT_FIELD, true, &nonNegativeFloat},
                    {"ic_var",           mp::FLOAT_FIELD, true, &nonNegativeFloat},
                    {"lambda",           mp::FLOAT_FIELD, true, &nonNegativeFloat},
                    {"critical_level_1", mp::FLOAT_FIELD, true, &nonNegativeFloat},
                    {"critical_level_2", mp::FLOAT_FIELD, true, &nonNegativeFloat}
            };
            std::vector<mp::iniConfiguration::field_info> passive_env_check;
        }


        std::pair<bool, std::string> validate_config(mp::iniConfiguration & c) {
          std::stringstream errors;
          bool is_error = false;
          auto names = c.getListOfSections();
          decltype(names) cells_active;
          decltype(names) cells_passive;
          decltype(names) environments_active;
          decltype(names) environments_passive;
          for (auto & it : names){
            if ("CELL_" == it.substr(0,5)){
              if (c.hasField(it, "type") && c.isString(it, "type")){
                auto type_val = c.getStringValue(it, "type");
                if ("active" == type_val){
                  cells_active.push_back(it);
                  continue;
                }
                if ("passive" == type_val){
                  cells_passive.push_back(it);
                  continue;
                }
                errors << "Cell " << it << " has wrong cell type set. It should be \"passive\" or \"active\"" << std::endl;
                is_error = true;
              } else {
                errors << "Cell " << it << " has no cell type set. It should be \"type = passive\" or \"type = active\"" << std::endl;
                is_error = true;
              }
            }
            if ("ENV_" == it.substr(0,4)){
              if (c.hasField(it, "type") && c.isString(it, "type")){
                auto type_val = c.getStringValue(it, "type");
                if ("active" == type_val){
                  environments_active.push_back(it);
                  continue;
                }
                if ("passive" == type_val){
                  environments_passive.push_back(it);
                  continue;
                }
                errors << "Environment " << it << " has wrong cell type set. It should be \"passive\" or \"active\"" << std::endl;
                is_error = true;
              } else {
                errors << "Environment " << it << " has no cell type set. It should be \"type = passive\" or \"type = active\"" << std::endl;
                is_error = true;
              }
            }
          }
          cells_active.shrink_to_fit();
          cells_passive.shrink_to_fit();
          environments_active.shrink_to_fit();
          environments_passive.shrink_to_fit();
          std::vector<mp::iniConfiguration::section_info> section_to_check = {
                  {"GLOBAL", global_check}
          };
          for (auto & it : cells_active){
            section_to_check.push_back(std::make_pair(it, active_cell_check));
          }
          for (auto & it : cells_passive){
            section_to_check.push_back(std::make_pair(it, passive_cell_check));
          }
          for (auto & it : environments_active){
            section_to_check.push_back(std::make_pair(it, active_env_check));
          }
          for (auto & it : environments_passive){
            section_to_check.push_back(std::make_pair(it, passive_env_check));
          }

          is_error = is_error || c.hasFieldsWithType(section_to_check, errors);
          return std::make_pair(is_error, errors.str());
        };



    }
}