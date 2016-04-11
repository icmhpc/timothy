//
// Created by Grzegorz Bokota on 11.02.16.
//

#include <python2.7/Python.h>
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
            std::vector<mp::iniConfiguration::field_info> global_check = {
                    {"size_x",              mp::NUMBER_FIELD, true},
                    {"size_y",              mp::NUMBER_FIELD, true},
                    {"size_z",              mp::NUMBER_FIELD, true},
                    {"max_number_of_cells", mp::NUMBER_FIELD, true},
            };
            std::vector<mp::iniConfiguration::field_info> active_cell_check = {
                    {"type",     mp::STRING_FIELD, true},
                    {"g1_phase", mp::FLOAT_FIELD,  true},
                    {"s_phase",  mp::FLOAT_FIELD,  true},
                    {"g2_phase", mp::FLOAT_FIELD,  true},
                    {"m_phase",  mp::FLOAT_FIELD,  true},
                    {"pre_function", mp::STRING_FIELD, false},
                    {"pre_function_data", mp::STRUCT_FIELD, false},
                    {"post_function", mp::STRING_FIELD, false},
                    {"post_function_data", mp::STRUCT_FIELD, false}
            };

            std::vector<mp::iniConfiguration::field_info> passive_cell_check;
            std::vector<mp::iniConfiguration::field_info> active_env_check = {
                    {"dc",               mp::FLOAT_FIELD, true},
                    {"bc",               mp::FLOAT_FIELD, true},
                    {"ic_mean",          mp::FLOAT_FIELD, true},
                    {"ic_var",           mp::FLOAT_FIELD, true},
                    {"lambda",           mp::FLOAT_FIELD, true},
                    {"critical_level_1", mp::FLOAT_FIELD, true},
                    {"critical_level_2", mp::FLOAT_FIELD, true}
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