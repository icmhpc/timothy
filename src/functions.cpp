//
// Created by czaki on 14.04.16.
//

#include "functions.h"
#include <math.h>


double ballEnv(doubleVector3d pos, void * d){
  struct local_data {
      double x;
      double y;
      double z;
      double intensity;
  };
  local_data * data =(local_data *) d;
  double dist = sqrt(
          (pos.x - data->x) * (pos.x - data->x) +
                  (pos.y - data->y) * (pos.y - data->y) +
                  (pos.z - data->z) * (pos.z - data->z));
  return data->intensity/(dist * dist * dist);
}


functionalEnv getFunEnvFunction(std::string s) {
  static std::map<std::string, functionalEnv > functions = {{"bellEnv", ballEnv}};
  return functions.at(s);

}

cellFunction getCellFunction(std::string s) {
  static std::map<std::string, cellFunction > cell_functions;
  return cell_functions.at(s);
}
