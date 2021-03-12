#ifndef __CONFIG_H
#define __CONFIG_H 

#include<vector>
#include<string>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>

enum dataType {
    DOUBLE_TYPE,
    INT_TYPE,
    FLOAT_TYPE,
};

enum dataLocation {
    CELL_LOCATION,
    POINT_LOCATION
};

struct dataConfig {
	std::string dataName;
    dataType type;
    dataLocation location;
    int nComponents;

    union {
        vtkDataArray   *dataArray;
        vtkIntArray    *intArray;
        vtkFloatArray  *floatArray;
        vtkDoubleArray *doubleArray;
    };

};

std::vector<dataConfig> * parse_ini_config(const std::string &ini_filename, vtkDataSet *data);

#endif /* __CONFIG_H */
