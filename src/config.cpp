#include "config.h"
#include "ini.h"
#include <iostream>

#include <vtkCellData.h>
#include <vtkPointData.h>

std::vector<dataConfig> * parse_ini_config(const std::string &ini_filename, vtkDataSet * data) {

    std::ifstream t(ini_filename);
    std::string str((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

    ini_t *ini = ini_load(str.c_str(), nullptr);

    int num_sections = ini_section_count(ini);

	auto configs = new std::vector<dataConfig>();

	for (int i = 1; i < num_sections; i++) {

		int num_properties = ini_property_count(ini, i);

		auto section_name = std::string(ini_section_name(ini, i));

		dataConfig config;
		config.dataName = section_name;

		for(int j = 0; j < num_properties; j++) {

			auto property_name = std::string(ini_property_name(ini, i, j));
			auto property_value = std::string(ini_property_value(ini, i, j));

            if(property_name == "type") {
                if(property_value == "int") {
                    config.type = INT_TYPE;
                }
                else if(property_value == "double") {
                    config.type = DOUBLE_TYPE;
                }
                else if(property_value == "float") {
                    config.type = FLOAT_TYPE;
                }
                else {
                    std::cerr << "Invalid property value '" << property_value << "' for property type, on section [" << section_name << "]" << std::endl;
                    std::cerr << "Valid values are: int, float or double" << std::endl;
                    exit(EXIT_FAILURE);
                }
            } else if (property_name == "location") {
				if(property_value == "cell") {
					config.location = CELL_LOCATION;
				}
				else if(property_value == "point") {
					config.location = POINT_LOCATION;
				}
				else {
					std::cerr << "Invalid property value '" << property_value << "' for property location, on section [" << section_name << "]" << std::endl;
					std::cerr << "Valid values are: cell or point" << std::endl;
					exit(EXIT_FAILURE);
				}

			} else if(property_name == "n_components") {
				config.nComponents = std::stoi(property_value);

				if(config.nComponents <= 0) {
					std::cerr << "Invalid property value '" << property_value << "' for property n_components, on section [" << section_name << "]" << std::endl;
					std::cerr << "Valid values are: integers greater than 0" << std::endl;
					exit(EXIT_FAILURE);
				}
			} else {
				std::cerr << "Invalid property '" << property_name << "' on section [" << section_name << "]" << std::endl;
				exit(EXIT_FAILURE);
			}

		}

        //Loading the data arrays from the vtk file
        if(config.nComponents == 1) {
            switch (config.type) {
                case INT_TYPE:
                {
                    if (config.location == CELL_LOCATION) {
                        config.intArray = vtkIntArray::SafeDownCast(data->GetCellData()->GetAbstractArray(config.dataName.c_str()));
                    } else {
                        config.intArray = vtkIntArray::SafeDownCast(data->GetPointData()->GetAbstractArray(config.dataName.c_str()));
                    }
                }
                    break;

                case FLOAT_TYPE:
                {
                    if (config.location == CELL_LOCATION) {
                        config.floatArray = vtkFloatArray::SafeDownCast(data->GetCellData()->GetAbstractArray(config.dataName.c_str()));
                    } else {
                        config.floatArray = vtkFloatArray::SafeDownCast(data->GetPointData()->GetAbstractArray(config.dataName.c_str()));
                    }

                }
                    break;

                case DOUBLE_TYPE:
                {
                    if (config.location == CELL_LOCATION) {
                        config.doubleArray = vtkDoubleArray::SafeDownCast(data->GetCellData()->GetAbstractArray(config.dataName.c_str()));
                    } else {
                        config.doubleArray = vtkDoubleArray::SafeDownCast(data->GetPointData()->GetAbstractArray(config.dataName.c_str()));
                    }
                }
                    break;

            }
        }
        else {
            if (config.location == CELL_LOCATION) {
                config.dataArray = vtkDataArray::SafeDownCast(data->GetCellData()->GetAbstractArray(config.dataName.c_str()));
            } else {
                config.dataArray = vtkDataArray::SafeDownCast(data->GetPointData()->GetAbstractArray(config.dataName.c_str()));
            }
        }


		configs->push_back(config);
	}

	ini_destroy(ini);
    return configs;
}


