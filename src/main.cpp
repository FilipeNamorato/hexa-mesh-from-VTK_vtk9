#include <map>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkCenterOfMass.h>
#include <vtkDataSet.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkMergePoints.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPolyData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

#define INI_IMPLEMENTATION
#include "ini.h"

#include "config.h"

#include "ProgramOptions.hxx"
#include <iostream>

void usage(char *p_name) {
    std::cerr << "Usage: " << p_name << " [options]" << std::endl << std::endl;
    std::cerr << "Options:" <<  std::endl;
    std::cerr << "--input_mesh | -i, file VTU file containing the mesh to be converted (required)" << std::endl;
    std::cerr << "--dx | -x, X discretization of the converted mesh (same unit as the original mesh, required)" << std::endl;
    std::cerr << "--dy | -y, Y discretization of the converted mesh (same unit as the original mesh, required)" << std::endl;
    std::cerr << "--dz | -z, Z discretization of the converted mesh (same unit as the original mesh, required)" << std::endl;
    std::cerr << "--conversion_rate | -r, rate that needs to be applied to convert the original mesh size unit to micrometers (optional, default: 1)" << std::endl;
    std::cerr << "--output_file | -o, name of the converted mesh file (optional, default: converted_mesh.alg)" << std::endl;
    std::cerr << "--config_file | -c, name of .ini file that has the configuration to extract the arrays data from the VTU file (optional, default: not used)" << std::endl;
    std::cerr << "--help | -h. Shows this help and exit" << std::endl;
    exit(EXIT_FAILURE);

}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void print_progress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\rProgress - %3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

bool valid_extension(const std::string& filename) {

    std::string::size_type idx;

    idx = filename.rfind('.');

    if (idx != std::string::npos) {
        std::string extension = filename.substr(idx + 1);
        return extension == "vtu";
    } else {
        return false;
    }
}

int main(int argc, char *argv[]) {

	po::parser parser;
    auto& input_mesh_opt = parser["input_mesh"]
		.abbreviation('i')            
        .type(po::string);            

	auto& desired_dx_opt = parser["dx"]  
		.abbreviation('x')            
        .type(po::f64);            

	auto& desired_dy_opt = parser["dy"]  
		.abbreviation('y')            
        .type(po::f64);

	auto& desired_dz_opt = parser["dz"]  
		.abbreviation('z')            
        .type(po::f64);            

	auto& conversion_rate_opt = parser["conversion_rate"]
		.abbreviation('r')            
        .type(po::f64)
        .fallback(1.0);

	auto& output_file_opt = parser["output_file"]  
		.abbreviation('o')            
        .type(po::string)
		.fallback("converted_mesh.alg");  

	auto& config_file_opt = parser["config_file"]  
		.abbreviation('c')            
        .type(po::string);

    auto& help_opt = parser["help"]
            .abbreviation('h');

    parser(argc, argv);               // parses the command line arguments

    if(help_opt.was_set()) {
        usage(argv[0]);
    }

	bool error = !input_mesh_opt.available() || !desired_dx_opt.available() || !desired_dy_opt.available() || !desired_dz_opt.available();

    if(error) {
        std::cout << "Wrong number of arguments!" << std::endl << std::endl;
		usage(argv[0]);
	}

	auto input_file = input_mesh_opt.get().string;

    if (!valid_extension(input_file)) {
        cerr << "Invalid input file! Please convert your mesh to VTU format (.vtu extension)!" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "================================================================================" << endl;
    cout << "Conversion configuration:" << endl;
    cout << "================================================================================" << endl;
    cout << "Input mesh file: " << input_file << endl;
    cout << "Target dx: " << desired_dx_opt.get().f64 << endl;
    cout << "Target dy: " << desired_dy_opt.get().f64 << endl;
    cout << "Target dz: " << desired_dz_opt.get().f64 << endl;
    cout << "Conversion rate: " << conversion_rate_opt.get().f64 << endl;
    cout << "Output file: " << output_file_opt.get().string << endl;
    if(config_file_opt.was_set()) {
        cout << "Configuration file: " << config_file_opt.get().string << endl;
    }
    else {
        cout << "Data arrays will not be extracted from the input mesh" << endl;
    }
    cout << "================================================================================" << endl << endl;

	//TODO: print usage and print the input configurations on the screen

    std::vector<dataConfig> *configs = nullptr;

    std::cout << "Start reading mesh from " << input_file << std::endl;

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    reader->SetFileName(input_file.c_str());
    reader->Update();

    double dx = desired_dx_opt.get().f64;
    double dy = desired_dy_opt.get().f64; 
    double dz = desired_dz_opt.get().f64;

    double conversion_rate = conversion_rate_opt.get().f64; 

    ofstream converted_mesh;
    converted_mesh.open (output_file_opt.get().string);

    double bounds[6];

    reader->GetOutput()->GetBounds(bounds);

    double mesh_min_x = bounds[0], mesh_min_y = bounds[2], mesh_min_z = bounds[4];

    vtkDataSet *data;
	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints =  vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	//Translating the original mesh
	vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
	t->Translate(-mesh_min_x, -mesh_min_y, -mesh_min_z);

	vtkSmartPointer<vtkTransformFilter> tf = vtkSmartPointer<vtkTransformFilter>::New();

	cout << "Translating mesh to 0, 0, 0" << endl;
	tf->SetInputData(reader->GetOutput());
	tf->SetTransform(t);
	tf->Update();
	
	//This is used to extract the arrays in the VTU mesh
	data = vtkDataSet::SafeDownCast(tf->GetOutput());

	
	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	surfaceFilter->SetInputData(tf->GetOutput());
	surfaceFilter->Update();
	{
		vtkPolyData *meshData = surfaceFilter->GetOutput();

		selectEnclosedPoints->Initialize(meshData);
		selectEnclosedPoints->SetTolerance(0.0);
		meshData->GetBounds(bounds);
	}

    if(config_file_opt.available()) {
        configs = parse_ini_config(config_file_opt.get().string, data);
    }

    vtkSmartPointer<vtkPointLocator> pointLoc =
            vtkSmartPointer<vtkPointLocator>::New();

    pointLoc->SetDataSet(data);
    pointLoc->AutomaticOn();
    pointLoc->BuildLocator();

    vtkSmartPointer<vtkCellLocator> cellLocator =
            vtkSmartPointer<vtkCellLocator>::New();

    cellLocator->SetDataSet(data);
    cellLocator->BuildLocator();

    cout << "End reading mesh" << endl;

    double mesh_max_x = bounds[1], mesh_max_y = bounds[3], mesh_max_z = bounds[5];

    double min_x = 0, min_y = 0, min_z = 0;
    double max_x = mesh_max_x + dx, max_y = mesh_max_y + dy,
            max_z = mesh_max_z + dz;


    mesh_min_x = bounds[0], mesh_min_y = bounds[2], mesh_min_z = bounds[4];
    mesh_max_x = bounds[1], mesh_max_y = bounds[3], mesh_max_z = bounds[5];

    min_x = 0, min_y = 0, min_z = 0;
    max_x = mesh_max_x + dx, max_y = mesh_max_y + dy,
            max_z = mesh_max_z + dz;

    cout << "Mesh X from " << mesh_min_x << " to " << mesh_max_x << endl;
    cout << "Mesh Y from " << mesh_min_y << " to " << mesh_max_y << endl;
    cout << "Mesh Z from " << mesh_min_z << " to " << mesh_max_z << endl;

    cout << "X from " << min_x << " to " << max_x << endl;
    cout << "Y from " << min_y << " to " << max_y << endl;
    cout << "Z from " << min_z << " to " << max_z << endl;

    vtkSmartPointer<vtkMergePoints> pointLocator =
            vtkSmartPointer<vtkMergePoints>::New();

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> points = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> initialPoints = vtkSmartPointer<vtkPoints>::New();

    double aux_p[3];

    aux_p[0] = 0.0;
    aux_p[1] = 0.0;
    aux_p[2] = 0.0;
    initialPoints->InsertNextPoint(aux_p);

    aux_p[0] = max_x;
    aux_p[1] = 0.0;
    aux_p[2] = 0.0;
    initialPoints->InsertNextPoint(aux_p);

    aux_p[0] = 0.0;
    aux_p[1] = max_y;
    aux_p[2] = 0.0;
    initialPoints->InsertNextPoint(aux_p);

    aux_p[0] = max_x;
    aux_p[1] = max_y;
    aux_p[2] = 0.0;
    initialPoints->InsertNextPoint(aux_p);

    //////////////
    aux_p[0] = 0.0;
    aux_p[1] = 0.0;
    aux_p[2] = max_z;
    initialPoints->InsertNextPoint(aux_p);

    aux_p[0] = max_x;
    aux_p[1] = 0.0;
    aux_p[2] = max_z;
    initialPoints->InsertNextPoint(aux_p);

    aux_p[0] = 0.0;
    aux_p[1] = max_y;
    aux_p[2] = max_z;
    initialPoints->InsertNextPoint(aux_p);

    aux_p[0] = max_x;
    aux_p[1] = max_y;
    aux_p[2] = max_z;

    initialPoints->InsertNextPoint(aux_p);

    points->SetPoints(initialPoints);

    pointLocator->InitPointInsertion(points->GetPoints(), points->GetBounds());

    int total_points_x = (max_x - min_x) / dx;
    int total_points_y = (max_y - min_y) / dy;
    int total_points_z = (max_z - min_z) / dz;

    double aux1[3];
    double aux2[3];
    double aux3[3];
    double aux4[3];
    double aux5[3];
    double aux6[3];
    double aux7[3];
    double aux8[3];

    double halfl = dx / 2.0;
    double centerx = min_x + halfl, centery, centerz;

    double center_point[3];

    int total_points = total_points_x * total_points_y * total_points_z;
    int count = 1;

    for (int i = 0; i < total_points_x; ++i) {

        centery = min_y + halfl;

        for (int j = 0; j < total_points_y; ++j) {

            centerz = min_z + halfl;

            for (int z = 0; z < total_points_z; ++z) {

	            double percentage = ((double) count / (double) total_points);
				print_progress(percentage);

                count++;

                center_point[0] = centerx;
                center_point[1] = centery;
                center_point[2] = centerz;

                if (selectEnclosedPoints->IsInsideSurface(center_point)) {

                    static int add_counter = 0;

                    double closestPoint[3];  // the coordinates of the closest point will
                    // be returned here
                    double closestPointDist2;// the squared distance to the closest point
                    // will be returned here
                    vtkIdType cellId;        // the cell id of the cell containing the closest
                    // point will be returned here
                    vtkIdType pointId;
                    int subId;// this is rarely used (in triangle strips only, I believe)

                    cellLocator->FindClosestPoint(center_point, closestPoint, cellId, subId, closestPointDist2);

                    pointId = pointLoc->FindClosestPoint(center_point);

                    converted_mesh << centerx*conversion_rate << "," << centery*conversion_rate << "," << centerz *conversion_rate << ","
                                   << halfl * conversion_rate << "," << halfl * conversion_rate << "," << halfl * conversion_rate;

                    if(!configs || configs->empty()) {
                        converted_mesh << endl;
                    }
                    else {
                        for(auto & config : *configs) {

                            vtkIdType dataLocation;
                            double dataValue;

                            if(config.location == CELL_LOCATION) {
                                dataLocation = cellId;
                            }
                            else {
                                dataLocation = pointId;
                            }

                            if(config.nComponents == 1) {
                                switch (config.type) {
                                    case INT_TYPE:
                                        dataValue = (double) config.intArray->GetValue(dataLocation);
                                        break;
                                    case FLOAT_TYPE:
                                        dataValue = (double) config.floatArray->GetValue(dataLocation);
                                        break;

                                    case DOUBLE_TYPE:
                                        dataValue = (double) config.doubleArray->GetValue(dataLocation);
                                        break;
                                }

                                converted_mesh << "," << dataValue;
                            }
                            else {

                                double *components = config.dataArray->GetTuple(dataLocation);
                                for(int k = 0; k < config.nComponents; k++) {
                                    converted_mesh << "," << components[k];
                                }

                            }
                        }
                        converted_mesh << endl;
                    }

                    aux1[0] = (centerx - halfl);
                    aux1[1] = (centery - halfl);
                    aux1[2] = (centerz - halfl);
                    pointLocator->InsertUniquePoint(aux1, pointId);

                    aux2[0] = (centerx + halfl);
                    aux2[1] = (centery - halfl);
                    aux2[2] = (centerz - halfl);
                    pointLocator->InsertUniquePoint(aux2, pointId);

                    aux3[0] = (centerx + halfl);
                    aux3[1] = (centery + halfl);
                    aux3[2] = (centerz - halfl);
                    pointLocator->InsertUniquePoint(aux3, pointId);

                    aux4[0] = (centerx - halfl);
                    aux4[1] = (centery + halfl);
                    aux4[2] = (centerz - halfl);
                    pointLocator->InsertUniquePoint(aux4, pointId);

                    aux5[0] = (centerx - halfl);
                    aux5[1] = (centery - halfl);
                    aux5[2] = (centerz + halfl);
                    pointLocator->InsertUniquePoint(aux5, pointId);

                    aux6[0] = (centerx + halfl);
                    aux6[1] = (centery - halfl);
                    aux6[2] = (centerz + halfl);
                    pointLocator->InsertUniquePoint(aux6, pointId);

                    aux7[0] = (centerx + halfl);
                    aux7[1] = (centery + halfl);
                    aux7[2] = (centerz + halfl);
                    pointLocator->InsertUniquePoint(aux7, pointId);

                    aux8[0] = (centerx - halfl);
                    aux8[1] = (centery + halfl);
                    aux8[2] = (centerz + halfl);
                    pointLocator->InsertUniquePoint(aux8, pointId);

                    cells->InsertNextCell(8);
                    cells->InsertCellPoint(pointLocator->IsInsertedPoint(aux1));
                    cells->InsertCellPoint(pointLocator->IsInsertedPoint(aux2));
                    cells->InsertCellPoint(pointLocator->IsInsertedPoint(aux3));
                    cells->InsertCellPoint(pointLocator->IsInsertedPoint(aux4));
                    cells->InsertCellPoint(pointLocator->IsInsertedPoint(aux5));
                    cells->InsertCellPoint(pointLocator->IsInsertedPoint(aux6));
                    cells->InsertCellPoint(pointLocator->IsInsertedPoint(aux7));
                    cells->InsertCellPoint(pointLocator->IsInsertedPoint(aux8));
                }
                centerz += dz;
            }
            centery += dy;
        }
        centerx += dx;
    }

    converted_mesh.close();

    selectEnclosedPoints->Complete();
    vtkSmartPointer<vtkUnstructuredGrid> ug =
            vtkSmartPointer<vtkUnstructuredGrid>::New();

    ug->SetPoints(points->GetPoints());

    ug->SetCells(VTK_HEXAHEDRON, cells);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    writer->SetInputData(ug);
    writer->SetFileName("converted_mesh.vtu");
    writer->Write();

    cout << endl;
}
