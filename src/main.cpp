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

void usage(char *p_name) {
    cerr << "Usage: " << p_name
         << " vtk_mesh new_mesh_dx new_mesh_dy new_mesh_dz unit_conversion_rate"
            " output_filename  [extra_config.ini]"
         << endl;
    exit(1);
}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

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

void print_progress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\rProgress - %3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}


int main(int argc, char *argv[]) {

    if (!(argc == 7 || argc == 8)) {
        usage(argv[0]);
    }

    if (!valid_extension(argv[1])) {
        cerr << "Invalid input file! Please convert your mesh to VTU format (.vtu extension)!" << endl;
        exit(EXIT_FAILURE);
    }

    std::vector<dataConfig> *configs = nullptr;

    std::cout << "Start reading mesh from " << argv[1] << std::endl;

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    reader->SetFileName(argv[1]);
    reader->Update();

    double dx = strtod(argv[2], nullptr);
    double dy = strtod(argv[3], nullptr);
    double dz = strtod(argv[4], nullptr);

    double conversion_rate = strtod(argv[5], nullptr);

    ofstream converted_mesh;
    converted_mesh.open (argv[6]);

    double bounds[6];

    reader->GetOutput()->GetBounds(bounds);

    double mesh_min_x = bounds[0], mesh_min_y = bounds[2], mesh_min_z = bounds[4];

    vtkSmartPointer<vtkTransform> t =
            vtkSmartPointer<vtkTransform>::New();

    t->Translate(-mesh_min_x, -mesh_min_y, -mesh_min_z);

    vtkSmartPointer<vtkTransformFilter> tf =
            vtkSmartPointer<vtkTransformFilter>::New();


    cout << "Translating mesh to 0, 0, 0" << endl;
    tf->SetInputData(reader->GetOutput());
    tf->SetTransform(t);
    tf->Update();

    vtkDataSet *data = vtkDataSet::SafeDownCast(tf->GetOutput());

    if(argc == 8) {
        configs = parse_ini_config(argv[7], data);
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

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
            vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputData(tf->GetOutput());
    surfaceFilter->Update();

    vtkPolyData *meshData = surfaceFilter->GetOutput();
    vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints =
            vtkSmartPointer<vtkSelectEnclosedPoints>::New();

    selectEnclosedPoints->Initialize(meshData);
    selectEnclosedPoints->SetTolerance(0.0);

    meshData->GetBounds(bounds);

    cout << "End reading mesh from " << argv[1] << endl;

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

                            if(config.location == CELL) {
                                dataLocation = cellId;
                            }
                            else {
                                dataLocation = pointId;
                            }

                            if(config.nComponents == 1) {
                                switch (config.type) {
                                    case INT:
                                        dataValue = (double) config.intArray->GetValue(dataLocation);
                                        break;
                                    case FLOAT:
                                        dataValue = (double) config.floatArray->GetValue(dataLocation);
                                        break;

                                    case DOUBLE:
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
                centerz += dx;
            }
            centery += dy;
        }
        centerx += dz;
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
