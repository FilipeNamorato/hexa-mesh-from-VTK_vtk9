#include <vtkCellArray.h>
#include <vtkCellCenters.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkCenterOfMass.h>
#include <vtkCubeSource.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractSelection.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkMergePoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

void usage(char *p_name) {
	cout << "Usage: " << p_name << " vtk_mesh new_mesh_dx new_mesh_dy new_mesh_dz unit_conversion_rate output_filename" << endl; 
	exit(1);
}

int main(int argc, char *argv[]) {

	if(argc != 7) usage(argv[0]);

    cout << "Start reading mesh from " << argv[1] << endl;

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    reader->SetFileName(argv[1]);
    reader->Update();

	double dx = strtod(argv[2], nullptr);
	double dy = strtod(argv[3], nullptr);
	double dz = strtod(argv[4], nullptr);

	double conversion_rate = strtod(argv[5], nullptr);

	FILE *converted_mesh = fopen(argv[6],"w");

	std::string fibers_filename(argv[6]);
	fibers_filename.append(".fibers");

	FILE *fibers_file = fopen(fibers_filename.c_str(), "w");

    vtkDataSet *data = vtkDataSet::SafeDownCast(reader->GetOutput());

    vtkSmartPointer<vtkPointLocator> pointLoc = vtkSmartPointer<vtkPointLocator>::New();

    pointLoc->SetDataSet(data);
    pointLoc->AutomaticOn();
    pointLoc->BuildLocator();

	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();

    cellLocator->SetDataSet(data);
    cellLocator->BuildLocator();

    int numCells = data->GetNumberOfCells();

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputData(reader->GetOutput());
    surfaceFilter->Update();

    vtkPolyData *meshData = surfaceFilter->GetOutput();
    vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
    selectEnclosedPoints->Initialize(meshData);
    selectEnclosedPoints->SetTolerance(0.0);

    double bounds[6];

    meshData->GetBounds(bounds);
    
	cout << "End reading mesh from " << argv[1] << endl;

    double mesh_min_x = bounds[0], mesh_min_y = bounds[2], mesh_min_z = bounds[4];
    double mesh_max_x = bounds[1], mesh_max_y = bounds[3], mesh_max_z = bounds[5];

    double min_x = 0, min_y = 0, min_z = 0;
    double max_x = mesh_max_x + dx, max_y = mesh_max_y + dy, max_z = mesh_max_z + dz;

    cout << "Mesh X from " << mesh_min_x  << " to " << mesh_max_x << endl;
    cout << "Mesh Y from " << mesh_min_y  << " to " << mesh_max_y << endl;
    cout << "Mesh Z from " << mesh_min_z  << " to " << mesh_max_z << endl;

    cout << "X from " << min_x  << " to " << max_x << endl;
    cout << "Y from " << min_y  << " to " << max_y << endl;
    cout << "Z from " << min_z  << " to " << max_z << endl;

    vtkSmartPointer<vtkMergePoints> pointLocator = vtkSmartPointer<vtkMergePoints>::New();

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> points = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> initialPoints = vtkSmartPointer<vtkPoints>::New();


	/*
	 * THIS NEED TO BE CHANGED IN ORDER EXTRACT THE ARRAYS FROM THE VTK FILE. THE DATA ARRAYS CAN HAVE DIFFERENT NAMES
	 * AND CAN BE ASSOCIATED TO A POINT OR A CELL. 
	 * IF THE ARRAY IS ASSOCIATED WITH A POINT YOU NEED TO USE GetPointData(), USE GetCellData() OTHERWISE.
	 * VTK SUPPORTS TUPLE ARRAYS, IF THAT IS THE CASE, YOU HAVE TO USE A vtkDataArray.
	 * vtkFloatArray, vtkIntArray or vtkDoubleArray ARE USED WHEN THE COMPONENTS OF THE ARRAY ARE SCALARS AND NOT TUPLES.
	 * THE PARAMETER FOR THE GetAbstractArray FUNCTION IS THE NAME OF THE ARRAY IN THE VTK FILE. 
	 * YOU CAN CHECK THE NAMES OF THE ARRAYS AND IF THEY ARE CELL OR POINT DATA USING PARAVIEW.
	 */

	//THE EXAMPLE BELOW IS TO CONVERT A SPECIFIC FILE WITH THESE VECTORS ASSOCIATED WITH IT
	/*
	vtkDataArray  *fiber  = vtkDataArray::SafeDownCast(data->GetPointData()->GetAbstractArray("FIBRE"));
	vtkDataArray  *sheet  = vtkDataArray::SafeDownCast(data->GetPointData()->GetAbstractArray("SHEET"));
	vtkDataArray  *normal = vtkDataArray::SafeDownCast(data->GetPointData()->GetAbstractArray("NORMA"));
	vtkFloatArray *trans  = vtkFloatArray::SafeDownCast(data->GetPointData()->GetAbstractArray("TRANS"));
	vtkFloatArray *hcmre  = vtkFloatArray::SafeDownCast(data->GetCellData()->GetAbstractArray("HCMRE"));
	*/

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
    vtkIdType pointId;

    double centerx = min_x + halfl, centery, centerz;

    double center_point[3];
	
	int total_points = total_points_x * total_points_y * total_points_z;
	int count = 1;

    for (int i = 0; i < total_points_x; ++i) {
		
		centery = min_y + halfl;
		centerz = min_z + halfl;

        for (int j = 0; j < total_points_y; ++j) {

			centerz = min_z + halfl;

            for (int z = 0; z < total_points_z; ++z) {

				double percentage = ((double)count / (double)total_points)*100.0;

				if(	fabsf((int)percentage - percentage) < 0.0001 ) {
					printf("adding point %d of %d - %.2lf%% completed\n", count, total_points, percentage);	
				}
				count++;

                center_point[0] = centerx;
                center_point[1] = centery;
                center_point[2] = centerz;

                if (selectEnclosedPoints->IsInsideSurface(center_point)) {
				
					static int add_counter = 0;

                    //Find the closest points to TestPoint
                    double closestPoint[3];  //the coordinates of the closest point will be returned here
                    double closestPointDist2;//the squared distance to the closest point will be returned here
                    vtkIdType cellId;        //the cell id of the cell containing the closest point will be returned here
                    vtkIdType pointId;
                    int subId;               //this is rarely used (in triangle strips only, I believe)

                    cellLocator->FindClosestPoint(center_point, closestPoint, cellId, subId, closestPointDist2);
                    pointId = pointLoc->FindClosestPoint(center_point);


					//THE EXAMPLE BELOW GETS VALUES FROM A POINT DATA ARRAY AND A CELL DATE ARRAY RESPECTIVELY
					//float transValue = trans->GetValue(pointId);
					//float hcmreValue = hcmre->GetValue(cellId);

					//THIS PRINTF IS USED TO CUSTOMIZE HOW THE NEW MESH WILL BE SAVED.
					//THE FIRST SIX PARAMETERS ARE GEOMETRIC INFORMATION AND SHOULD NOT BE CHANGED
					//THE REST OF THE PARAMETERS CAN BE CHANGED AND THEY ARE USUALLY DATA 
					//THAT COME FROM THE DATA ARRAYS EXTRACTED FROM THE VTK FILE
                    
					//fprintf(converted_mesh, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lld,%lld\n",
                    //       (centerx) * conversion_rate,
                    //       (centery) * conversion_rate,
                    //       (centerz) * conversion_rate,
                    //       halfl * conversion_rate,
                    //       halfl * conversion_rate,
                    //       halfl * conversion_rate,
					//	   transValue,
					//	   hcmreValue,
					//	   cellId,
					//	   pointId);
					
					//THIS printf WILL SAVE THE GEOMETRIC DATA AND ALSO THE CELL ID OF THE CLOSEST CELL
					//AND THE POINT ID OF THE CLOSEST POINT TO THE CENTRE OF THE VOLUME THAT IS BEING CREATED
					//FOR THE NEW MESH
					fprintf(converted_mesh, "%lf,%lf,%lf,%lf,%lf,%lf,%lld,%lld\n",
                           (centerx) * conversion_rate,
                           (centery) * conversion_rate,
                           (centerz) * conversion_rate,
                           halfl * conversion_rate,
                           halfl * conversion_rate,
                           halfl * conversion_rate,
						   cellId,
						   pointId);					

					//THE EXAMPLE BELOW IS USED TO GET THE FIBRE, SHEET AND NORMAL INFORMATION FROM AN EXAMPLE MESH
					//double *fiber_data  = fiber->GetTuple(pointId);
					//double *sheet_data  = sheet->GetTuple(pointId);
					//double *normal_data = normal->GetTuple(pointId);

					//HERE WE ARE SAVING THE CONDUCTIVITY DATA IN A SEPARATE FILE WITH A .fibers EXTENSION
					//fprintf(fibers_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", fiber_data[0], fiber_data[1], fiber_data[2], sheet_data[0], 
					//		sheet_data[1], sheet_data[2], normal_data[0], normal_data[1], normal_data[2]);

					//NOTING NEEDS TO BE CHANGE FROM HERE
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

	fclose(fibers_file);
	fclose(converted_mesh);

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
}
