#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkXMLPolyDataWriter.h>

#include <vector>
#include <string>

void FindAllData(vtkPolyData* polydata);

int main(int argc, char *argv[])
{
	vtkSmartPointer<vtkPolyData> polydata =
		vtkSmartPointer<vtkPolyData>::New();

	if (argc < 2)
	{
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->Update();
		vtkSmartPointer<vtkXMLPolyDataWriter> writer =
			vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetFileName("test.vtp");
		writer->SetInputConnection(sphereSource->GetOutputPort());
		writer->Write();

		polydata = sphereSource->GetOutput();
	}
	else
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader =
			vtkSmartPointer<vtkXMLPolyDataReader>::New();
		reader->SetFileName(argv[1]);
		reader->Update();
		polydata = reader->GetOutput();
	}

	FindAllData(polydata);
	system("pause");
	return EXIT_SUCCESS;
}

void FindAllData(vtkPolyData* polydata)
{
	std::cout << "Normals: " << polydata->GetPointData()->GetNormals() << std::endl;

	vtkIdType numberOfPointArrays = polydata->GetPointData()->GetNumberOfArrays();
	std::cout << "Number of PointData arrays: " << numberOfPointArrays << std::endl;

	vtkIdType numberOfCellArrays = polydata->GetCellData()->GetNumberOfArrays();
	std::cout << "Number of CellData arrays: " << numberOfCellArrays << std::endl;

	std::cout << "Type table/key: " << std::endl;;
	//more values can be found in <VTK_DIR>/Common/vtkSetGet.h
	std::cout << VTK_UNSIGNED_CHAR << " unsigned char" << std::endl;
	std::cout << VTK_UNSIGNED_INT << " unsigned int" << std::endl;
	std::cout << VTK_FLOAT << " float" << std::endl;
	std::cout << VTK_DOUBLE << " double" << std::endl;

	for (vtkIdType i = 0; i < numberOfPointArrays; i++)
	{
		// The following two lines are equivalent
		//arrayNames.push_back(polydata->GetPointData()->GetArray(i)->GetName());
		//arrayNames.push_back(polydata->GetPointData()->GetArrayName(i));
		int dataTypeID = polydata->GetPointData()->GetArray(i)->GetDataType();
		std::cout << "Array " << i << ": " << polydata->GetPointData()->GetArrayName(i)
			<< " (type: " << dataTypeID << ")" << std::endl;
	}

	for (vtkIdType i = 0; i < numberOfCellArrays; i++)
	{
		// The following two lines are equivalent
		//polydata->GetPointData()->GetArray(i)->GetName();
		//polydata->GetPointData()->GetArrayName(i);
		int dataTypeID = polydata->GetCellData()->GetArray(i)->GetDataType();
		std::cout << "Array " << i << ": " << polydata->GetCellData()->GetArrayName(i)
			<< " (type: " << dataTypeID << ")" << std::endl;
	}
}