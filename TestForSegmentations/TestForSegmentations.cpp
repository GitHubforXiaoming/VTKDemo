#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkTriangleFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkIdTypeArray.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkProperty.h>
#include <vtkExtractEdges.h>

#include <list>

void GetNeighbor(vtkSmartPointer<vtkTriangleFilter> triangleFilter, std::list<vtkIdType>& neighbors, vtkIdType cellId)
{
	vtkSmartPointer<vtkIdList> cellPointIds =
		vtkSmartPointer<vtkIdList>::New();
	triangleFilter->GetOutput()->GetCellPoints(cellId, cellPointIds);

	/*For each vertice of the cell, we calculate which cells uses that point.
	So if we make this, for each vertice, we have all the neighbors.
	In the case we use ''cellPointIds'' as a parameter of ''GeteCellNeighbors'',
	we will obtain an empty set. Because the only cell that is using that set of points
	is the current one. That is why we have to make each vertice at time.*/

	for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkIdList> idList =
			vtkSmartPointer<vtkIdList>::New();
		idList->InsertNextId(cellPointIds->GetId(i));

		//get the neighbors of the cell
		vtkSmartPointer<vtkIdList> neighborCellIds =
			vtkSmartPointer<vtkIdList>::New();

		triangleFilter->GetOutput()->GetCellNeighbors(cellId, idList, neighborCellIds);

		for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++)
		{
			neighbors.push_back(neighborCellIds->GetId(j));
		}

	}
}

void GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, std::list<vtkIdType>& neighbors, int id)
{
	vtkSmartPointer<vtkIdList> connectedVertices =
		vtkSmartPointer<vtkIdList>::New();

	//get all cells that vertex 'id' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(id, cellIdList);


	/*cout << "Vertex " << id << " is used in cells ";
	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
		cout << cellIdList->GetId(i) << ", ";
	}
	cout << endl;*/


	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
		//cout << "cell id " << i << " : " << cellIdList->GetId(i) << endl;

		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

		//cout << "End points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;

		if (pointIdList->GetId(0) != id)
		{
			//cout << "Connected to " << pointIdList->GetId(0) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(0));
		}
		else
		{
			//cout << "Connected to " << pointIdList->GetId(1) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(1));
		}
	}

	for (vtkIdType j = 0; j < connectedVertices->GetNumberOfIds(); j++)
	{
		neighbors.push_back(connectedVertices->GetId(j));
	}
}

int main(int argc, char *argv[])
{
	vtkSmartPointer<vtkPolyData> data;
	if (argc > 1)
	{
		vtkSmartPointer<vtkSTLReader> reader =
			vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(argv[1]);
		reader->Update();

		data = reader->GetOutput();
	}
	else
	{
		// Create a sphere
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetCenter(0.0, 0.0, 0.0);
		sphereSource->SetRadius(0.5);
		sphereSource->Update();

		data = sphereSource->GetOutput();
	}

	vtkSmartPointer<vtkTriangleFilter> triangleFilter =
		vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(data);
	triangleFilter->Update();
	vtkSmartPointer<vtkExtractEdges> extractEdges =
		vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInputConnection(triangleFilter->GetOutputPort());
	extractEdges->Update();

	vtkSmartPointer<vtkPolyData> mesh = extractEdges->GetOutput();

	// Find all cells connected to point 0

	for (vtkIdType i = 0; i < data->GetNumberOfPoints(); i++)
	{
		std::list<vtkIdType> neighbors;
		GetConnectedVertices(mesh, neighbors, i);
		std::cout << "---------------------------" << endl;
		std::cout << i << ": Point neighbor ids are: " << std::endl;

		for (std::list<vtkIdType>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
		{
			std::cout << " " << *it1;
		}
		std::cout << std::endl;
		std::cout << "size is " << neighbors.size() << endl;
		std::cout << "---------------------------" << endl;
	}

	system("pause");
	return EXIT_SUCCESS;
}