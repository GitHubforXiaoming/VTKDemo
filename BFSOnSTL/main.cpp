#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkTriangleFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkIdTypeArray.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkProperty.h>
#include <vtkSTLReader.h>
#include <vtkPolyData.h>

#include <list>
#include <queue>
#include <set>

void GetNeighborIds(vtkIdType cellId,vtkSmartPointer<vtkTriangleFilter> triangleFilter, std::set<vtkIdType>& neighbors)
{
	vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
	triangleFilter->GetOutput()->GetCellPoints(cellId, cellPointIds);

	for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkIdList> idList =
			vtkSmartPointer<vtkIdList>::New();

		//add one of the edge points
		idList->InsertNextId(cellPointIds->GetId(i));

		//add the other edge point
		if (i + 1 == cellPointIds->GetNumberOfIds())
		{
			idList->InsertNextId(cellPointIds->GetId(0));
		}
		else
		{
			idList->InsertNextId(cellPointIds->GetId(i + 1));
		}
		//get the neighbors of the cell
		vtkSmartPointer<vtkIdList> neighborCellIds =
			vtkSmartPointer<vtkIdList>::New();

		triangleFilter->GetOutput()->GetCellNeighbors(cellId, idList, neighborCellIds);

		for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++)
		{
			neighbors.insert(neighborCellIds->GetId(j));
		}
	}
}

int main(int argc, char *argv[])
{

	// Create a sphere
	vtkSmartPointer<vtkPolyData> data;
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	if (argc > 1)
	{
		reader->SetFileName("1.stl");
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


	std::string inputFilename = argv[1];

	reader->SetFileName(inputFilename.c_str());
	reader->Update();

	data = reader->GetOutput();
	int size = data->GetNumberOfCells();
	int* visited = new int[size];
	for (int i = 0; i < size; i++)
	{
		visited[i] = 0;
	}

	vtkSmartPointer<vtkTriangleFilter> triangleFilter =
		vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(data);
	triangleFilter->Update();

	int n = 1;
	vtkIdType cellId = 0;
	std::vector<vtkIdType> cellIds;
	cellIds.push_back(cellId);

	std::queue<vtkIdType> q;
	visited[cellId] = 1;
	q.push(cellId);
	while (!q.empty())
	{
		cellId = q.front();
		q.pop();
		std::set<vtkIdType> neighbors;
		GetNeighborIds(cellId, triangleFilter, neighbors);
		for (std::set<vtkIdType>::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); neighbor++)
		{
			if (visited[*neighbor] == 0)
			{
				//std::cout << n++ << ":" << *neighbor << endl;;
				visited[*neighbor] = 1;
				q.push(*neighbor);
				cellIds.push_back(*neighbor);
			}
		}
	}
	/*vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(data);
	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetEdgeColor(0, 0, 0);
	actor->GetProperty()->EdgeVisibilityOn();*/

	vtkSmartPointer<vtkDataSetMapper> neighborCellsMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();

	// Create a dataset with the neighbor cells
	{
		vtkSmartPointer<vtkIdTypeArray> ids =
			vtkSmartPointer<vtkIdTypeArray>::New();
		ids->SetNumberOfComponents(1);
		for (std::vector<vtkIdType>::iterator it1 = cellIds.begin(); it1 != cellIds.end(); it1++)
		{
			ids->InsertNextValue(*it1);
		}

		vtkSmartPointer<vtkSelectionNode> selectionNode =
			vtkSmartPointer<vtkSelectionNode>::New();
		selectionNode->SetFieldType(vtkSelectionNode::CELL);
		selectionNode->SetContentType(vtkSelectionNode::INDICES);
		selectionNode->SetSelectionList(ids);

		vtkSmartPointer<vtkSelection> selection =
			vtkSmartPointer<vtkSelection>::New();
		selection->AddNode(selectionNode);

		vtkSmartPointer<vtkExtractSelection> extractSelection =
			vtkSmartPointer<vtkExtractSelection>::New();
		extractSelection->SetInputData(0, data);
#if VTK_MAJOR_VERSION <= 5
		extractSelection->SetInput(1, selection);
#else
		extractSelection->SetInputData(1, selection);
#endif
		extractSelection->Update();

		neighborCellsMapper->SetInputConnection(extractSelection->GetOutputPort());

	}

	vtkSmartPointer<vtkActor> neighborCellsActor =
		vtkSmartPointer<vtkActor>::New();
	neighborCellsActor->SetMapper(neighborCellsMapper);
	neighborCellsActor->GetProperty()->SetColor(0, 1, 0);

	// Visualize
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(data);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->AddActor(neighborCellsActor);
	renderer->SetBackground(.3, .6, .3); // Background color green

	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}