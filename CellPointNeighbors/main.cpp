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

#include <list>

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

	// Find all cells connected to point 0
	vtkIdType cellId = 0;

	vtkSmartPointer<vtkIdList> cellPointIds =
		vtkSmartPointer<vtkIdList>::New();
	triangleFilter->GetOutput()->GetCellPoints(cellId, cellPointIds);

	// neighbor cells may be listed multiple times
	// use std::set instead of std::list if you want a unique list of neighbors
	std::list<vtkIdType> neighbors;

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

	std::cout << "Point neighbor ids are: " << std::endl;

	for (std::list<vtkIdType>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
	{
		std::cout << " " << *it1;
	}
	std::cout << std::endl;


	vtkSmartPointer<vtkDataSetMapper> sphereMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	sphereMapper->SetInputData(data);
	vtkSmartPointer<vtkActor> sphereActor =
		vtkSmartPointer<vtkActor>::New();
	sphereActor->SetMapper(sphereMapper);
	sphereActor->GetProperty()->SetEdgeColor(0, 0, 0);
	sphereActor->GetProperty()->EdgeVisibilityOn();

	vtkSmartPointer<vtkDataSetMapper> mainCellMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();

	vtkSmartPointer<vtkDataSetMapper> neighborCellsMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();

	// Create a dataset with the cell of interest
	{
		vtkSmartPointer<vtkIdTypeArray> ids =
			vtkSmartPointer<vtkIdTypeArray>::New();
		ids->SetNumberOfComponents(1);
		ids->InsertNextValue(cellId);

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

		mainCellMapper->SetInputConnection(extractSelection->GetOutputPort());

	}

	vtkSmartPointer<vtkActor> mainCellActor =
		vtkSmartPointer<vtkActor>::New();
	mainCellActor->SetMapper(mainCellMapper);
	mainCellActor->GetProperty()->SetColor(1, 0, 0);

	// Create a dataset with the neighbor cells
	{
		vtkSmartPointer<vtkIdTypeArray> ids =
			vtkSmartPointer<vtkIdTypeArray>::New();
		ids->SetNumberOfComponents(1);
		for (std::list<vtkIdType>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
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
	neighborCellsActor->GetProperty()->SetColor(1, 1, 0);

	//Create a renderer, render window, and interactor
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	//Add the actors to the scene
	renderer->AddActor(sphereActor);
	renderer->AddActor(mainCellActor);
	renderer->AddActor(neighborCellsActor);
	renderer->SetBackground(.3, .2, .1); // Background color dark red

	//Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}