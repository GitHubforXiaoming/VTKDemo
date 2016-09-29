#include <vtkVersion.h>
#include <vtkSmartPointer.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>

int main(int, char *[])
{
	vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->Update();

	std::cout << "There are " << sphereSource->GetOutput()->GetNumberOfPoints()
		<< " input points." << std::endl;
	std::cout << "There are " << sphereSource->GetOutput()->GetNumberOfCells()
		<< " input cells." << std::endl;

	vtkSmartPointer<vtkIdTypeArray> ids =
		vtkSmartPointer<vtkIdTypeArray>::New();
	ids->SetNumberOfComponents(1);

	// Specify that we want to extract cells 10 through 19
	for (unsigned int i = 10; i < 20; i++)
	{
		ids->InsertNextValue(i);
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
	extractSelection->SetInputConnection(0, sphereSource->GetOutputPort());
#if VTK_MAJOR_VERSION <= 5
	extractSelection->SetInput(1, selection);
#else
	extractSelection->SetInputData(1, selection);
#endif
	extractSelection->Update();

	// In selection
	vtkSmartPointer<vtkUnstructuredGrid> selected =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected->ShallowCopy(extractSelection->GetOutput());

	std::cout << "There are " << selected->GetNumberOfPoints()
		<< " points in the selection." << std::endl;
	std::cout << "There are " << selected->GetNumberOfCells()
		<< " cells in the selection." << std::endl;

	// Get points that are NOT in the selection
	selectionNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1); //invert the selection
	extractSelection->Update();

	vtkSmartPointer<vtkUnstructuredGrid> notSelected =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	notSelected->ShallowCopy(extractSelection->GetOutput());

	std::cout << "There are " << notSelected->GetNumberOfPoints()
		<< " points NOT in the selection." << std::endl;
	std::cout << "There are " << notSelected->GetNumberOfCells()
		<< " cells NOT in the selection." << std::endl;


	vtkSmartPointer<vtkProperty> backfaces =
		vtkSmartPointer<vtkProperty>::New();
	backfaces->SetColor(1, 0, 0);

	vtkSmartPointer<vtkDataSetMapper> inputMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	inputMapper->SetInputConnection(sphereSource->GetOutputPort());
	vtkSmartPointer<vtkActor> inputActor =
		vtkSmartPointer<vtkActor>::New();
	inputActor->SetMapper(inputMapper);
	inputActor->SetBackfaceProperty(backfaces);

	vtkSmartPointer<vtkDataSetMapper> selectedMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	selectedMapper->SetInputConnection(selected->GetProducerPort());
#else
	selectedMapper->SetInputData(selected);
#endif
	vtkSmartPointer<vtkActor> selectedActor =
		vtkSmartPointer<vtkActor>::New();
	selectedActor->SetMapper(selectedMapper);
	selectedActor->SetBackfaceProperty(backfaces);

	vtkSmartPointer<vtkDataSetMapper> notSelectedMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	notSelectedMapper->SetInputConnection(notSelected->GetProducerPort());
#else
	notSelectedMapper->SetInputData(notSelected);
#endif
	vtkSmartPointer<vtkActor> notSelectedActor =
		vtkSmartPointer<vtkActor>::New();
	notSelectedActor->SetMapper(notSelectedMapper);
	notSelectedActor->SetBackfaceProperty(backfaces);

	// There will be one render window
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(900, 300);

	// And one interactor
	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	// Define viewport ranges
	// (xmin, ymin, xmax, ymax)
	double leftViewport[4] = { 0.0, 0.0, 0.5, 1.0 };
	double centerViewport[4] = { 0.33, 0.0, .66, 1.0 };
	double rightViewport[4] = { 0.66, 0.0, 1.0, 1.0 };

	// Create a camera for all renderers
	vtkSmartPointer<vtkCamera> camera =
		vtkSmartPointer<vtkCamera>::New();

	// Setup the renderers
	vtkSmartPointer<vtkRenderer> leftRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderWindow->AddRenderer(leftRenderer);
	leftRenderer->SetViewport(leftViewport);
	leftRenderer->SetBackground(.6, .5, .4);
	leftRenderer->SetActiveCamera(camera);

	vtkSmartPointer<vtkRenderer> centerRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderWindow->AddRenderer(centerRenderer);
	centerRenderer->SetViewport(centerViewport);
	centerRenderer->SetBackground(.3, .1, .4);
	centerRenderer->SetActiveCamera(camera);

	vtkSmartPointer<vtkRenderer> rightRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderWindow->AddRenderer(rightRenderer);
	rightRenderer->SetViewport(rightViewport);
	rightRenderer->SetBackground(.4, .5, .6);
	rightRenderer->SetActiveCamera(camera);

	leftRenderer->AddActor(inputActor);
	centerRenderer->AddActor(selectedActor);
	rightRenderer->AddActor(notSelectedActor);

	leftRenderer->ResetCamera();

	renderWindow->Render();
	interactor->Start();

	return EXIT_SUCCESS;
}