#include <vtkVersion.h>

#include <vtkSphereSource.h>
#include <vtkStripper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkHardwareSelector.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractSelection.h>
#include <vtkSelection.h>
#include <vtkProperty.h>

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static KeyPressInteractorStyle* New();
	vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

	virtual void OnKeyPress()
	{
		// Get the keypress
		std::string key = this->Interactor->GetKeySym();

		// "s" for "s"elect
		if (key.compare("x") == 0)
		{
			vtkSmartPointer<vtkHardwareSelector> selector =
				vtkSmartPointer<vtkHardwareSelector>::New();
			selector->SetRenderer(this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
			int* temp = this->Interactor->GetRenderWindow()->GetSize();
			unsigned int windowSize[4];
			windowSize[0] = temp[2];
			windowSize[1] = temp[3];
			windowSize[2] = temp[0];
			windowSize[3] = temp[1];
			/*
			for(unsigned int i = 0; i < 4; i++)
			{
			windowSize[i] = temp[i];
			}
			*/
			selector->SetArea(windowSize);
			selector->SetFieldAssociation(vtkDataObject::FIELD_ASSOCIATION_CELLS);
			vtkSelection* selection = selector->Select();
			std::cout << "Selection has " << selection->GetNumberOfNodes() << " nodes." << std::endl;

			vtkSmartPointer<vtkExtractSelection> extractSelection =
				vtkSmartPointer<vtkExtractSelection>::New();
#if VTK_MAJOR_VERSION <= 5
			extractSelection->SetInput(0, this->Data);
			extractSelection->SetInput(1, selection);
#else
			extractSelection->SetInputData(0, this->Data);
			extractSelection->SetInputData(1, selection);
#endif
			extractSelection->Update();

			vtkSmartPointer<vtkDataSetMapper> mapper =
				vtkSmartPointer<vtkDataSetMapper>::New();
			mapper->SetInputConnection(extractSelection->GetOutputPort());

			vtkSmartPointer<vtkActor> actor =
				vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(1, 0, 0);
			this->Renderer->AddActor(actor);

		}

		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
	}

	vtkPolyData* Data;
	vtkRenderer* Renderer;

};
vtkStandardNewMacro(KeyPressInteractorStyle);

int main(int, char *[])
{
	//convert triangles to triangle strips
	vtkSmartPointer<vtkSphereSource> sphereSource1 =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource1->Update();

	std::cout << "Number of cells before stripping: " << sphereSource1->GetOutput()->GetNumberOfCells() << std::endl;

	vtkSmartPointer<vtkStripper> stripper =
		vtkSmartPointer<vtkStripper>::New();
	stripper->SetInputConnection(sphereSource1->GetOutputPort());
	stripper->Update();

	std::cout << "Number of cells after stripping: " << stripper->GetOutput()->GetNumberOfCells() << std::endl;

	//convert vtkUnstructuredGrid to vtkPolyData(1)
	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid1 =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
		vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	surfaceFilter->SetInput(unstructuredGrid);
#else
	surfaceFilter->SetInputData(unstructuredGrid1);
#endif
	surfaceFilter->Update();

	vtkPolyData* polydata1 = surfaceFilter->GetOutput();

	std::cout << "(1)Output has " << polydata1->GetNumberOfPoints() << " points." << std::endl;
	//convert vtkUnstructuredGrid to vtkPolyData(2)
	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkGeometryFilter> geometryFilter =
		vtkSmartPointer<vtkGeometryFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	geometryFilter->SetInput(unstructuredGrid);
#else
	geometryFilter->SetInputData(unstructuredGrid2);
#endif
	geometryFilter->Update();

	vtkPolyData* polydata2 = geometryFilter->GetOutput();

	std::cout << "(2)Output has " << polydata2->GetNumberOfPoints() << " points." << std::endl;

	// Create a sphere
	vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(5.0);
	sphereSource->Update();

	// Create a mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(sphereSource->GetOutputPort());

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// Create a renderer, render window, and interactor
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetMultiSamples(0); // Turn off anti-aliasing
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<KeyPressInteractorStyle> style =
		vtkSmartPointer<KeyPressInteractorStyle>::New();
	style->Renderer = renderer;
	renderWindowInteractor->SetInteractorStyle(style);
	style->SetCurrentRenderer(renderer);
	style->Data = sphereSource->GetOutput();

	// Add the actor to the scene
	renderer->AddActor(actor);
	renderer->SetBackground(.3, .6, .3); // Background color green

	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}