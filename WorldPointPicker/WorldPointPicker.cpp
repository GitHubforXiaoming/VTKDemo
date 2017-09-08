#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>

// Define interaction style
class MouseInteractorStylePP : public vtkInteractorStyleTrackballCamera
{
public:
	static MouseInteractorStylePP* New();
	vtkTypeMacro(MouseInteractorStylePP, vtkInteractorStyleTrackballCamera);

	virtual void OnLeftButtonDown()
	{
		int* pos = this->Interactor->GetEventPosition();
		std::cout << "Picking pixel: " << pos[0] << " " << pos[1] << std::endl;
		vtkSmartPointer<vtkPointPicker> pointPicker = vtkSmartPointer<vtkPointPicker>::New();
		pointPicker->SetTolerance(0.00005);
		pointPicker->Pick(pos[0], pos[1], 0, this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());

		std::cout << "Cell id is: " << pointPicker->GetPointId() << std::endl;

		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);
		std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;
		// Forward events
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	}

};

vtkStandardNewMacro(MouseInteractorStylePP);

int main(int, char *[])
{
	vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->Update();

	vtkSmartPointer<vtkPointPicker> pointPicker =
		vtkSmartPointer<vtkPointPicker>::New();

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
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetPicker(pointPicker);
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<MouseInteractorStylePP> style =
		vtkSmartPointer<MouseInteractorStylePP>::New();
	renderWindowInteractor->SetInteractorStyle(style);

	// Add the actor to the scene
	renderer->AddActor(actor);
	renderer->SetBackground(1, 1, 1); // Background color white

	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();
	return EXIT_SUCCESS;
}