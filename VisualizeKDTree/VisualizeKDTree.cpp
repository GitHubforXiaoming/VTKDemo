#include <vtkVersion.h>
#include <vtkPointSource.h>
#include <vtkKdTreePointLocator.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkSliderWidget.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkWidgetEvent.h>
#include <vtkCallbackCommand.h>
#include <vtkWidgetEventTranslator.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkProperty.h>
#include <vtkMath.h>

class vtkSliderCallback : public vtkCommand
{
public:
	static vtkSliderCallback *New()
	{
		return new vtkSliderCallback;
	}
	vtkSliderCallback() :KdTree(0), Level(0), PolyData(0), Renderer(0){}

	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkSliderWidget *sliderWidget =
			reinterpret_cast<vtkSliderWidget*>(caller);
		this->Level = vtkMath::Round(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
		this->KdTree->GenerateRepresentation(this->Level, this->PolyData);
		this->Renderer->Render();
	}

	vtkKdTreePointLocator* KdTree;
	int Level;
	vtkPolyData* PolyData;
	vtkRenderer* Renderer;
};

int main(int, char *[])
{
	// Create a point cloud
	vtkSmartPointer<vtkPointSource> pointSource =
		vtkSmartPointer<vtkPointSource>::New();
	pointSource->SetRadius(4);
	pointSource->SetNumberOfPoints(1000);
	pointSource->Update();

	vtkSmartPointer<vtkPolyDataMapper> pointsMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	pointsMapper->SetInputConnection(pointSource->GetOutputPort());

	vtkSmartPointer<vtkActor> pointsActor =
		vtkSmartPointer<vtkActor>::New();
	pointsActor->SetMapper(pointsMapper);
	pointsActor->GetProperty()->SetInterpolationToFlat();

	int maxLevel = 20;
	// Create the tree
	vtkSmartPointer<vtkKdTreePointLocator> kdTree =
		vtkSmartPointer<vtkKdTreePointLocator>::New();
	kdTree->SetDataSet(pointSource->GetOutput());
	kdTree->AutomaticOff();
	kdTree->SetMaxLevel(maxLevel);
	kdTree->BuildLocator();

	// Initialize the representation
	vtkSmartPointer<vtkPolyData> polydata =
		vtkSmartPointer<vtkPolyData>::New();
	kdTree->GenerateRepresentation(0, polydata);

	vtkSmartPointer<vtkPolyDataMapper> octreeMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	octreeMapper->SetInputConnection(polydata->GetProducerPort());
#else
	octreeMapper->SetInputData(polydata);
#endif

	vtkSmartPointer<vtkActor> octreeActor =
		vtkSmartPointer<vtkActor>::New();
	octreeActor->SetMapper(octreeMapper);
	octreeActor->GetProperty()->SetInterpolationToFlat();
	octreeActor->GetProperty()->SetRepresentationToWireframe();

	// A renderer and render window
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	// An interactor
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actors to the scene
	renderer->AddActor(pointsActor);
	renderer->AddActor(octreeActor);

	// Render an image (lights and cameras are created automatically)
	renderWindow->Render();

	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep =
		vtkSmartPointer<vtkSliderRepresentation2D>::New();
	sliderRep->SetMinimumValue(0);
	sliderRep->SetMaximumValue(kdTree->GetLevel());
	sliderRep->SetValue(0);
	sliderRep->SetTitleText("Level");
	sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	sliderRep->GetPoint1Coordinate()->SetValue(.2, .2);
	sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	sliderRep->GetPoint2Coordinate()->SetValue(.8, .2);
	sliderRep->SetSliderLength(0.075);
	sliderRep->SetSliderWidth(0.05);
	sliderRep->SetEndCapLength(0.05);

	vtkSmartPointer<vtkSliderWidget> sliderWidget =
		vtkSmartPointer<vtkSliderWidget>::New();
	sliderWidget->SetInteractor(renderWindowInteractor);
	sliderWidget->SetRepresentation(sliderRep);
	sliderWidget->SetAnimationModeToAnimate();
	sliderWidget->EnabledOn();

	vtkSmartPointer<vtkSliderCallback> callback =
		vtkSmartPointer<vtkSliderCallback>::New();
	callback->KdTree = kdTree;
	callback->PolyData = polydata;
	callback->Renderer = renderer;

	sliderWidget->AddObserver(vtkCommand::InteractionEvent, callback);

	renderWindowInteractor->Initialize();
	renderWindow->Render();

	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}