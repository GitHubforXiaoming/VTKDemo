#include <vtkSmartPointer.h>
#include <vtkVersion.h>

#include <vtkWarpVector.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkDataSetAttributes.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

int main(int argc, char *argv[])
{
	double scale = 1.0;
	vtkSmartPointer<vtkPolyData> inputPolyData;
	// If a file is provided, use it, otherwise generate a sphere
	if (argc > 1)
	{
		vtkSmartPointer<vtkSTLReader> reader =
			vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(argv[1]);
		reader->Update();
		inputPolyData = reader->GetOutput();
		if (argc > 2)
		{
			scale = atof(argv[2]);
		}
	}
	else
	{
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetPhiResolution(15);
		sphereSource->SetThetaResolution(15);
		sphereSource->Update();
		inputPolyData = sphereSource->GetOutput();
	}

	vtkSmartPointer<vtkCleanPolyData> clean =
		vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	clean->SetInput(inputPolyData);
#else
	clean->SetInputData(inputPolyData);
#endif

	// Generate normals
	vtkSmartPointer<vtkPolyDataNormals> normals =
		vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputConnection(clean->GetOutputPort());
	normals->SplittingOff();

	// Warp using the normals
	vtkSmartPointer<vtkWarpVector> warp =
		vtkSmartPointer<vtkWarpVector>::New();
	warp->SetInputConnection(normals->GetOutputPort());
	warp->SetInputArrayToProcess(0, 0, 0,
		vtkDataObject::FIELD_ASSOCIATION_POINTS,
		vtkDataSetAttributes::NORMALS);
	warp->SetScaleFactor(scale);

	// Visualize the original and warped models
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(warp->GetOutputPort());
	mapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> warpedActor =
		vtkSmartPointer<vtkActor>::New();
	warpedActor->SetMapper(mapper);

	vtkSmartPointer<vtkPolyDataMapper> originalMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	originalMapper->SetInputConnection(normals->GetOutputPort());
	originalMapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> originalActor =
		vtkSmartPointer<vtkActor>::New();
	originalActor->SetMapper(originalMapper);
	originalActor->GetProperty()->SetInterpolationToFlat();

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(640, 480);

	// Create a camera for all renderers
	vtkSmartPointer<vtkCamera> camera =
		vtkSmartPointer<vtkCamera>::New();

	// Define viewport ranges
	// (xmin, ymin, xmax, ymax)
	double leftViewport[4] = { 0.0, 0.0, 0.5, 1.0 };
	double rightViewport[4] = { 0.5, 0.0, 1.0, 1.0 };

	// Setup both renderers
	vtkSmartPointer<vtkRenderer> leftRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	leftRenderer->SetViewport(leftViewport);
	leftRenderer->SetBackground(.6, .5, .4);
	leftRenderer->SetActiveCamera(camera);

	vtkSmartPointer<vtkRenderer> rightRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	rightRenderer->SetViewport(rightViewport);
	rightRenderer->SetBackground(.4, .5, .6);
	rightRenderer->SetActiveCamera(camera);

	leftRenderer->AddActor(originalActor);
	rightRenderer->AddActor(warpedActor);

	rightRenderer->ResetCamera();

	renderWindow->AddRenderer(rightRenderer);
	renderWindow->AddRenderer(leftRenderer);

	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	renderWindow->Render();
	interactor->Start();

	return EXIT_SUCCESS;
}