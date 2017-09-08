#include <vtkSmartPointer.h>

#include <vtkImageCanvasSource2D.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTexture.h>

#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>

int main(int argc, char *argv[])
{
	// Create an image
	vtkSmartPointer<vtkImageCanvasSource2D> imageSource =
		vtkSmartPointer<vtkImageCanvasSource2D>::New();
	imageSource->SetScalarTypeToUnsignedChar();
	imageSource->SetExtent(0, 20, 0, 20, 0, 0);
	imageSource->SetNumberOfScalarComponents(3);
	imageSource->SetDrawColor(127, 255, 100);
	imageSource->FillBox(0, 20, 0, 20);
	imageSource->SetDrawColor(20, 20, 20);
	imageSource->DrawSegment(0, 0, 19, 19);
	imageSource->DrawSegment(19, 0, 0, 19);
	imageSource->Update();

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
		vtkSmartPointer<vtkPlaneSource> sphereSource =
			vtkSmartPointer<vtkPlaneSource>::New();
		/*sphereSource->SetThetaResolution(30);
		sphereSource->SetPhiResolution(15);*/
		sphereSource->Update();
		data = sphereSource->GetOutput();
	}


	// Create a plane

	// Apply the texture
	vtkSmartPointer<vtkTexture> texture =
		vtkSmartPointer<vtkTexture>::New();
	texture->SetInputConnection(imageSource->GetOutputPort());

	vtkSmartPointer<vtkPolyDataMapper> planeMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	planeMapper->SetInputData(data);

	vtkSmartPointer<vtkActor> texturedPlane =
		vtkSmartPointer<vtkActor>::New();
	texturedPlane->SetMapper(planeMapper);
	texturedPlane->SetTexture(texture);

	// Visualize the textured plane
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(texturedPlane);
	renderer->SetBackground(1, 1, 1); // Background color white
	renderer->ResetCamera();

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderWindow->Render();

	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}