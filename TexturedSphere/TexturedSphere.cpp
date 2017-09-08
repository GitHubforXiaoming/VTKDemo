#include <vtkSmartPointer.h>

#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkTransformTextureCoords.h>
#include <vtkTexture.h>
#include <vtkTextureMapToSphere.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageReader.h>
#include <vtkTexturedSphereSource.h>
#include <vtkSTLReader.h>
#include <vtkSphereSource.h>
#include <vtkJPEGReader.h>
#include <vtkTextureMapToCylinder.h>
#include <vtkTextureMapToSphere.h>
#include <vtkTextureMapToPlane.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkProperty.h>
#include <vtkPlaneSource.h>

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0]
			<< " texture(.png)"
			<< " [translate]" << std::endl;
		return EXIT_FAILURE;
	}
	vtkSmartPointer<vtkPolyData> data;
	if (argc > 2)
	{
		vtkSmartPointer<vtkSTLReader> reader =
			vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(argv[2]);
		reader->Update();
		data = reader->GetOutput();
	}
	else
	{
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetThetaResolution(30);
		sphereSource->SetPhiResolution(15);
		sphereSource->Update();
		data = sphereSource->GetOutput(); 
	}

	double translate[3];
	translate[0] = 4.0;
	translate[1] = 4.0;
	translate[2] = 1.0;
	
	// Create a sphere with texture coordinates
	// Create a plane
	vtkSmartPointer<vtkPlaneSource> plane =
		vtkSmartPointer<vtkPlaneSource>::New();
	plane->SetCenter(0.0, 0.0, 0.0);
	plane->SetNormal(0.0, 0.0, 1.0);

	// Read texture file
	vtkSmartPointer<vtkJPEGReader> reader =
		vtkSmartPointer<vtkJPEGReader>::New();
		reader->SetFileName(argv[1]);

	// Create texture
	vtkSmartPointer<vtkTexture> texture =
		vtkSmartPointer<vtkTexture>::New();
	texture->SetInputConnection(reader->GetOutputPort());
	texture->InterpolateOn();
	texture->Update();

	vtkSmartPointer<vtkTextureMapToPlane> textureMap = vtkSmartPointer<vtkTextureMapToPlane>::New();
	textureMap->SetInputData(data);
	textureMap->Update();

	vtkSmartPointer<vtkTransformTextureCoords> ttc = vtkSmartPointer<vtkTransformTextureCoords>::New();
	ttc->SetInputConnection(textureMap->GetOutputPort());
	ttc->SetScale(translate);

	vtkSmartPointer<vtkPointData> pointData = textureMap->GetOutput()->GetPointData();
	vtkSmartPointer<vtkDataArray> coords = pointData->GetTCoords();
	std::cout << coords->GetNumberOfComponents() << endl;
	for (vtkIdType i = 0; i < coords->GetNumberOfTuples(); i++)
	{
		std::cout << i << ":(" << coords->GetTuple(i)[0] << "," << coords->GetTuple(i)[1] << ")" << endl;
	}
	 
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(textureMap->GetOutputPort());
	//mapper->SetInputData(data);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->SetTexture(texture);

	//vtkSmartPointer<vtkImageData> imageData = actor->GetTexture()->GetImageDataInput(0);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->SetBackground(.1, .2, .3);

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renderWindow);

	renderWindow->Render();
	renWinInteractor->Start();

	return EXIT_SUCCESS;
}