#include <vtkPolyData.h>
#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkConeSource.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		std::cout << "Required parameters: filename.stl" << std::endl;
		return EXIT_FAILURE;
	}

	std::string filename = argv[1];

	vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->Update();
	vtkSmartPointer<vtkSphereSource> sphereSource2 =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource2->SetCenter(0, 0, 1);
	sphereSource2->Update();

	vtkSmartPointer<vtkSTLWriter> stlWriter =
		vtkSmartPointer<vtkSTLWriter>::New();
	stlWriter->SetFileName(filename.c_str());
	stlWriter->SetInputConnection(sphereSource->GetOutputPort());
	stlWriter->SetInputConnection(sphereSource2->GetOutputPort());
	stlWriter->Write();

	// Read and display for verification
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(sphereSource->GetOutputPort());

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
	renderer->SetBackground(.3, .6, .3); // Background color green

	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}