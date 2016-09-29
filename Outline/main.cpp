#include <vtkVersion.h>
#include <vtkSphereSource.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOutlineFilter.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>

int main(int, char *[])
{

	//Triangle Area
	//setup points (geometry)
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	points->InsertNextPoint(0.0, 0.0, 0.0);
	points->InsertNextPoint(1.0, 0.0, 0.0);
	points->InsertNextPoint(1.0, 1.0, 0.0);
	points->InsertNextPoint(0.0, 1.0, 0.0);

	//create a triangle on the three points in the polydata
	vtkSmartPointer<vtkTriangle> triangle1 =
		vtkSmartPointer<vtkTriangle>::New();
	triangle1->GetPointIds()->SetId(0, 0);
	triangle1->GetPointIds()->SetId(1, 1);
	triangle1->GetPointIds()->SetId(2, 2);

	vtkSmartPointer<vtkTriangle> triangle2 =
		vtkSmartPointer<vtkTriangle>::New();
	triangle2->GetPointIds()->SetId(0, 2);
	triangle2->GetPointIds()->SetId(1, 3);
	triangle2->GetPointIds()->SetId(2, 0);

	//add the triangles to the list of triangles
	vtkSmartPointer<vtkCellArray> triangles =
		vtkSmartPointer<vtkCellArray>::New();
	triangles->InsertNextCell(triangle1);
	triangles->InsertNextCell(triangle2);

	//create a polydata object
	vtkSmartPointer<vtkPolyData> polydata =
		vtkSmartPointer<vtkPolyData>::New();

	//add the geometry and topology to the polydata
	polydata->SetPoints(points);
	polydata->SetPolys(triangles);

	for (vtkIdType i = 0; i < polydata->GetNumberOfCells(); i++)
	{
		vtkCell* cell = polydata->GetCell(0);

		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(cell);
		double p0[3];
		double p1[3];
		double p2[3];
		triangle->GetPoints()->GetPoint(0, p0);
		std::cout << "p0: " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
		triangle->GetPoints()->GetPoint(1, p1);
		std::cout << "p1: " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
		triangle->GetPoints()->GetPoint(2, p2);
		std::cout << "p2: " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;

		double area = vtkTriangle::TriangleArea(p0, p1, p2);

		std::cout << "area of triangle " << i << ": " << area << std::endl;
	}

	// Create a sphere
	vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(5.0);
	sphereSource->Update();

	vtkPolyData* sphere = sphereSource->GetOutput();
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(sphere);
#else
	mapper->SetInputData(sphere);
#endif
	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// Create the outline
	vtkSmartPointer<vtkOutlineFilter> outline =
		vtkSmartPointer<vtkOutlineFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	outline->SetInput(sphere);
#else
	outline->SetInputData(sphere);
#endif
	vtkSmartPointer<vtkPolyDataMapper> outlineMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	outlineMapper->SetInputConnection(outline->GetOutputPort());
	vtkSmartPointer<vtkActor> outlineActor =
		vtkSmartPointer<vtkActor>::New();
	outlineActor->SetMapper(outlineMapper);
	outlineActor->GetProperty()->SetColor(0, 0, 0);

	// Setup the window
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actors to the scene
	renderer->AddActor(actor);
	renderer->AddActor(outlineActor);
	renderer->SetBackground(1, 1, 1); // Background color white

	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}