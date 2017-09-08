#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkPLYReader.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkSphereSource.h>
#include <vtkMaskPoints.h>
#include <vtkProperty.h>
#include <vtkRegularPolygonSource.h>


using namespace std;

//²âÊÔÎÄ¼þ£º../data/fran_cut.vtk

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << argv[0] << " *.vtk" << std::endl;
		return EXIT_FAILURE;
	}

	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(argv[1]);
	reader->Update();

	vtkSmartPointer<vtkPolyDataNormals> normFilter =
		vtkSmartPointer<vtkPolyDataNormals>::New();
	normFilter->SetInputData(reader->GetOutput());
	normFilter->SetComputePointNormals(1);
	normFilter->SetComputeCellNormals(0);
	normFilter->SetAutoOrientNormals(1);
	normFilter->SetSplitting(0);
	normFilter->Update();

	vtkSmartPointer<vtkDataArray> normals = normFilter->GetOutput()->GetPointData()->GetNormals();

	for (vtkIdType i = 0; i < normals->GetNumberOfTuples(); i++)
	{
		cout << "normal:(" << normals->GetTuple(i)[0] << "," << normals->GetTuple(i)[1] << "," << normals->GetTuple(i)[2] << endl;
	}

	vtkSmartPointer<vtkMaskPoints> mask =
		vtkSmartPointer<vtkMaskPoints>::New();
	mask->SetInputData(normFilter->GetOutput());
	//mask->SetMaximumNumberOfPoints(300);
	mask->SetOnRatio(5);
	//mask->RandomModeOn();

	vtkSmartPointer<vtkRegularPolygonSource> arrow =
		vtkSmartPointer<vtkRegularPolygonSource>::New();
	arrow->SetNumberOfSides(160);
	arrow->SetRadius(0.3);

	vtkSmartPointer<vtkGlyph3D> glyph =
		vtkSmartPointer<vtkGlyph3D>::New();
	glyph->SetSourceConnection(arrow->GetOutputPort());
	glyph->SetInputConnection(mask->GetOutputPort());
	glyph->SetVectorModeToUseNormal();
	glyph->SetScaleFactor(3);
	glyph->SetColorModeToColorByVector();
	glyph->SetScaleModeToScaleByVector();
	//glyph->OrientOn();
	glyph->Update();

	vtkSmartPointer<vtkPolyDataMapper> originMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	originMapper->SetInputData(reader->GetOutput());

	vtkSmartPointer<vtkActor> originActor =
		vtkSmartPointer<vtkActor>::New();
	originActor->SetMapper(originMapper);

	vtkSmartPointer<vtkPolyDataMapper> normedMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	normedMapper->SetInputData(normFilter->GetOutput());

	vtkSmartPointer<vtkActor> normedActor =
		vtkSmartPointer<vtkActor>::New();
	normedActor->SetMapper(normedMapper);

	vtkSmartPointer<vtkPolyDataMapper> glyphMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	glyphMapper->SetInputData(glyph->GetOutput());

	vtkSmartPointer<vtkActor> glyphActor =
		vtkSmartPointer<vtkActor>::New();
	glyphActor->SetMapper(glyphMapper);
	glyphActor->GetProperty()->SetColor(1., 0., 0.);

	vtkSmartPointer<vtkRenderer> glyphRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	glyphRenderer->AddActor(glyphActor);
	glyphRenderer->AddActor(normedActor);
	glyphRenderer->ResetCamera();
	glyphRenderer->SetBackground(1.0, 1.0, 1.0);

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(glyphRenderer);
	renderWindow->Render();
	renderWindow->SetWindowName("PolyDataNormal");

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}