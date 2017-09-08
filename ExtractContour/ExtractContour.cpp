#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>

#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkExtractSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCellCenters.h>
#include <vtkDoubleArray.h>
#include <vtkTable.h>
#include <vtkPCAStatistics.h>
#include <vtkMath.h>
#include <vtkLine.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCubeAxesActor.h>
#include <vtkTextProperty.h>

#include <set>

using namespace std;
	

int main()
{
	int num = 0;
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName("1.stl");
	reader->Update();
	vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
	vtkSmartPointer<vtkPolyData> source = reader->GetOutput();

	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalGenerator->SetInputData(data);
	normalGenerator->ComputePointNormalsOff();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->Update();
	data = normalGenerator->GetOutput();

	/*
	 *获取三角面片的中点，利用这些点进行PCA降维，并计算出物体的中心点
	 *为了定位物体的中心，并分析三角面片的向量与某个方向的夹角
	 */
	vtkSmartPointer<vtkDataArray> normals = data->GetCellData()->GetNormals();
	vtkSmartPointer<vtkCellCenters> centers =
		vtkSmartPointer<vtkCellCenters>::New();
	vtkSmartPointer<vtkDoubleArray> centerPoints =
		vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	centers->SetInputData(data);
	centers->VertexCellsOn();
	centers->Update();
	for (vtkIdType i = 0; i < centers->GetOutput()->GetNumberOfPoints(); i++)
	{
		double point[3];
		centers->GetOutput()->GetPoint(i, point);
		point[2] = 0;
		points->InsertNextPoint(point);
	}
	
	vtkSmartPointer<vtkCubeAxesActor> cubeAxesActor =
		vtkSmartPointer<vtkCubeAxesActor>::New();
	cubeAxesActor->SetBounds(data->GetBounds());
	
	cubeAxesActor->GetTitleTextProperty(0)->SetColor(1.0, 0.0, 0.0);
	cubeAxesActor->GetLabelTextProperty(0)->SetColor(1.0, 0.0, 0.0);

	cubeAxesActor->GetTitleTextProperty(1)->SetColor(0.0, 1.0, 0.0);
	cubeAxesActor->GetLabelTextProperty(1)->SetColor(0.0, 1.0, 0.0);

	cubeAxesActor->GetTitleTextProperty(2)->SetColor(0.0, 0.0, 1.0);
	cubeAxesActor->GetLabelTextProperty(2)->SetColor(0.0, 0.0, 1.0);

	cubeAxesActor->DrawXGridlinesOn();
	cubeAxesActor->DrawYGridlinesOn();
	cubeAxesActor->DrawZGridlinesOn();
	cubeAxesActor->SetGridLineLocation(VTK_GRID_LINES_FURTHEST);

	cubeAxesActor->XAxisMinorTickVisibilityOff();
	cubeAxesActor->YAxisMinorTickVisibilityOff();
	cubeAxesActor->ZAxisMinorTickVisibilityOff();


	/*根据一定的阈值进行体数据的着色*/
	vtkSmartPointer<vtkIdTypeArray> ids =
		vtkSmartPointer<vtkIdTypeArray>::New();
	ids->SetNumberOfComponents(1);

	for (vtkIdType i = 0; i < normals->GetNumberOfTuples(); i++)
	{
		double normal[3];
		normals->GetTuple(i, normal);
		double cosine = (normal[2]) / sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
		if (abs(cosine) < 0.6)
		{
			ids->InsertNextValue(i);
			//cout << normal[0] << "," << normal[1] << "," << normal[2] << "\t" << cosine << endl;
			num++;
		}
	}
	cout << num << endl;
	vtkSmartPointer<vtkExtractSelection> extractSelection =
		vtkSmartPointer<vtkExtractSelection>::New();
	vtkSmartPointer<vtkSelectionNode> selectionNode =
		vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(ids);

	vtkSmartPointer<vtkSelection> selection =
		vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);

	extractSelection->SetInputData(0, source);
	extractSelection->SetInputData(1, selection);
	extractSelection->Update();

	// Setup actor and mapper
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	vtkSmartPointer<vtkDataSetMapper> contourMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputData(source);
	contourMapper->SetInputConnection(extractSelection->GetOutputPort());

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkActor> contourActor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	contourActor->SetMapper(contourMapper);
	contourActor->GetProperty()->SetColor(112 / 255.0, 24 / 255.0, 111 / 255.0);

	// Setup render window, renderer, and interactor
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetGlobalWarningDisplay(0);
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	cubeAxesActor->SetCamera(renderer->GetActiveCamera());
	renderer->AddActor(actor);
	renderer->AddActor(contourActor);
	renderer->AddActor(cubeAxesActor);
	renderer->SetBackground(.1, .2, .3);
	renderer->GetActiveCamera()->Azimuth(30);
	renderer->GetActiveCamera()->Elevation(30);
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}