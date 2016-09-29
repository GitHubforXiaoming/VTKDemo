#include <algorithm>

#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkPlane.h>
#include <vtkPlaneSource.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkCellCenters.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkUnsignedCharArray.h>

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

using namespace std;

double planeNormal[3] = { 0.0, 0.0, 1.0 };
double planeCenter[3] = { 0.0, 0.0, 1.0 };

//��ȡÿ��cell������
void getNormals(vtkSmartPointer<vtkPolyData> polyData)
{
	vtkSmartPointer<vtkDoubleArray> normalsArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	normalsArray->SetNumberOfComponents(3);
	normalsArray->SetNumberOfTuples(polyData->GetNumberOfPoints());
	polyData->GetCellData()->SetNormals(normalsArray);
	vtkSmartPointer<vtkDoubleArray> cellNormalsRetrieve =
		vtkDoubleArray::SafeDownCast(polyData->GetCellData()->GetNormals());
	
	if (cellNormalsRetrieve)
	{
		cout << "There are " << cellNormalsRetrieve->GetNumberOfTuples() << " cells normals" << endl;
		for (vtkIdType i = 0; i < cellNormalsRetrieve->GetNumberOfTuples(); i++)
		{
			double n[3];
			cellNormalsRetrieve->GetTuple(i, n);
			
			cout << "Cell normal" << i << ": (" << n[0] << ", " << n[1] << ", " << n[2] << ")" << endl; 
 		}
	}
	else
	{
		cout << "No cell normals" << endl;
	}
}

//��ȡÿ��cell�����ĵ�
vtkSmartPointer<vtkCellCenters> getCenterPoints(vtkSmartPointer<vtkPolyData> polyData)
{
	int pointsNum = vtkMath::Random(0, 10);
	vtkSmartPointer<vtkCellCenters> centersFilter =
		vtkSmartPointer<vtkCellCenters>::New();
	centersFilter->SetInputData(polyData);
	centersFilter->VertexCellsOn();
	centersFilter->Update();
	cout << "һ��" << centersFilter->GetOutput()->GetNumberOfPoints() << "cells" << endl;
	for (vtkIdType i = 0; i < centersFilter->GetOutput()->GetNumberOfPoints(); i++)
	{
		double point[3];
		centersFilter->GetOutput()->GetPoint(i, point);

		//cout << "Point " << i << ":(" << point[0] << ", " << point[1] << ", " << point[2] << ")" << endl;
	}
	return centersFilter;
}

//Ϊÿ��cellȾɫ
vtkSmartPointer<vtkUnsignedCharArray> dyeForCells(int numbersOfCells)
{
	// Create cell data
	vtkMath::RandomSeed(8775070); // for reproducibility
	vtkSmartPointer<vtkUnsignedCharArray> cellData =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	cellData->SetNumberOfComponents(3);
	cellData->SetNumberOfTuples(numbersOfCells);
	for (int i = 0; i < numbersOfCells; i++)
	{
		float rgb[3];
		rgb[0] = vtkMath::Random(64, 255);
		rgb[1] = vtkMath::Random(64, 255);
		rgb[2] = vtkMath::Random(64, 255);
		cellData->InsertTuple(i, rgb);
	}
	return cellData;
}

//�ҳ�һ����Ŀ����С����
double* getMinDistance(vtkSmartPointer<vtkFloatArray> distances)
{
	double minDistance[20];
	int length = distances->GetNumberOfTuples();
	vector<double> data;
	for (int i = 0; i < length; i++)
	{
		data.push_back(distances->GetTuple(i)[0]);
	}
	sort(data.begin(), data.end());
	cout << "��С���루ǰ20��";
	for (int n = 0; n < 20; n++)
	{
		minDistance[n] = data[n];
		if (!(n % 5))
		{
			cout << endl;
		}
		cout << n << "��" << minDistance[n] << "\t";
	}
	cout << endl;
	return minDistance;
}

int main(int, char *[])
{
	// Create a plane
	vtkSmartPointer<vtkPlaneSource> planeSource =
		vtkSmartPointer<vtkPlaneSource>::New();
	planeSource->SetCenter(planeCenter);
	planeSource->SetNormal(planeNormal);
	planeSource->Update();

	// Create a sphere
	vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(0.5);
	sphereSource->Update();

	vtkSmartPointer<vtkCellCenters> cellCenter =
		getCenterPoints(sphereSource->GetOutput());
	vtkSmartPointer<vtkPolyData> linesPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> lines =
		vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	// compute distances to each point
	vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance =
		vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicitPolyDataDistance->SetInput(planeSource->GetOutput());
	vtkSmartPointer<vtkFloatArray> distances =
		vtkSmartPointer<vtkFloatArray>::New();
	distances->SetNumberOfComponents(1);
	distances->SetName("Distances");
	for (vtkIdType i = 0; i < cellCenter->GetOutput()->GetNumberOfPoints(); i++)
	{
		double pointOnSphere[3];
		double pointOnPlane[3];
		cellCenter->GetOutput()->GetPoint(i, pointOnSphere);
		pointOnPlane[0] = pointOnSphere[0];
		pointOnPlane[1] = pointOnSphere[1];
		pointOnPlane[2] = 1.0;
		//��������صĵ㼯�浽һ��������
		points->InsertNextPoint(pointOnSphere);
		points->InsertNextPoint(pointOnPlane);
		//�������Բ�ϵ�������Ƭ�е㵽ƽ��ľ���
		float distance = implicitPolyDataDistance->EvaluateFunction(pointOnSphere);
		if (distance < 0)
		{
			distance = -distance;
		}
		distances->InsertNextValue(distance);
		if (!(i % 5))
		{
			cout << endl;
		}
		cout << "Distance:" << distance << "\t"; 
		
	}
	cout << endl;

	getMinDistance(distances);
	
	linesPolyData->SetPoints(points);
	//����ֱ��
	for (vtkIdType j = 0; j < linesPolyData->GetNumberOfPoints();)
	{
		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		//cout << "point1:" << points->GetPoint(j)[0] << "," << points->GetPoint(j)[1] << "," << points->GetPoint(j)[2] << endl;
		//cout << "point2:" << points->GetPoint(j + 1)[0] << "," << points->GetPoint(j + 1)[1] << "," << points->GetPoint(j + 1)[2] << endl;
		line->GetPointIds()->SetId(0, j);
		line->GetPointIds()->SetId(1, j + 1);
		j += 2;
		lines->InsertNextCell(line);
	}
	linesPolyData->SetLines(lines);
	linesPolyData->GetCellData()->SetScalars(dyeForCells(cellCenter->GetOutput()->GetNumberOfPoints()));
	
	// Create a mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> planeMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	vtkSmartPointer<vtkPolyDataMapper> sphereMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	vtkSmartPointer<vtkPolyDataMapper> centerMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	vtkSmartPointer<vtkPolyDataMapper> vertexMapper = 
		vtkSmartPointer<vtkPolyDataMapper>::New();

#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(plane);
#else
	planeMapper->SetInputData(planeSource->GetOutput());
	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
	centerMapper->SetInputConnection(cellCenter->GetOutputPort());
	vertexMapper->SetInputData(linesPolyData);
#endif

	vtkSmartPointer<vtkActor> planeActor =
		vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkActor> sphereActor =
		vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkActor> centerActor =
		vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkActor> vertexActor =
		vtkSmartPointer<vtkActor>::New();
	planeActor->SetMapper(planeMapper);
	sphereActor->SetMapper(sphereMapper);
	centerActor->SetMapper(centerMapper);
	vertexActor->SetMapper(vertexMapper);

	// Create a renderer, render window and interactor
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetGlobalWarningDisplay(0);
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actors to the scene
	renderer->AddActor(planeActor);
	renderer->AddActor(sphereActor);
	renderer->AddActor(centerActor);
	renderer->AddActor(vertexActor);
	renderer->SetBackground(.2, .3, .4); // Background color dark blue

	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}