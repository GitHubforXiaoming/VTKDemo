// The source
#include <vtkSuperquadricSource.h>
#include <vtkParametricTorus.h>
#include <vtkParametricRandomHills.h>
#include <vtkParametricFunctionSource.h>
#include <vtkCleanPolyData.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSphereSource.h>
#include <vtkPlaneSource.h>
#include <vtkPointSource.h>
#include <vtkRegularPolygonSource.h>
#include <vtkElevationFilter.h>
#include <vtkCellCenters.h>
#include <vtkExtractSelection.h>
#include <vtkUnstructuredGrid.h>
#include <vtkModifiedBSPTree.h>
// 
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkCellLocator.h>
#include <vtkExtractEdges.h>
// Curvatures
#include <vtkClipPolyData.h>
#include <vtkCurvatures.h>
#include <vtkImageData.h>
#include <vtkLine.h>
#include <vtkPlane.h>
#include <vtkImplicitBoolean.h>
#include <vtkImplicitPolyDataDistance.h>
// For annotating
#include <vtkVariantArray.h>
// Lookup table
#include <vtkColorSeries.h>
#include <vtkLookupTable.h>
// For glyphing
#include <vtkReverseSense.h>
#include <vtkMaskPoints.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkTransform.h>
// For contouring
#include <vtkBandedPolyDataContourFilter.h>
#include <vtkTextProperty.h>
// For Reader
#include <vtkBMPReader.h>
#include <vtkJPEGReader.h>
#include <vtkSTLReader.h>
// For Array
#include <vtkDoubleArray.h>
// Mappers, actors, renderers etc.
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageViewer2.h>
#include <vtkImageMapper3D.h>
#include <vtkTextMapper.h>
#include <vtkAxesActor.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkImageActor.h>
#include <vtkCubeAxesActor.h>
#include <vtkProperty.h>
#include <vtkTextProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWin32OpenGLRenderWindow.h>
#include <vtkWin32RenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <list>
#include <set>
#include <queue>
#include <algorithm>
#include <vector>

using namespace std;

#define square(a) (a) * (a)
#define VTK_NEW(type, name)\
	vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

void GetPointNeighborIds(vtkIdType pointId, vtkSmartPointer<vtkPolyData> mesh, set<vtkIdType>& neighbors)
{
	vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();
	//get all cells that vertex 'id' is a part of
	vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(pointId, cellIdList);

	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
		//cout << "id " << i << " : " << cellIdList->GetId(i) << endl;

		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

		//cout << "End points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;

		if (pointIdList->GetId(0) != pointId)
		{
			//cout << "Connected to " << pointIdList->GetId(0) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(0));
		}
		else
		{
			//cout << "Connected to " << pointIdList->GetId(1) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(1));
		}
	}
	for (vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
	{
		neighbors.insert(connectedVertices->GetId(i));
	}
}


void GetElasticityCoefficient(vector<double>& coefficients, vtkSmartPointer<vtkPolyData> data)
{
	// 获得网格数据中的所有边
	VTK_NEW(vtkExtractEdges, edges);
	edges->SetInputData(data);
	edges->Update();

	// 遍历所有的边，并计算每条边的弹性系数
	for (vtkIdType i = 0; i < edges->GetOutput()->GetNumberOfCells(); i++)
	{
		double A[3], B[3], C[3], D[3]; // A B C D四个点的坐标
		double l, l1, l2, l3, l4; // AB，AC，BC，AD，BD的长度
		double area1, area2; // 两个三角形的面积
		set<vtkIdType> neighborsC, neighborsD;
		vtkSmartPointer<vtkLine> line = vtkLine::SafeDownCast(edges->GetOutput()->GetCell(i));

		// 获得A、B、C、D四个点
		vtkIdType idA, idB, idC, idD;
		set<vtkIdType> ids;
		vector<vtkIdType> _ids;
		idA = line->GetPointIds()->GetId(0);
		idB = line->GetPointIds()->GetId(1);
		GetPointNeighborIds(idA, edges->GetOutput(), neighborsC);
		GetPointNeighborIds(idB, edges->GetOutput(), neighborsD);
		set_intersection(neighborsC.begin(), neighborsC.end(), neighborsD.begin(), neighborsD.end(), inserter(ids, ids.begin()));
		for (set<vtkIdType>::iterator it = ids.begin(); it != ids.end(); it++)
		{
			_ids.push_back(*it);
		}
		idC = _ids[0];
		idD = _ids[1];
		data->GetPoint(idA, A);
		data->GetPoint(idB, B);
		data->GetPoint(idC, C);
		data->GetPoint(idD, D);

		// 计算弹性系数
		double coefficient;
		l = vtkMath::Distance2BetweenPoints(A, B);
		l1 = vtkMath::Distance2BetweenPoints(A, C);
		l2 = vtkMath::Distance2BetweenPoints(B, C);
		l3 = vtkMath::Distance2BetweenPoints(A, D);
		l4 = vtkMath::Distance2BetweenPoints(B, D);
		area1 = vtkTriangle::TriangleArea(A, B, C);
		area2 = vtkTriangle::TriangleArea(A, B, D);
		coefficient = (square(l1) + square(l2) - square(l)) / area1
			+ (square(l3) + square(l4) - square(l)) / area2;
		coefficients.push_back(coefficient);
	}
}

int main(int argc, char *argv[])
{

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
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetThetaResolution(30);
		sphereSource->SetPhiResolution(15);
		sphereSource->Update();
		data = sphereSource->GetOutput();
	}
	vector<double> coefficents;
	GetElasticityCoefficient(coefficents, data);
	for (unsigned int i = 0; i < coefficents.size(); i++)
		cout << coefficents[i] << endl;
	system("pause");
	return EXIT_SUCCESS;
}