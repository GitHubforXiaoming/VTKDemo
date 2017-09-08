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
#include <vtkViewTheme.h>
// For glyphing
#include <vtkReverseSense.h>
#include <vtkMaskPoints.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkTransform.h>
#include <vtkGraphToGlyphs.h>
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
#include <vtkMutableUndirectedGraph.h>
#include <vtkCircularLayoutStrategy.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkGraphLayoutView.h>
#include <vtkIntArray.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderedGraphRepresentation.h>

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

void GetGraphFromMesh(vtkSmartPointer<vtkMutableUndirectedGraph>& graph, vtkSmartPointer<vtkPolyData> data)
{
	vector<vtkIdType> cellIds;
	int size = data->GetNumberOfPoints();
	vtkIdType* vertexs = new vtkIdType[size];
	int* visited = new int[size];
	for (int i = 0; i < size; i++)
	{
		visited[i] = 0;
	}
	for (int i = 0; i < size; i++)
	{
		vertexs[i] = graph->AddVertex();
	}

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(data);
	triangleFilter->Update();

	vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInputConnection(triangleFilter->GetOutputPort());
	extractEdges->Update();

	vtkSmartPointer<vtkPolyData> mesh = extractEdges->GetOutput();

	for (vtkIdType i = 0; i < data->GetNumberOfPoints(); i++)
	{
		set<vtkIdType> neighbors;
		GetPointNeighborIds(i, mesh, neighbors);
		for (set<vtkIdType>::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); neighbor++)
		{
			if (visited[*neighbor] == 0)
			{
				graph->AddEdge(vertexs[i], vertexs[*neighbor]);
				visited[i] = 1;
			}
		}
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



	vtkSmartPointer<vtkMutableUndirectedGraph> g =
		vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	GetGraphFromMesh(g, data);

	// Create the color array
	vtkSmartPointer<vtkIntArray> edgeColors =
		vtkSmartPointer<vtkIntArray>::New();
	edgeColors->SetNumberOfComponents(1);
	edgeColors->SetName("Color");

	vtkSmartPointer<vtkLookupTable> lookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	lookupTable->SetNumberOfTableValues(data->GetNumberOfPoints());

	for (vtkIdType i = 0; i < data->GetNumberOfPoints(); i++)
	{
		double rgb[3];
		rgb[0] = vtkMath::Random(0, 255);
		rgb[1] = vtkMath::Random(0, 255);
		rgb[2] = vtkMath::Random(0, 255);
		edgeColors->InsertNextValue(i);
		lookupTable->SetTableValue(i, rgb);
	}
	lookupTable->Build();
	// Add the color array to the graph
	g->GetEdgeData()->AddArray(edgeColors);
	vtkSmartPointer<vtkGraphLayoutView> graphLayoutView =
		vtkSmartPointer<vtkGraphLayoutView>::New();
	graphLayoutView->AddRepresentationFromInput(g);
	graphLayoutView->SetEdgeColorArrayName("Color");
	graphLayoutView->ColorEdgesOn();

	vtkSmartPointer<vtkViewTheme> theme =
		vtkSmartPointer<vtkViewTheme>::New();
	theme->SetPointLookupTable(lookupTable);

	graphLayoutView->ApplyViewTheme(theme);
	vtkRenderedGraphRepresentation::SafeDownCast(graphLayoutView->GetRepresentation())->SetGlyphType(vtkGraphToGlyphs::CIRCLE);
	graphLayoutView->ResetCamera();
	graphLayoutView->Render();
	graphLayoutView->GetInteractor()->Start();
	return EXIT_SUCCESS;
}