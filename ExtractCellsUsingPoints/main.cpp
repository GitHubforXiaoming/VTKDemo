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
#include <vtkInformation.h>
#include <vtkDataSetSurfaceFilter.h>

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

using namespace std;

void CurvatureThreshold(vtkSmartPointer<vtkPolyData> data, vector<double>& curvatureArray)
{
	vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();

	curvaturesFilter->SetInputData(data);
	curvaturesFilter->SetCurvatureTypeToMean();
	curvaturesFilter->Update();

	//获得每个网格数据的曲率
	vtkSmartPointer<vtkDataArray> curvatures =
		static_cast<vtkSmartPointer<vtkDataArray>>(curvaturesFilter->GetOutput()->GetPointData()->GetAttribute(0));
	for (vtkIdType i = 0; i < curvatures->GetNumberOfTuples(); i++)
	{
		curvatureArray.push_back(curvatures->GetTuple(i)[0]);
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



	std::cout << "There are " << data->GetNumberOfPoints()
		<< " input points." << std::endl;
	std::cout << "There are " << data->GetNumberOfCells()
		<< " input cells." << std::endl;

	vtkSmartPointer<vtkIdTypeArray> ids =
		vtkSmartPointer<vtkIdTypeArray>::New();
	ids->SetNumberOfComponents(1);

	// Set values
	double maxCurvature, minCurvature;
	double threshold;
	vector<double> curvatureArrays;
	CurvatureThreshold(data, curvatureArrays);

	maxCurvature = *max_element(curvatureArrays.begin(), curvatureArrays.end());
	minCurvature = *min_element(curvatureArrays.begin(), curvatureArrays.end());
	threshold = minCurvature + (maxCurvature - minCurvature) * 4 / 10.0;
	for (vtkIdType i = 0; i < data->GetNumberOfPoints(); i++)
	{
		if (curvatureArrays[i] > threshold)
		{
			ids->InsertNextValue(i);
		}
	}

	vtkSmartPointer<vtkSelectionNode> selectionNode =
		vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::POINT);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(ids);
	selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(), 1);

	vtkSmartPointer<vtkSelection> selection =
		vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);

	vtkSmartPointer<vtkExtractSelection> extractSelection =
		vtkSmartPointer<vtkExtractSelection>::New();

	extractSelection->SetInputData(0, data);
#if VTK_MAJOR_VERSION <= 5
	extractSelection->SetInput(1, selection);
#else
	extractSelection->SetInputData(1, selection);
#endif
	extractSelection->Update();

	// In selection
	vtkSmartPointer<vtkUnstructuredGrid> selected =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected->ShallowCopy(extractSelection->GetOutput());

	std::cout << "There are " << selected->GetNumberOfPoints()
		<< " points in the selection." << std::endl;
	std::cout << "There are " << selected->GetNumberOfCells()
		<< " cells in the selection." << std::endl;

	// Get points that are NOT in the selection
	selectionNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1); //invert the selection
	extractSelection->Update();

	vtkSmartPointer<vtkUnstructuredGrid> notSelected =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	notSelected->ShallowCopy(extractSelection->GetOutput());

	std::cout << "There are " << notSelected->GetNumberOfPoints()
		<< " points NOT in the selection." << std::endl;
	std::cout << "There are " << notSelected->GetNumberOfCells()
		<< " cells NOT in the selection." << std::endl;

	vtkSmartPointer<vtkDataSetMapper> inputMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	inputMapper->SetInputData(data);
	vtkSmartPointer<vtkActor> inputActor =
		vtkSmartPointer<vtkActor>::New();
	inputActor->SetMapper(inputMapper);

	vtkSmartPointer<vtkDataSetMapper> selectedMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilterSelected =
		vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	surfaceFilterSelected->SetInputData(selected);
	surfaceFilterSelected->Update();
#if VTK_MAJOR_VERSION <= 5
	selectedMapper->SetInputConnection(selected->GetProducerPort());
#else
	selectedMapper->SetInputData(surfaceFilterSelected->GetOutput());
#endif
	vtkSmartPointer<vtkActor> selectedActor =
		vtkSmartPointer<vtkActor>::New();
	selectedActor->SetMapper(selectedMapper);

	vtkSmartPointer<vtkDataSetMapper> notSelectedMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilterNotSelected =
		vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	surfaceFilterNotSelected->SetInputData(notSelected);
	surfaceFilterNotSelected->Update();
	vtkSmartPointer<vtkMutableUndirectedGraph> g =
		vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	GetGraphFromMesh(g, surfaceFilterSelected->GetOutput());

	// Create the color array
	vtkSmartPointer<vtkIntArray> edgeColors =
		vtkSmartPointer<vtkIntArray>::New();
	edgeColors->SetNumberOfComponents(1);
	edgeColors->SetName("Color");

	vtkSmartPointer<vtkLookupTable> lookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	lookupTable->SetNumberOfTableValues(surfaceFilterSelected->GetOutput()->GetNumberOfPoints());

	for (vtkIdType i = 0; i < surfaceFilterSelected->GetOutput()->GetNumberOfPoints(); i++)
	{
		double rgb[3];
		rgb[0] = vtkMath::Random(64, 255) / 255.0;
		rgb[1] = vtkMath::Random(64, 255) / 255.0;
		rgb[2] = vtkMath::Random(64, 255) / 255.0;
		cout << "(" << rgb[0] << "," << rgb[1] << "," << rgb[2] << ")" << endl;
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