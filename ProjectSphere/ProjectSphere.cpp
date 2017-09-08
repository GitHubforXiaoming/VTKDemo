#include <vtkMinimalStandardRandomSequence.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
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

int main(int vtkNotUsed(argc), char *vtkNotUsed(argv)[])
{
	vtkSmartPointer<vtkMinimalStandardRandomSequence> randomSequence
		= vtkSmartPointer<vtkMinimalStandardRandomSequence>::New();
	randomSequence->SetSeed(1);

	vtkSmartPointer<vtkParametricFunctionSource> parametricFunctionSource
		= vtkSmartPointer<vtkParametricFunctionSource>::New();
	parametricFunctionSource->SetUResolution(64);
	parametricFunctionSource->SetVResolution(64);
	parametricFunctionSource->SetWResolution(64);
	parametricFunctionSource->SetScalarModeToNone();
	parametricFunctionSource->GenerateTextureCoordinatesOff();

	parametricFunctionSource->SetOutputPointsPrecision(vtkAlgorithm::SINGLE_PRECISION);

	vtkSmartPointer<vtkParametricEllipsoid> parametricEllipsoid
		= vtkSmartPointer<vtkParametricEllipsoid>::New();

	randomSequence->Next();
	double xRadius = randomSequence->GetValue();
	parametricEllipsoid->SetXRadius(xRadius);

	randomSequence->Next();
	double yRadius = randomSequence->GetValue();
	parametricEllipsoid->SetYRadius(yRadius);

	randomSequence->Next();
	double zRadius = randomSequence->GetValue();
	parametricEllipsoid->SetZRadius(zRadius);

	parametricFunctionSource->SetParametricFunction(parametricEllipsoid);

	parametricFunctionSource->Update();

	vtkSmartPointer<vtkPolyData> polyData = parametricFunctionSource->GetOutput();
	vtkSmartPointer<vtkPoints> points = polyData->GetPoints();

	if (points->GetDataType() != VTK_FLOAT)
	{
		return EXIT_FAILURE;
	}

	parametricFunctionSource->SetOutputPointsPrecision(vtkAlgorithm::DOUBLE_PRECISION);

	randomSequence->Next();
	xRadius = randomSequence->GetValue();
	parametricEllipsoid->SetXRadius(xRadius);

	randomSequence->Next();
	yRadius = randomSequence->GetValue();
	parametricEllipsoid->SetYRadius(yRadius);

	randomSequence->Next();
	zRadius = randomSequence->GetValue();
	parametricEllipsoid->SetZRadius(zRadius);

	parametricFunctionSource->SetParametricFunction(parametricEllipsoid);

	parametricFunctionSource->Update();

	polyData = parametricFunctionSource->GetOutput();
	points = polyData->GetPoints();

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputData(polyData);
	actor->SetMapper(mapper);

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

	if (points->GetDataType() != VTK_DOUBLE)
	{
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

