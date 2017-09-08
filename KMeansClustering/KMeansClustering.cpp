#include <vtkVersion.h>

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPoints.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkKMeansStatistics.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSTLReader.h>
#include <vtkSphereSource.h>

#include <sstream>

int main(int argc, char* argv[])
{
	// Create 2 clusters, one near (0,0,0) and the other near (3,3,3)

	vtkSmartPointer<vtkPolyData> data;
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	if (argc > 1)
	{
		reader->SetFileName(argv[1]);
		reader->Update();
		data = reader->GetOutput();
	}
	else
	{
		vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
		sphere->SetRadius(5.0);
		sphere->SetCenter(.0, .0, .0);
		sphere->Update();
		data = sphere->GetOutput();
	}
	vtkSmartPointer<vtkPoints> points = data->GetPoints();

	// Get the points into the format needed for KMeans
	vtkSmartPointer<vtkTable> inputData =
		vtkSmartPointer<vtkTable>::New();

	for (int c = 0; c < 3; ++c)
	{
		std::stringstream colName;
		colName << "coord " << c;
		vtkSmartPointer<vtkDoubleArray> doubleArray =
			vtkSmartPointer<vtkDoubleArray>::New();
		doubleArray->SetNumberOfComponents(1);
		doubleArray->SetName(colName.str().c_str());
		doubleArray->SetNumberOfTuples(points->GetNumberOfPoints());


		for (int r = 0; r < points->GetNumberOfPoints(); ++r)
		{
			double p[3];
			points->GetPoint(r, p);

			doubleArray->SetValue(r, p[c]);
		}

		inputData->AddColumn(doubleArray);
	}


	vtkSmartPointer<vtkKMeansStatistics> kMeansStatistics =
		vtkSmartPointer<vtkKMeansStatistics>::New();

#if VTK_MAJOR_VERSION <= 5
	kMeansStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, inputData);
#else
	kMeansStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, inputData);
#endif
	kMeansStatistics->SetColumnStatus(inputData->GetColumnName(0), 1);
	kMeansStatistics->SetColumnStatus(inputData->GetColumnName(1), 1);
	kMeansStatistics->SetColumnStatus(inputData->GetColumnName(2), 1);
	//kMeansStatistics->SetColumnStatus( "Testing", 1 );
	kMeansStatistics->RequestSelectedColumns();
	kMeansStatistics->SetAssessOption(true);
	kMeansStatistics->SetDefaultNumberOfClusters(4);
	kMeansStatistics->Update();

	// Display the results
	kMeansStatistics->GetOutput()->Dump();

	vtkSmartPointer<vtkIntArray> clusterArray =
		vtkSmartPointer<vtkIntArray>::New();
	clusterArray->SetNumberOfComponents(1);
	clusterArray->SetName("ClusterId");

	for (int r = 0; r < kMeansStatistics->GetOutput()->GetNumberOfRows(); r++)
	{
		vtkVariant v = kMeansStatistics->GetOutput()->GetValue(r, kMeansStatistics->GetOutput()->GetNumberOfColumns() - 1);
		std::cout << "Point " << r << " is in cluster " << v.ToInt() << std::endl;
		clusterArray->InsertNextValue(v.ToInt());
	}

	/*

	// Output the cluster centers

	vtkMultiBlockDataSet* outputMetaDS = vtkMultiBlockDataSet::SafeDownCast( kMeansStatistics->GetOutputDataObject( vtkStatisticsAlgorithm::OUTPUT_MODEL ) );
	vtkSmartPointer<vtkTable> outputMeta = vtkTable::SafeDownCast( outputMetaDS->GetBlock( 0 ) );
	//vtkSmartPointer<vtkTable> outputMeta = vtkTable::SafeDownCast( outputMetaDS->GetBlock( 1 ) );
	vtkDoubleArray* coord0 = vtkDoubleArray::SafeDownCast(outputMeta->GetColumnByName("coord 0"));
	vtkDoubleArray* coord1 = vtkDoubleArray::SafeDownCast(outputMeta->GetColumnByName("coord 1"));
	vtkDoubleArray* coord2 = vtkDoubleArray::SafeDownCast(outputMeta->GetColumnByName("coord 2"));

	for(unsigned int i = 0; i < coord0->GetNumberOfTuples(); ++i)
	{
	std::cout << coord0->GetValue(i) << " " << coord1->GetValue(i) << " " << coord2->GetValue(i) << std::endl;
	}
	*/

	vtkSmartPointer<vtkPolyData> polydata =
		vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	//polydata->GetPointData()->AddArray(clusterArray);
	polydata->GetPointData()->SetScalars(clusterArray);

	// Display
	vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
		vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	glyphFilter->SetInputConnection(polydata->GetProducerPort());
#else
	glyphFilter->SetInputData(polydata);
#endif
	glyphFilter->Update();

	// Create a mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(glyphFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(3);

	// Create a renderer, render window, and interactor
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actor to the scene
	renderer->AddActor(actor);
	//renderer->SetBackground(.3, .6, .3); // Background color green

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	renderWindowInteractor->SetInteractorStyle(style);

	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}