#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkWindowToImageFilter.h>
#include <vtkContourFilter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkXMLPolyDataReader.h>
#include <vtkSTLReader.h>


int main(int argc, char *argv[])
{
	vtkSmartPointer<vtkPolyData> data3d;

	if (argc > 1)
	{
		vtkSmartPointer<vtkSTLReader> reader =
			vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(argv[1]);
		reader->Update();

		data3d = reader->GetOutput();
	}
	else
	{
		vtkSmartPointer<vtkSphereSource> source =
			vtkSmartPointer<vtkSphereSource>::New();
		source->SetCenter(0.0, 0.0, 5.0);
		source->SetRadius(2.0);
		source->SetPhiResolution(20.);
		source->SetThetaResolution(20.);
		source->Update();

		data3d = source->GetOutput();
	}

	double bounds_data[6], center_data[3];
	data3d->GetBounds(bounds_data);
	data3d->GetCenter(center_data);

	cout << "Xmin: " << bounds_data[0] << endl;
	cout << "Xmax: " << bounds_data[1] << endl;
	cout << "Ymin: " << bounds_data[2] << endl;
	cout << "Ymax: " << bounds_data[3] << endl;
	cout << "Zmin: " << bounds_data[4] << endl;
	cout << "Zmax: " << bounds_data[5] << endl;

	cout << "center (" << center_data[0] << ", " << center_data[1] << ", " << center_data[2] << ")" << endl;

	// Black and white scene with the data in order to print the view
	vtkSmartPointer<vtkPolyDataMapper> mapper_data =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper_data->SetInput(data3d);
#else
	mapper_data->SetInputData(data3d);
#endif

	vtkSmartPointer<vtkActor> actor_data =
		vtkSmartPointer<vtkActor>::New();
	actor_data->SetMapper(mapper_data);
	actor_data->GetProperty()->SetColor(0, 0, 0);

	vtkSmartPointer<vtkRenderer> tmp_rend =
		vtkSmartPointer<vtkRenderer>::New();
	tmp_rend->SetBackground(1, 1, 1);

	tmp_rend->AddActor(actor_data);
	tmp_rend->ResetCamera();
	tmp_rend->GetActiveCamera()->SetParallelProjection(1);

	vtkSmartPointer<vtkRenderWindow> tmp_rW =
		vtkSmartPointer<vtkRenderWindow>::New();
	tmp_rW->SetOffScreenRendering(1);
	tmp_rW->AddRenderer(tmp_rend);

	tmp_rW->Render();

	// Get a print of the window
	vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
		vtkSmartPointer<vtkWindowToImageFilter>::New();
	windowToImageFilter->SetInput(tmp_rW);
	windowToImageFilter->SetMagnification(2); //image quality
	windowToImageFilter->Update();

	// Extract the silhouette corresponding to the black limit of the image
	vtkSmartPointer<vtkContourFilter> ContFilter = vtkSmartPointer<vtkContourFilter>::New();
	ContFilter->SetInputConnection(windowToImageFilter->GetOutputPort());
	ContFilter->SetValue(0, 255);
	ContFilter->Update();

	// Make the contour coincide with the data.
	vtkSmartPointer<vtkPolyData> contour = ContFilter->GetOutput();

	double bounds_contour[6], center_contour[3];
	double trans_x = 0., trans_y = 0., trans_z = 0., ratio_x = 0., ratio_y = 0.;
	contour->GetBounds(bounds_contour);

	cout << "Xmin-contour: " << bounds_contour[0] << endl;
	cout << "Xmax-contour: " << bounds_contour[1] << endl;
	cout << "Ymin-contour: " << bounds_contour[2] << endl;
	cout << "Ymax-contour: " << bounds_contour[3] << endl;
	cout << "Zmin-contour: " << bounds_contour[4] << endl;
	cout << "Zmax-contour: " << bounds_contour[5] << endl;

	ratio_x = (bounds_data[1] - bounds_data[0]) / (bounds_contour[1] - bounds_contour[0]);
	ratio_y = (bounds_data[3] - bounds_data[2]) / (bounds_contour[3] - bounds_contour[2]);

	cout << "ratio-x: " << ratio_x << endl;
	cout << "ratio-y: " << ratio_y << endl;

	// Rescale the contour so that it shares the same bounds as the
	// input data
	vtkSmartPointer<vtkTransform>transform1 =
		vtkSmartPointer<vtkTransform>::New();
	transform1->Scale(ratio_x, ratio_y, 1.);

	vtkSmartPointer<vtkTransformPolyDataFilter> tfilter1 =
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	tfilter1->SetInput(contour);
#else
	tfilter1->SetInputData(contour);
#endif
	tfilter1->SetTransform(transform1);
	tfilter1->Update();

	contour = tfilter1->GetOutput();

	// Translate the contour so that it shares the same center as the
	// input data
	contour->GetCenter(center_contour);
	cout << "center (" << center_contour[0] << ", " << center_contour[1] << ", " << center_contour[2] << ")" << endl;
	trans_x = center_data[0] - center_contour[0];
	trans_y = center_data[1] - center_contour[1];
	trans_z = center_data[2] - center_contour[2];

	vtkSmartPointer<vtkTransform> transform2 =
		vtkSmartPointer<vtkTransform>::New();
	transform2->Translate(trans_x, trans_y, trans_z);

	vtkSmartPointer<vtkTransformPolyDataFilter> tfilter2 =
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	tfilter2->SetInput(contour);
#else
	tfilter2->SetInputData(contour);
#endif
	tfilter2->SetTransform(transform2);
	tfilter2->Update();

	contour = tfilter2->GetOutput();

	// Render the result : Input data + resulting silhouette

	// Updating the color of the data
	actor_data->GetProperty()->SetColor(0.9, 0.9, 0.8);

	// Create a mapper and actor of the silhouette
	vtkSmartPointer<vtkPolyDataMapper> mapper_contour =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper_contour->SetInput(contour);
#else
	mapper_contour->SetInputData(contour);
#endif

	vtkSmartPointer<vtkActor> actor_contour =
		vtkSmartPointer<vtkActor>::New();
	actor_contour->SetMapper(mapper_contour);
	actor_contour->GetProperty()->SetLineWidth(2.);

	// 2 renderers and a render window
	vtkSmartPointer<vtkRenderer> renderer1 =
		vtkSmartPointer<vtkRenderer>::New();
	renderer1->AddActor(actor_data);

	vtkSmartPointer<vtkRenderer> renderer2 =
		vtkSmartPointer<vtkRenderer>::New();
	renderer2->AddActor(actor_contour);

	vtkSmartPointer<vtkRenderWindow> renderwindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderwindow->SetSize(400, 400);

	renderwindow->AddRenderer(renderer1);
	renderer1->SetViewport(0., 0., 0.5, 1.);

	renderwindow->AddRenderer(renderer2);
	renderer2->SetViewport(0.5, 0., 1., 1.);

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	vtkSmartPointer<vtkRenderWindowInteractor> iren =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();

	iren->SetRenderWindow(renderwindow);
	iren->SetInteractorStyle(style);

	renderwindow->Render();
	iren->Start();

	return EXIT_SUCCESS;
}