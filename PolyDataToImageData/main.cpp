#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkSphereSource.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkSTLReader.h>

#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

/**
* This program generates a sphere (closed surface, vtkPolyData) and converts it into volume
* representation (vtkImageData) where the foreground voxels are 1 and the background voxels are
* 0. Internally vtkPolyDataToImageStencil is utilized. The resultant image is saved to disk
* in metaimage file format (SphereVolume.mhd).
*/
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
		sphereSource->SetRadius(20);
		sphereSource->SetPhiResolution(30);
		sphereSource->SetThetaResolution(30);
		sphereSource->Update();
		data = sphereSource->GetOutput();
	}

	vtkSmartPointer<vtkImageData> whiteImage =
		vtkSmartPointer<vtkImageData>::New();
	double bounds[6];
	data->GetBounds(bounds);
	double spacing[3]; // desired volume spacing
	spacing[0] = 0.5;
	spacing[1] = 0.5;
	spacing[2] = 0.5;
	whiteImage->SetSpacing(spacing);

	// compute dimensions
	int dim[3];
	for (int i = 0; i < 3; i++)
	{
		dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
	}
	whiteImage->SetDimensions(dim);
	whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

	double origin[3];
	origin[0] = bounds[0] + spacing[0] / 2;
	origin[1] = bounds[2] + spacing[1] / 2;
	origin[2] = bounds[4] + spacing[2] / 2;
	whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
	whiteImage->SetScalarTypeToUnsignedChar();
	whiteImage->AllocateScalars();
#else
	whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
#endif
	// fill the image with foreground voxels:
	unsigned char inval = 255;
	unsigned char outval = 0;
	vtkIdType count = whiteImage->GetNumberOfPoints();
	for (vtkIdType i = 0; i < count; ++i)
	{
		whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
	}

	// polygonal data --> image stencil:
	vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =
		vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
	pol2stenc->SetInput(pd);
#else
	pol2stenc->SetInputData(data);
#endif
	pol2stenc->SetOutputOrigin(origin);
	pol2stenc->SetOutputSpacing(spacing);
	pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
	pol2stenc->Update();

	// cut the corresponding white image and set the background:
	vtkSmartPointer<vtkImageStencil> imgstenc =
		vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
	imgstenc->SetInput(whiteImage);
	imgstenc->SetStencil(pol2stenc->GetOutput());
#else
	imgstenc->SetInputData(whiteImage);
	imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
	imgstenc->ReverseStencilOff();
	imgstenc->SetBackgroundValue(outval);
	imgstenc->Update();

	vtkSmartPointer<vtkMetaImageWriter> writer =
		vtkSmartPointer<vtkMetaImageWriter>::New();
	writer->SetFileName("SphereVolume.mhd");
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(imgstenc->GetOutput());
#else
	writer->SetInputData(imgstenc->GetOutput());
#endif
	writer->Write();

	vtkSmartPointer<vtkMetaImageReader> reader =
		vtkSmartPointer<vtkMetaImageReader>::New();
	reader->SetFileName("SphereVolume.mhd");
	reader->Update();

	// Visualize
	vtkSmartPointer<vtkImageActor> actor =
		vtkSmartPointer<vtkImageActor>::New();
	actor->GetMapper()->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->ResetCamera();

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<vtkInteractorStyleImage> style =
		vtkSmartPointer<vtkInteractorStyleImage>::New();

	renderWindowInteractor->SetInteractorStyle(style);

	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindowInteractor->Initialize();

	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}