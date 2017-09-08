#include <vtkVersion.h>
#include "vtkImageData.h"
#include "vtkImageShiftScale.h"
#include "vtkFastSplatter.h"
#include "vtkImageViewer2.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkSTLReader.h"
#include "vtkSphereSource.h"

#include <cmath>

const int SPLAT_IMAGE_SIZE = 100;

int main(int argc, char *argv[])
{
	// For the purposes of this example we'll build the splat image by
	// hand.

	vtkSmartPointer<vtkImageData> splatImage =
		vtkSmartPointer<vtkImageData>::New();
	splatImage->SetDimensions(SPLAT_IMAGE_SIZE, SPLAT_IMAGE_SIZE, 1);
#if VTK_MAJOR_VERSION <= 5
	splatImage->SetScalarTypeToFloat();
	splatImage->SetNumberOfScalarComponents(1);
	splatImage->AllocateScalars();
#else
	splatImage->AllocateScalars(VTK_FLOAT, 1);
#endif
	for (int i = 0; i < SPLAT_IMAGE_SIZE; ++i)
	{
		for (int j = 0; j < SPLAT_IMAGE_SIZE; ++j)
		{
			double xCoord = 1 - fabs((i - SPLAT_IMAGE_SIZE / 2)
				/ (SPLAT_IMAGE_SIZE / 2.0));
			double yCoord = 1 - fabs((j - SPLAT_IMAGE_SIZE / 2)
				/ (SPLAT_IMAGE_SIZE / 2.0));

			splatImage->SetScalarComponentFromDouble(i, j, 0, 0,
				xCoord * yCoord);
		}
	}

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


	vtkSmartPointer<vtkFastSplatter> splatter =
		vtkSmartPointer<vtkFastSplatter>::New();
#if VTK_MAJOR_VERSION <= 5
	splatter->SetInput(splatPoints);
	splatter->SetInput(1, splatImage);
#else
	splatter->SetInputData(data);
	splatter->SetInputData(1, splatImage);
#endif
	splatter->SetOutputDimensions(2 * SPLAT_IMAGE_SIZE,
		2 * SPLAT_IMAGE_SIZE,
		1);

	// The image viewers and writers are only happy with unsigned char
	// images.  This will convert the floats into that format.
	vtkSmartPointer<vtkImageShiftScale> resultScale =
		vtkSmartPointer<vtkImageShiftScale>::New();
	resultScale->SetOutputScalarTypeToUnsignedChar();
	resultScale->SetShift(0);
	resultScale->SetScale(255);
	resultScale->SetInputConnection(splatter->GetOutputPort());

	splatter->Update();
	resultScale->Update();

	// Set up a viewer for the image.  vtkImageViewer and
	// vtkImageViewer2 are convenient wrappers around vtkActor2D,
	// vtkImageMapper, vtkRenderer, and vtkRenderWindow.  All you need
	// to supply is the interactor and hooray, Bob's your uncle.
	vtkSmartPointer<vtkImageViewer2> imageViewer =
		vtkSmartPointer<vtkImageViewer2>::New();
	imageViewer->SetInputConnection(resultScale->GetOutputPort());
	imageViewer->SetColorLevel(127);
	imageViewer->SetColorWindow(255);

	vtkSmartPointer<vtkRenderWindowInteractor> iren =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	imageViewer->SetupInteractor(iren);

	imageViewer->Render();
	imageViewer->GetRenderer()->ResetCamera();

	imageViewer->Render();
	iren->Start();

	return EXIT_SUCCESS;
}