#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkOctreePointLocator.h>

int main(int, char *[])
{
	// Create some random points
	vtkSmartPointer<vtkPointSource> pointSource =
		vtkSmartPointer<vtkPointSource>::New();
	pointSource->SetNumberOfPoints(10);
	pointSource->Update();

	// Create the tree
	vtkSmartPointer<vtkOctreePointLocator> octree =
		vtkSmartPointer<vtkOctreePointLocator>::New();
	octree->SetDataSet(pointSource->GetOutput());
	octree->BuildLocator();

	// Find the k closest points to (0,0,0)
	vtkIdType k = 2;
	double testPoint[3] = { 0.0, 0.0, 0.0 };
	vtkSmartPointer<vtkIdList> result =
		vtkSmartPointer<vtkIdList>::New();

	octree->FindClosestNPoints(k, testPoint, result);

	for (vtkIdType i = 0; i < k; i++)
	{
		vtkIdType point_ind = result->GetId(i);
		double p[3];
		pointSource->GetOutput()->GetPoint(point_ind, p);
		std::cout << "Closest point " << i << ": Point "
			<< point_ind << ": (" << p[0] << ", "
			<< p[1] << ", " << p[2] << ")" << std::endl;
	}
	system("pause");
	return EXIT_SUCCESS;
}