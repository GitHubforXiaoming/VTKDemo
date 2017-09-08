#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/legacy/legacy.hpp>

using namespace cv;
using namespace std;


void GetImageFromMatrixGray(Mat& image, int** pixels, int size)
{
	for (int i = 0; i < size; i++)
	{
		uchar* data = image.ptr<uchar>(i);
		for (int j = 0; j < size; j++)
		{
			data[j] = pixels[i][j];
		}
	}
}

int main()
{
	const int size = 5002;
	int** matrix = new int*[size];
	for (int i = 0; i < size; i++)
		matrix[i] = new int[size];
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j || i == 2 * j)
				matrix[i][j] = 0;
			else
				matrix[i][j] = 255;
		}
	}

	Mat out = Mat(size, size, CV_8UC1);
	GetImageFromMatrixGray(out, matrix, size); 
	for (int i = 0; i++; i < size)
	{
		delete matrix[i];
	}
	delete matrix;
	//cout << out;
	imwrite("1.jpg", out);

	

	return 0;
}