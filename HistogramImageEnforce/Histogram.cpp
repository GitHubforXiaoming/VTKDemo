#include <opencv2/opencv.hpp>

using namespace cv;

int ImageStretchByHistogram(IplImage *src1, IplImage *dst1)
/*************************************************
Function:      ͨ��ֱ��ͼ�任����ͼ����ǿ����ͼ��Ҷȵ���ֵ���쵽0-255
src1:               ��ͨ���Ҷ�ͼ��
dst1:              ͬ����С�ĵ�ͨ���Ҷ�ͼ��
*************************************************/
{
	assert(src1->width == dst1->width);
	double p[256], p1[256], num[256];

	memset(p, 0, sizeof(p));
	memset(p1, 0, sizeof(p1));
	memset(num, 0, sizeof(num));
	int height = src1->height;
	int width = src1->width;
	long wMulh = height * width;

	//statistics  
	for (int x = 0; x<src1->width; x++)
	{
		for (int y = 0; y<src1->height; y++){
			uchar v = ((uchar*)(src1->imageData + src1->widthStep*y))[x];
			num[v]++;
		}
	}
	//calculate probability  
	for (int i = 0; i<256; i++)
	{
		p[i] = num[i] / wMulh;
	}

	//p1[i]=sum(p[j]);  j<=i;  
	for (int i = 0; i<256; i++)
	{
		for (int k = 0; k <= i; k++)
			p1[i] += p[k];
	}

	// histogram transformation  
	for (int x = 0; x<src1->width; x++)
	{
		for (int y = 0; y<src1->height; y++){
			uchar v = ((uchar*)(src1->imageData + src1->widthStep*y))[x];
			((uchar*)(dst1->imageData + dst1->widthStep*y))[x] = p1[v] * 255 + 0.5;
		}
	}
	return 0;
}


int main()
{
	IplImage* src = cvLoadImage("2.bmp", 0);
	cvShowImage("��ǿǰ", src);
	IplImage* dst = cvCreateImage(cvGetSize(src), src->depth, src->nChannels);
	ImageStretchByHistogram(src, dst);
	cvShowImage("��ǿ��", dst);
	waitKey(0);
}