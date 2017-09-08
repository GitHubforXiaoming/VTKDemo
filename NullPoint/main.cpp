#include "Utils.h"
using namespace std;

void DrawCircleOnMesh(vtkSmartPointer<vtkRenderer> renderer, vtkSmartPointer<vtkPolyData> source, int scheme)
{
	Utils utils;

	VTK_NEW(vtkPolyData, data);
	VTK_NEW(vtkDoubleArray, centers);

	data = utils.CellNormalsGenerator(source);													//����ÿ��cell������
	vtkSmartPointer<vtkDataArray> normals = data->GetCellData()->GetNormals();					//���ÿ��cell������
	utils.GetCenterPointsOfCells(source, centers);												//���ÿ��cell��������Ƭ��������

	VTK_NEW(vtkModifiedBSPTree, bspTree);														//���㶨λ��
	VTK_NEW(vtkCellLocator, cellLocator);														//cell��λ��
	set<vtkIdType> curCellIds;																	//��ǰ��ѡ�е�������Ƭ��id(�ļ���)

	bspTree->SetDataSet(source);
	bspTree->BuildLocator();
	cellLocator->SetDataSet(source);
	cellLocator->BuildLocator();

	/*�����ÿ��cell���ڵ�cell��ID*/
	set<vtkIdType> coverageIds;															//�洢Բ���ϵĵ���ӳ�䵽������Ƭ��id

	
	VTK_NEW(vtkColorSeries, colorSeries);
	int size = data->GetNumberOfCells();
	int numColors = colorSeries->GetNumberOfColors();
	vector<vector<double>> colorRange;
	colorSeries->SetColorScheme(scheme);
	
	for (int i = 0; i < numColors; i++)
	{
		vector<double> range(3);
		vtkColor3ub color = colorSeries->GetColor(i);
		range[0] = static_cast<double> (color[0]) / 255;
		range[1] = static_cast<double> (color[1]) / 255;
		range[2] = static_cast<double> (color[2]) / 255;

		colorRange.push_back(range);
	}
	for (vtkIdType i = 0; i < size; i++)
	{
		double curNormal[3];																//���������ʾ
		double curCenter[3];																//��ǰ������Ƭ�����ĵ�
		double circleCenter[3];																//Բ���ϵ�Բ��

		normals->GetTuple(i, curNormal);											//��õ�ǰ������Ƭ�ķ�����
		centers->GetTuple(i, curCenter);											//��õ�ǰ������Ƭ����������

		for (int layerNum = 0; layerNum < 3; layerNum++)
		{
			if (i % 100 == 0)
			{
				/*��������ֵ���������괴��Բ*/
				double radius = 0.5 * (layerNum + 1);
				VTK_NEW(vtkPoints, points);															//���ڻ���ֱ�ߵĵ㼯
				VTK_NEW(vtkPolyData, lines);
				VTK_NEW(vtkRegularPolygonSource, circle);									//��ǰ������Ƭ�Ϸ���Բ��

				utils.GenerateCircle(curNormal, curCenter, circle, radius);					//���ݷ�����������������������Ƭ�Ϸ�����Բ��
				circleCenter[0] = circle->GetCenter()[0];
				circleCenter[1] = circle->GetCenter()[1];
				circleCenter[2] = circle->GetCenter()[2];


				/*��ȡÿ����cell���ڵ�cell����������*/
				VTK_NEW(vtkDoubleArray, selectedCenters);									//ÿ����cell���ڵ�cell����������
				VTK_NEW(vtkDoubleArray, selectedNormals);									//ÿ����cell���ڵ�cell�ķ�����

				curCellIds.insert(i);


				/*�������*/
				vector<double> distances;													//������Ƭ�����ĵ㵽Բ��ľ���

				/*����Բ���Ͼ��ȵĵ㵽������Ƭ�ľ���*/
				VTK_NEW(vtkPoints, intersectPoints);
				VTK_NEW(vtkIdList, intersectCells);
				utils.CalcDistanceCircleToPoint(bspTree, cellLocator, intersectPoints, intersectCells,
					coverageIds, curNormal, circleCenter, points, distances, radius);
			

				VTK_NEW(vtkPolyDataMapper, circleMapper);
				VTK_NEW(vtkActor, circleActor);
				circleMapper->SetInputData(circle->GetOutput());
				circleActor->SetMapper(circleMapper);
				circleActor->GetProperty()->SetColor(colorRange[i % numColors][0], colorRange[i % numColors][1], colorRange[i % numColors][2]);
				renderer->AddActor(circleActor);

				/*����ֱ��*/
				utils.DrawLines(lines, points);
				lines->GetCellData()->SetScalars(utils.DyeForCells(lines->GetNumberOfLines()));

				VTK_NEW(vtkPolyDataMapper, linesMapper);
				VTK_NEW(vtkActor, linesActor);
				linesMapper->SetInputData(lines);

				linesActor->SetMapper(linesMapper);
				renderer->AddActor(linesActor);
			}
		}
	}
	
	/*�������ĸ��������뵽��Ⱦ����*/
	VTK_NEW(vtkPolyDataMapper, sourceMapper);
	

	VTK_NEW(vtkActor, sourceActor);

	sourceMapper->SetInputData(source);
	

	sourceMapper->SetInputData(source);
	sourceActor->SetMapper(sourceMapper);
	sourceActor->GetProperty()->SetOpacity(0.5);
	sourceActor->GetProperty()->EdgeVisibilityOn();
	sourceActor->GetProperty()->SetEdgeColor(171 / 255.0, 224 / 255.0, 255 / 255.0);
	sourceActor->GetProperty()->SetColor(45 / 255.0, 179 / 255.0, 160 / 255.0);

	renderer->AddActor(sourceActor);
	
}

int main(int argc, char* argv[])
{
	int scheme = 61;
	vtkSmartPointer<vtkPolyData> data;
	VTK_NEW(vtkSTLReader, reader);
	if (argc > 1)
	{
		reader->SetFileName(argv[1]);
		reader->Update();
		data = reader->GetOutput();
	}
	else
	{
		VTK_NEW(vtkSphereSource, sphere);
		sphere->SetRadius(5.0);
		sphere->SetCenter(.0, .0, .0);
		sphere->Update();
		data = sphere->GetOutput();
	}

	if (argc > 2)
	{
		scheme = atoi(argv[2]);
	}

	VTK_NEW(vtkRenderWindowInteractor, interactor);
	VTK_NEW(vtkRenderWindow, renderwindow);
	VTK_NEW(vtkInteractorStyleTrackballCamera, style);
	VTK_NEW(vtkRenderer, renderer);
	DrawCircleOnMesh(renderer, data, scheme);
	renderwindow->AddRenderer(renderer);
	interactor->SetRenderWindow(renderwindow);

	

	renderwindow->Render();
	interactor->SetInteractorStyle(style);
	interactor->Initialize();
	interactor->Start();

	return EXIT_SUCCESS;

}