#include "Utils.h"
using namespace std;

void DrawCircleOnMesh(vtkSmartPointer<vtkRenderer> renderer, vtkSmartPointer<vtkPolyData> source, int scheme)
{
	Utils utils;

	VTK_NEW(vtkPolyData, data);
	VTK_NEW(vtkDoubleArray, centers);

	data = utils.CellNormalsGenerator(source);													//生成每个cell的向量
	vtkSmartPointer<vtkDataArray> normals = data->GetCellData()->GetNormals();					//获得每个cell的向量
	utils.GetCenterPointsOfCells(source, centers);												//获得每个cell（三角面片）的中心

	VTK_NEW(vtkModifiedBSPTree, bspTree);														//交点定位树
	VTK_NEW(vtkCellLocator, cellLocator);														//cell定位器
	set<vtkIdType> curCellIds;																	//当前所选中的三角面片的id(的集合)

	bspTree->SetDataSet(source);
	bspTree->BuildLocator();
	cellLocator->SetDataSet(source);
	cellLocator->BuildLocator();

	/*获得与每个cell相邻的cell的ID*/
	set<vtkIdType> coverageIds;															//存储圆周上的点所映射到三角面片的id

	
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
		double curNormal[3];																//向量坐标表示
		double curCenter[3];																//当前三角面片的中心点
		double circleCenter[3];																//圆面上的圆心

		normals->GetTuple(i, curNormal);											//获得当前三角面片的法向量
		centers->GetTuple(i, curCenter);											//获得当前三角面片的中心坐标

		for (int layerNum = 0; layerNum < 3; layerNum++)
		{
			if (i % 100 == 0)
			{
				/*根据向量值和中心坐标创建圆*/
				double radius = 0.5 * (layerNum + 1);
				VTK_NEW(vtkPoints, points);															//用于绘制直线的点集
				VTK_NEW(vtkPolyData, lines);
				VTK_NEW(vtkRegularPolygonSource, circle);									//当前三角面片上方的圆面

				utils.GenerateCircle(curNormal, curCenter, circle, radius);					//根据法向量和中心坐标在三角面片上方绘制圆面
				circleCenter[0] = circle->GetCenter()[0];
				circleCenter[1] = circle->GetCenter()[1];
				circleCenter[2] = circle->GetCenter()[2];


				/*获取每个与cell相邻的cell的中心坐标*/
				VTK_NEW(vtkDoubleArray, selectedCenters);									//每个与cell相邻的cell的中心坐标
				VTK_NEW(vtkDoubleArray, selectedNormals);									//每个与cell相邻的cell的法向量

				curCellIds.insert(i);


				/*计算距离*/
				vector<double> distances;													//三角面片的中心点到圆面的距离

				/*计算圆周上均匀的点到三角面片的距离*/
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

				/*绘制直线*/
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
	
	/*将创建的各种网格入到渲染器中*/
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