// The source
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
#include <vtkMutableUndirectedGraph.h>
#include <vtkCircularLayoutStrategy.h>
#include <vtkDataSetAttributes.h>
#include <vtkGraphLayoutView.h>
#include <vtkMutableUndirectedGraph.h>
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
// For glyphing
#include <vtkReverseSense.h>
#include <vtkMaskPoints.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkTransform.h>
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


#include <list>
#include <set>
#include <queue>
#include <algorithm>
#include <math.h>

#define square(a) (a) * (a)
#define VTK_NEW(type, name)\
	vtkSmartPointer<type> name = vtkSmartPointer<type>::New()
#define VTK_DEF(type, name)\
	vtkSmartPointer<type> name
#define VTK_INSTANCE(type, name)\
	name = vtkSmartPointer<type>::New()

#define LAYER_NUM 3
#define RADIUS_DELTA 0.25
#define RADIUS_EXTERNAL 0.5
#define RADIUS_INTERNAL 0.25
#define DELTA 10000
#define POINT_NUM 8
#define ANGLE_ZONE 8
#define PI vtkMath::Pi()

typedef unsigned short ushort;
enum FILE_OPTION { DISTANCE = 1, COEFFICENT = 2 };

#pragma once

using namespace std;

class Utils
{

private:
	string filePath;
	vtkSmartPointer<vtkPolyData> originData;							//原始体数据	

public:

	Utils();

	~Utils(){};

	void SetSourceData(vtkSmartPointer<vtkPolyData> originData);

	/**************************************************************************************************
	// 函数名: GetPolyDataFromPathSTL
	// 描述: 获取(.stl)文件
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: vtkPolyData
	// 参数类型: string --- 文件的路径
	**************************************************************************************************/
	vtkSmartPointer<vtkPolyData> GetPolyDataFromPathSTL(string path);

	void OpenAndDisplay(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkRenderer> renderer);

	string GetFileName();

	/**************************************************************************************************
	// 函数名: GetCellNeighborIds
	// 描述: 根据某个三角面片的ID，获取该三角面片的邻居ID
	// 作者: xiaoming
	// 时间: 2017/04/05
	// 返回值: void
	// 参数类型: vtkIdType cellId --- 当前三角面片的id
				vtkSmartPointer<vtkTriangleFilter> triangleFilter --- 三角面的过滤器
				set<vtkIdType>& neighbors --- 返回的邻居三角面片的id集合
	**************************************************************************************************/
	void GetCellNeighborIds(vtkIdType cellId, vtkSmartPointer<vtkTriangleFilter> triangleFilter, set<vtkIdType>& neighbors);



	/**************************************************************************************************
	// 函数名: GetPointNeighborIds
	// 描述: 获取与每个cell中的每个点相邻的cell的id
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: vtkIdType, vtkPolyData, vtkIdType
	//         --- 邻居的ids，输入数据，指定的cell的id
	**************************************************************************************************/
	void GetPointNeighborIds(vtkIdType pointId, vtkSmartPointer<vtkPolyData> mesh, set<vtkIdType>& connectedVertices);






	/**************************************************************************************************
	// 函数名: GetCellNormals
	// 描述: 获取每个cell的向量
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: vtkPolyData, vtkDoubleArray
	//		   --- 输入的数据，存储向量的数组
	**************************************************************************************************/
	void GetCellNormals(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkDataArray> &normals);





	/**************************************************************************************************
	// 函数名: GetCenterPointsOfCells
	// 描述: 获取每个cell的中心点
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: vtkPolyData, vtkCellCenters
	//		   --- 输入的数据，中心点的容器
	**************************************************************************************************/
	void GetCenterPointsOfCells(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkDoubleArray> &centersPoint);





	/**************************************************************************************************
	// 函数名: CreateDataSetFromIds
	// 描述: 根据id列表创建相应dataset
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: vtkPolyData, vtkExtractSelection, vtkIdTypeArray
	//		   --- 原始数据，根据ids提取出来的数据，输入的id列表
	**************************************************************************************************/
	void CreateDataSetFromIds(vtkSmartPointer<vtkExtractSelection> &extractSelection, vtkSmartPointer<vtkPolyData> source, set<vtkIdType> idSets);





	/**************************************************************************************************
	// 函数名: CreateDataSetFromIdList
	// 描述: 根据id列表创建相应dataset
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: vtkPolyData, vtkExtractSelection, vtkIdTypeArray
	//		   --- 原始数据，根据ids提取出来的数据，输入的id列表
	**************************************************************************************************/
	void CreateDataSetFromIdList(vtkSmartPointer<vtkPolyData> originData, vtkSmartPointer<vtkExtractSelection> &extractSelection, vtkSmartPointer<vtkIdList> idList);






	/**************************************************************************************************
	// 函数名: CalcDistancePointToPlane
	// 描述: 计算点到圆面的距离
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: vtkPlane, vtkPolyData, vtkDoubleArray， vtkPoints
	//		   --- 指定的平面，经过筛选后的数据集，存储距离的double类型的数组，存储平面数据集与源数据集上对应的点（一对一存储）
	**************************************************************************************************/
	void CalcDistancePointToCircle(vtkSmartPointer<vtkRegularPolygonSource> circle, double center[3], 
		vtkSmartPointer<vtkDoubleArray> centerArray, vtkSmartPointer<vtkDoubleArray> normalArray,
		vector<double> &distances, vtkSmartPointer<vtkPoints> &points);






	/**************************************************************************************************
	// 函数名: GenerateCircle
	// 描述: 绘制平面
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: double [], double[], vtkPlaneSource
	//		   --- 平面的法向量，与平面的原点相关，生成的平面数据集
	**************************************************************************************************/
	void GenerateCircle(double normal[3], double originPoint[3], vtkSmartPointer<vtkRegularPolygonSource> &circle, double radius);






	/**************************************************************************************************
	// 函数名: GenerateSphere
	// 描述: 绘制平面
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: double [], double[], vtkPlaneSource
	//		   --- 平面的法向量，与平面的原点相关，生成的平面数据集
	**************************************************************************************************/
	void GenerateSphere(double normal[3], double originPoint[3], vtkSmartPointer<vtkSphereSource> &sphere);






	/**************************************************************************************************
	// 函数名: GenerateCirclePoints
	// 描述: 绘制平面
	// 作者: xiaoming
	// 时间: 2016/11/03
	// 返回值: void
	// 参数类型: double [], double[], vtkPlaneSource
	//		   --- 平面的法向量，与平面的原点相关，生成的平面数据集
	**************************************************************************************************/
	void GenerateCirclePoints(double normal[3], double originPoint[3], vtkSmartPointer<vtkPointSource> &circlePoints, double radius);






	/**************************************************************************************************
	// 函数名: DrawLines
	// 描述: 根据提供的点集绘制线
	// 作者: xiaoming
	// 时间: 2016/10/09
	// 返回值: void
	// 参数类型: vtkPolyData, vtkPoints
	//		   --- 	存储绘制线数据集，根据所提供的点绘制线
	**************************************************************************************************/
	void DrawLines(vtkSmartPointer<vtkPolyData> &lineData, vtkSmartPointer<vtkPoints> points);







	/**************************************************************************************************
	// 函数名: CellNormalsGenerator
	// 描述: 生成体数据的法向量
	// 作者: xiaoming
	// 时间: 2016/10/15
	// 返回值: vtkSmartPointer<vtkPolyData>
	**************************************************************************************************/
	vtkSmartPointer<vtkPolyData> CellNormalsGenerator(vtkSmartPointer<vtkPolyData> data);






	/**************************************************************************************************
	// 函数名: DyeForCells
	// 描述: 为每个体数据的cell染色
	// 作者: xiaoming
	// 时间: 2016/10/15
	// 返回值: vtkSmartPointer<vtkUnsignedCharArray>
	**************************************************************************************************/
	vtkSmartPointer<vtkUnsignedCharArray> DyeForCells(int numbersOfCells);





	/**************************************************************************************************
	// 函数名: RandomColors
	// 描述: 为每个体数据的cell“随机”染色
	// 作者: xiaoming
	// 时间: 2017/04/05
	// 返回值: void
	**************************************************************************************************/
	void RandomColors(int con, double* rgb);




	/**************************************************************************************************
	// 函数名: GetSelectedCellCenters
	// 描述: 获得所选中的三角面片
	// 作者: xiaoming
	// 时间: 2016/10/15
	// 返回值: void
	**************************************************************************************************/
	void GetSelectedCellCenters(vtkSmartPointer<vtkDoubleArray> &selectedCenters, vtkSmartPointer<vtkDoubleArray> centers,
		vtkSmartPointer<vtkDoubleArray> &selectedNormals, vtkSmartPointer<vtkDataArray> normals, set<vtkIdType> ids);




	/**************************************************************************************************
	// 函数名: DrawPoints
	// 描述: 绘制点集
	// 作者: xiaoming
	// 时间: 2016/10/15
	// 返回值: void
	**************************************************************************************************/
	void DrawPoints(vtkSmartPointer<vtkPolyData> &pointData, vtkSmartPointer<vtkDoubleArray> randomPoints);





	/**************************************************************************************************
	// 函数名: GenerateRandomPoints
	// 描述: 在圆周上生成位置均匀的30个点（存在问题）
	// 作者: xiaoming
	// 时间: 2016/11/06
	// 返回值: void
	//		   file:///D:/projects/FragmentMosaic/DataProcess/三维坐标系中圆周的参数方程.gif
	**************************************************************************************************/
	void GenerateBalancePoints(vtkSmartPointer<vtkDoubleArray> &balancePoints, double normal[3], double center[3], double radius);






	/**************************************************************************************************
	// 函数名: CalcDistanceCircleToPoint
	// 描述: 计算圆周上均匀的点到体数据中三角面片的距离
	// 作者: xiaoming
	// 时间: 2016/10/15
	// 返回值: void
	**************************************************************************************************/
	void CalcDistanceCircleToPoint(vtkSmartPointer<vtkModifiedBSPTree> bspTree, vtkSmartPointer<vtkCellLocator> cellLocator, vtkSmartPointer<vtkPoints> &intersectPoints, vtkSmartPointer<vtkIdList> &intersectCells, 
		set<vtkIdType> &ids, double normal[3], double center[3], vtkSmartPointer<vtkPoints> &points, vector<double> &distances, double radius);






	/**************************************************************************************************
	// 函数名: CalcDistanceCircleToPoint
	// 描述: 计算圆周上均匀的点到体数据中三角面片的距离
	// 作者: xiaoming
	// 时间: 2016/10/15
	// 返回值: void
	**************************************************************************************************/
	void CalcDistanceSphereToPoint(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkSphereSource> sphere,
		set<vtkIdType> &ids, double normal[3], vector<double> &distances, vtkSmartPointer<vtkPoints> &points);





	/**************************************************************************************************
	// 函数名: GetAllFormatFiles
	// 描述: 根据文件夹路径，获取当前文件夹下的（.stl）文件
	// 作者: xiaoming
	// 时间: 2017/04/05
	// 返回值: void
	// 参数类型: string path --- 文件夹路径
				vector<string>& files --- 文件夹下的所有复合要求的要求的文件名
				string format --- 用于过滤（.stl）文件
	**************************************************************************************************/
	void GetAllFormatFiles(string path, vector<string>& files, string format);




	/**************************************************************************************************
	// 函数名: ReplacePost
	// 描述: 更改文件的后缀
	// 作者: xiaoming
	// 时间: 2017/04/05
	// 返回值: void
	// 参数类型: string& s1, const string& s2, const string& s3 --- 原字符串，原字符串待替换的部分，替换的字符串
	**************************************************************************************************/
	void ReplacePost(string& s1, const string& s2, const string& s3);




	void CurvatureThreshold(vtkSmartPointer<vtkPolyData> data, vector<double>& curvatureArray);

	double GetVectorAngle(double *curVector, double *neiVector);




	/**************************************************************************************************
	// 函数名: TransformVectorToPixel
	// 描述: 将一个向量转换成一个像素值
	// 作者: xiaoming
	// 时间: 2017/04/05
	// 返回值: void
	// 参数类型: vector<double> distances --- 待转换的向量
	**************************************************************************************************/
	short TransformVectorToPixel(vector<double> distances);


	void GetViewport(unsigned int size, vector<double*>& viewports);




};