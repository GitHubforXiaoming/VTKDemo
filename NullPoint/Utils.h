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
	vtkSmartPointer<vtkPolyData> originData;							//ԭʼ������	

public:

	Utils();

	~Utils(){};

	void SetSourceData(vtkSmartPointer<vtkPolyData> originData);

	/**************************************************************************************************
	// ������: GetPolyDataFromPathSTL
	// ����: ��ȡ(.stl)�ļ�
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: vtkPolyData
	// ��������: string --- �ļ���·��
	**************************************************************************************************/
	vtkSmartPointer<vtkPolyData> GetPolyDataFromPathSTL(string path);

	void OpenAndDisplay(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkRenderer> renderer);

	string GetFileName();

	/**************************************************************************************************
	// ������: GetCellNeighborIds
	// ����: ����ĳ��������Ƭ��ID����ȡ��������Ƭ���ھ�ID
	// ����: xiaoming
	// ʱ��: 2017/04/05
	// ����ֵ: void
	// ��������: vtkIdType cellId --- ��ǰ������Ƭ��id
				vtkSmartPointer<vtkTriangleFilter> triangleFilter --- ������Ĺ�����
				set<vtkIdType>& neighbors --- ���ص��ھ�������Ƭ��id����
	**************************************************************************************************/
	void GetCellNeighborIds(vtkIdType cellId, vtkSmartPointer<vtkTriangleFilter> triangleFilter, set<vtkIdType>& neighbors);



	/**************************************************************************************************
	// ������: GetPointNeighborIds
	// ����: ��ȡ��ÿ��cell�е�ÿ�������ڵ�cell��id
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: vtkIdType, vtkPolyData, vtkIdType
	//         --- �ھӵ�ids���������ݣ�ָ����cell��id
	**************************************************************************************************/
	void GetPointNeighborIds(vtkIdType pointId, vtkSmartPointer<vtkPolyData> mesh, set<vtkIdType>& connectedVertices);






	/**************************************************************************************************
	// ������: GetCellNormals
	// ����: ��ȡÿ��cell������
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: vtkPolyData, vtkDoubleArray
	//		   --- ��������ݣ��洢����������
	**************************************************************************************************/
	void GetCellNormals(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkDataArray> &normals);





	/**************************************************************************************************
	// ������: GetCenterPointsOfCells
	// ����: ��ȡÿ��cell�����ĵ�
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: vtkPolyData, vtkCellCenters
	//		   --- ��������ݣ����ĵ������
	**************************************************************************************************/
	void GetCenterPointsOfCells(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkDoubleArray> &centersPoint);





	/**************************************************************************************************
	// ������: CreateDataSetFromIds
	// ����: ����id�б�����Ӧdataset
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: vtkPolyData, vtkExtractSelection, vtkIdTypeArray
	//		   --- ԭʼ���ݣ�����ids��ȡ���������ݣ������id�б�
	**************************************************************************************************/
	void CreateDataSetFromIds(vtkSmartPointer<vtkExtractSelection> &extractSelection, vtkSmartPointer<vtkPolyData> source, set<vtkIdType> idSets);





	/**************************************************************************************************
	// ������: CreateDataSetFromIdList
	// ����: ����id�б�����Ӧdataset
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: vtkPolyData, vtkExtractSelection, vtkIdTypeArray
	//		   --- ԭʼ���ݣ�����ids��ȡ���������ݣ������id�б�
	**************************************************************************************************/
	void CreateDataSetFromIdList(vtkSmartPointer<vtkPolyData> originData, vtkSmartPointer<vtkExtractSelection> &extractSelection, vtkSmartPointer<vtkIdList> idList);






	/**************************************************************************************************
	// ������: CalcDistancePointToPlane
	// ����: ����㵽Բ��ľ���
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: vtkPlane, vtkPolyData, vtkDoubleArray�� vtkPoints
	//		   --- ָ����ƽ�棬����ɸѡ������ݼ����洢�����double���͵����飬�洢ƽ�����ݼ���Դ���ݼ��϶�Ӧ�ĵ㣨һ��һ�洢��
	**************************************************************************************************/
	void CalcDistancePointToCircle(vtkSmartPointer<vtkRegularPolygonSource> circle, double center[3], 
		vtkSmartPointer<vtkDoubleArray> centerArray, vtkSmartPointer<vtkDoubleArray> normalArray,
		vector<double> &distances, vtkSmartPointer<vtkPoints> &points);






	/**************************************************************************************************
	// ������: GenerateCircle
	// ����: ����ƽ��
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: double [], double[], vtkPlaneSource
	//		   --- ƽ��ķ���������ƽ���ԭ����أ����ɵ�ƽ�����ݼ�
	**************************************************************************************************/
	void GenerateCircle(double normal[3], double originPoint[3], vtkSmartPointer<vtkRegularPolygonSource> &circle, double radius);






	/**************************************************************************************************
	// ������: GenerateSphere
	// ����: ����ƽ��
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: double [], double[], vtkPlaneSource
	//		   --- ƽ��ķ���������ƽ���ԭ����أ����ɵ�ƽ�����ݼ�
	**************************************************************************************************/
	void GenerateSphere(double normal[3], double originPoint[3], vtkSmartPointer<vtkSphereSource> &sphere);






	/**************************************************************************************************
	// ������: GenerateCirclePoints
	// ����: ����ƽ��
	// ����: xiaoming
	// ʱ��: 2016/11/03
	// ����ֵ: void
	// ��������: double [], double[], vtkPlaneSource
	//		   --- ƽ��ķ���������ƽ���ԭ����أ����ɵ�ƽ�����ݼ�
	**************************************************************************************************/
	void GenerateCirclePoints(double normal[3], double originPoint[3], vtkSmartPointer<vtkPointSource> &circlePoints, double radius);






	/**************************************************************************************************
	// ������: DrawLines
	// ����: �����ṩ�ĵ㼯������
	// ����: xiaoming
	// ʱ��: 2016/10/09
	// ����ֵ: void
	// ��������: vtkPolyData, vtkPoints
	//		   --- 	�洢���������ݼ����������ṩ�ĵ������
	**************************************************************************************************/
	void DrawLines(vtkSmartPointer<vtkPolyData> &lineData, vtkSmartPointer<vtkPoints> points);







	/**************************************************************************************************
	// ������: CellNormalsGenerator
	// ����: ���������ݵķ�����
	// ����: xiaoming
	// ʱ��: 2016/10/15
	// ����ֵ: vtkSmartPointer<vtkPolyData>
	**************************************************************************************************/
	vtkSmartPointer<vtkPolyData> CellNormalsGenerator(vtkSmartPointer<vtkPolyData> data);






	/**************************************************************************************************
	// ������: DyeForCells
	// ����: Ϊÿ�������ݵ�cellȾɫ
	// ����: xiaoming
	// ʱ��: 2016/10/15
	// ����ֵ: vtkSmartPointer<vtkUnsignedCharArray>
	**************************************************************************************************/
	vtkSmartPointer<vtkUnsignedCharArray> DyeForCells(int numbersOfCells);





	/**************************************************************************************************
	// ������: RandomColors
	// ����: Ϊÿ�������ݵ�cell�������Ⱦɫ
	// ����: xiaoming
	// ʱ��: 2017/04/05
	// ����ֵ: void
	**************************************************************************************************/
	void RandomColors(int con, double* rgb);




	/**************************************************************************************************
	// ������: GetSelectedCellCenters
	// ����: �����ѡ�е�������Ƭ
	// ����: xiaoming
	// ʱ��: 2016/10/15
	// ����ֵ: void
	**************************************************************************************************/
	void GetSelectedCellCenters(vtkSmartPointer<vtkDoubleArray> &selectedCenters, vtkSmartPointer<vtkDoubleArray> centers,
		vtkSmartPointer<vtkDoubleArray> &selectedNormals, vtkSmartPointer<vtkDataArray> normals, set<vtkIdType> ids);




	/**************************************************************************************************
	// ������: DrawPoints
	// ����: ���Ƶ㼯
	// ����: xiaoming
	// ʱ��: 2016/10/15
	// ����ֵ: void
	**************************************************************************************************/
	void DrawPoints(vtkSmartPointer<vtkPolyData> &pointData, vtkSmartPointer<vtkDoubleArray> randomPoints);





	/**************************************************************************************************
	// ������: GenerateRandomPoints
	// ����: ��Բ��������λ�þ��ȵ�30���㣨�������⣩
	// ����: xiaoming
	// ʱ��: 2016/11/06
	// ����ֵ: void
	//		   file:///D:/projects/FragmentMosaic/DataProcess/��ά����ϵ��Բ�ܵĲ�������.gif
	**************************************************************************************************/
	void GenerateBalancePoints(vtkSmartPointer<vtkDoubleArray> &balancePoints, double normal[3], double center[3], double radius);






	/**************************************************************************************************
	// ������: CalcDistanceCircleToPoint
	// ����: ����Բ���Ͼ��ȵĵ㵽��������������Ƭ�ľ���
	// ����: xiaoming
	// ʱ��: 2016/10/15
	// ����ֵ: void
	**************************************************************************************************/
	void CalcDistanceCircleToPoint(vtkSmartPointer<vtkModifiedBSPTree> bspTree, vtkSmartPointer<vtkCellLocator> cellLocator, vtkSmartPointer<vtkPoints> &intersectPoints, vtkSmartPointer<vtkIdList> &intersectCells, 
		set<vtkIdType> &ids, double normal[3], double center[3], vtkSmartPointer<vtkPoints> &points, vector<double> &distances, double radius);






	/**************************************************************************************************
	// ������: CalcDistanceCircleToPoint
	// ����: ����Բ���Ͼ��ȵĵ㵽��������������Ƭ�ľ���
	// ����: xiaoming
	// ʱ��: 2016/10/15
	// ����ֵ: void
	**************************************************************************************************/
	void CalcDistanceSphereToPoint(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkSphereSource> sphere,
		set<vtkIdType> &ids, double normal[3], vector<double> &distances, vtkSmartPointer<vtkPoints> &points);





	/**************************************************************************************************
	// ������: GetAllFormatFiles
	// ����: �����ļ���·������ȡ��ǰ�ļ����µģ�.stl���ļ�
	// ����: xiaoming
	// ʱ��: 2017/04/05
	// ����ֵ: void
	// ��������: string path --- �ļ���·��
				vector<string>& files --- �ļ����µ����и���Ҫ���Ҫ����ļ���
				string format --- ���ڹ��ˣ�.stl���ļ�
	**************************************************************************************************/
	void GetAllFormatFiles(string path, vector<string>& files, string format);




	/**************************************************************************************************
	// ������: ReplacePost
	// ����: �����ļ��ĺ�׺
	// ����: xiaoming
	// ʱ��: 2017/04/05
	// ����ֵ: void
	// ��������: string& s1, const string& s2, const string& s3 --- ԭ�ַ�����ԭ�ַ������滻�Ĳ��֣��滻���ַ���
	**************************************************************************************************/
	void ReplacePost(string& s1, const string& s2, const string& s3);




	void CurvatureThreshold(vtkSmartPointer<vtkPolyData> data, vector<double>& curvatureArray);

	double GetVectorAngle(double *curVector, double *neiVector);




	/**************************************************************************************************
	// ������: TransformVectorToPixel
	// ����: ��һ������ת����һ������ֵ
	// ����: xiaoming
	// ʱ��: 2017/04/05
	// ����ֵ: void
	// ��������: vector<double> distances --- ��ת��������
	**************************************************************************************************/
	short TransformVectorToPixel(vector<double> distances);


	void GetViewport(unsigned int size, vector<double*>& viewports);




};