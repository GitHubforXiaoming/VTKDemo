#include "Utils.h"
//#include "FragmentMosaicFrmView.h"

Utils::Utils()
{
	filePath = "";
}

string Utils::GetFileName()
{
	return filePath;
}

void Utils::SetSourceData(vtkSmartPointer<vtkPolyData> originData)
{
	this->originData = originData;
}

vtkSmartPointer<vtkPolyData> Utils::GetPolyDataFromPathSTL(string path)
{
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(path.c_str());
	reader->Update();
	filePath = path;
	vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
	return data;
}

void Utils::OpenAndDisplay(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkRenderer> renderer)
{
	VTK_NEW(vtkPolyData, source); 
	VTK_NEW(vtkPolyDataMapper, sourceMapper);
	VTK_NEW(vtkActor, sourceActor);

	sourceMapper->SetInputData(data);
	sourceActor->SetMapper(sourceMapper);
	sourceActor->GetProperty()->SetOpacity(0.1);
	sourceActor->GetProperty()->EdgeVisibilityOn();
	sourceActor->GetProperty()->SetEdgeColor(171 / 255.0, 224 / 255.0, 255 / 255.0);
	sourceActor->GetProperty()->SetColor(45 / 255.0, 179 / 255.0, 160 / 255.0);

	renderer->AddActor(sourceActor);
}

void Utils::GetCenterPointsOfCells(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkDoubleArray> &centersPoint)
{
	VTK_NEW(vtkCellCenters, centers);
	centers->SetInputData(data);
	centers->VertexCellsOn();
	centers->Update();

	centersPoint->SetNumberOfComponents(3);
	
	for (vtkIdType i = 0; i < centers->GetOutput()->GetNumberOfPoints(); i++)
	{
		double point[3];
		centers->GetOutput()->GetPoint(i, point);
		centersPoint->InsertNextTuple(point);
	}
}

void Utils::GetCellNormals(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkDataArray> &normals)
{
	normals = data->GetCellData()->GetNormals();
	if (normals)
	{
		return;
	}
}

void Utils::GenerateCircle(double normal[3], double originPoint[3], vtkSmartPointer<vtkRegularPolygonSource> &circle, double radius)
{
	circle->SetNormal(normal);
	circle->SetNumberOfSides(500);
	circle->SetRadius(radius);
	//求圆的圆心坐标
	double delta = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
	double x = normal[0] / delta + originPoint[0];
	double y = normal[1] / delta + originPoint[1];
	double z = normal[2] / delta + originPoint[2];

	circle->SetCenter(x, y, z);
	circle->Update();
}

vtkSmartPointer<vtkPolyData> Utils::CellNormalsGenerator(vtkSmartPointer<vtkPolyData> data)
{
	// Generate normals
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
#if VTK_MAJOR_VERSION <= 5
	normalGenerator->SetInput(polydata);
#else
	normalGenerator->SetInputData(data);
#endif
	normalGenerator->ComputePointNormalsOff();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->Update();
	
	// Optional settings
	/*normalGenerator->SetFeatureAngle(0.1);
	normalGenerator->SetSplitting(1);
	normalGenerator->SetConsistency(0);
	normalGenerator->SetAutoOrientNormals(0);
	normalGenerator->SetComputePointNormals(1);
	normalGenerator->SetComputeCellNormals(0);
	normalGenerator->SetFlipNormals(0);
	normalGenerator->SetNonManifoldTraversal(1);*/
	
	data = normalGenerator->GetOutput();
	return data;
}

void Utils::GetCellNeighborIds(vtkIdType cellId, vtkSmartPointer<vtkTriangleFilter> triangleFilter, set<vtkIdType>& neighbors)
{
	vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
	triangleFilter->GetOutput()->GetCellPoints(cellId, cellPointIds);

	for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkIdList> idList =
			vtkSmartPointer<vtkIdList>::New();

		//add one of the edge points
		idList->InsertNextId(cellPointIds->GetId(i));

		//add the other edge point
		if (i + 1 == cellPointIds->GetNumberOfIds())
		{
			idList->InsertNextId(cellPointIds->GetId(0));
		}
		else
		{
			idList->InsertNextId(cellPointIds->GetId(i + 1));
		}
		//get the neighbors of the cell
		vtkSmartPointer<vtkIdList> neighborCellIds =
			vtkSmartPointer<vtkIdList>::New();

		triangleFilter->GetOutput()->GetCellNeighbors(cellId, idList, neighborCellIds);

		for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++)
		{
			neighbors.insert(neighborCellIds->GetId(j));
		}
	}
}

void Utils::GetPointNeighborIds(vtkIdType pointId, vtkSmartPointer<vtkPolyData> mesh, set<vtkIdType>& neighbors)
{
	vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();
	//get all cells that vertex 'id' is a part of
	vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(pointId, cellIdList);

	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
		//cout << "id " << i << " : " << cellIdList->GetId(i) << endl;

		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

		//cout << "End points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;

		if (pointIdList->GetId(0) != pointId)
		{
			//cout << "Connected to " << pointIdList->GetId(0) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(0));
		}
		else
		{
			//cout << "Connected to " << pointIdList->GetId(1) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(1));
		}
	}
	for (vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
	{
		neighbors.insert(connectedVertices->GetId(i));
	}
}

void Utils::CreateDataSetFromIds(vtkSmartPointer<vtkExtractSelection> &extractSelection, vtkSmartPointer<vtkPolyData> source, set<vtkIdType> idSets)
{
	vtkSmartPointer<vtkIdTypeArray> ids =
		vtkSmartPointer<vtkIdTypeArray>::New();
	ids->SetNumberOfComponents(1);
	set<vtkIdType>::iterator it;
	for (it = idSets.begin(); it != idSets.end(); it++)
	{
		ids->InsertNextValue(*it);
	}
	vtkSmartPointer<vtkSelectionNode> selectionNode =
		vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(ids);

	vtkSmartPointer<vtkSelection> selection =
		vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);

	extractSelection->SetInputData(0, source);
	extractSelection->SetInputData(1, selection);
	extractSelection->Update();
}

void Utils::DrawLines(vtkSmartPointer<vtkPolyData> &linesData, vtkSmartPointer<vtkPoints> points)
{
	linesData->SetPoints(points);
	vtkSmartPointer<vtkCellArray> lines =
		vtkSmartPointer<vtkCellArray>::New();
	//绘制直线
	int size = points->GetNumberOfPoints();
	for (vtkIdType j = 0; j < size;)
	{
		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		line->GetPointIds()->SetId(0, j);
		line->GetPointIds()->SetId(1, j + 1);
		j += 2;
		lines->InsertNextCell(line);
	}
	linesData->SetLines(lines);
}

vtkSmartPointer<vtkUnsignedCharArray> Utils::DyeForCells(int numbersOfCells)
{
	// Create cell data
	vtkMath::RandomSeed(8775070); // for reproducibility
	vtkSmartPointer<vtkUnsignedCharArray> cellDataArray =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	cellDataArray->SetNumberOfComponents(3);
	cellDataArray->SetNumberOfTuples(numbersOfCells);
	for (int i = 0; i < numbersOfCells; i++)
	{
		float rgb[3];
		rgb[0] = vtkMath::Random(64, 255);
		rgb[1] = vtkMath::Random(64, 255);
		rgb[2] = vtkMath::Random(64, 255);
		cellDataArray->InsertTuple(i, rgb);
	}
	return cellDataArray;
}

void Utils::GetSelectedCellCenters(vtkSmartPointer<vtkDoubleArray> &selectedCenters, vtkSmartPointer<vtkDoubleArray> centers,
	vtkSmartPointer<vtkDoubleArray> &selectedNormals, vtkSmartPointer<vtkDataArray> normals, set<vtkIdType> ids)
{
	set<vtkIdType>::iterator it;
	selectedCenters->SetNumberOfComponents(3);
	for (it = ids.begin(); it != ids.end(); it++)
	{
		selectedCenters->InsertNextTuple(centers->GetTuple(*it));
		selectedNormals->InsertNextTuple(normals->GetTuple(*it));
	}
}


void Utils::GenerateBalancePoints(vtkSmartPointer<vtkDoubleArray> &balancePoints, double normal[3], double center[3], double radius)
{
	double alpha = normal[0];
	double beta = normal[1];
	double gamma = normal[2];
	double alphaPrime = sqrt(1 - square(alpha));
	double betaPrime = sqrt(1 - square(beta));
	double gammaPrime = sqrt(1 - square(gamma));
	balancePoints->SetNumberOfComponents(3);
	balancePoints->SetNumberOfTuples(POINT_NUM);
	for (int i = 0; i < POINT_NUM; i++)
	{
		double t = 2 * i * PI / POINT_NUM;
		double coordinate[3];
		coordinate[0] = center[0] + radius * (beta / gammaPrime) * cos(t) + radius * alpha * (gamma / gammaPrime) * sin(t);
		coordinate[1] = center[1] - radius * (alpha / gammaPrime) * cos(t) + radius * beta * (gamma / gammaPrime) * sin(t);
		coordinate[2] = center[2] - radius * gammaPrime * sin(t);
		balancePoints->SetTuple(i, coordinate);
	}
}

void Utils::CalcDistanceCircleToPoint(vtkSmartPointer<vtkModifiedBSPTree> bspTree, vtkSmartPointer<vtkCellLocator> cellLocator, vtkSmartPointer<vtkPoints> &intersectPoints, vtkSmartPointer<vtkIdList> &intersectCells,
	set<vtkIdType> &ids, double normal[3], double center[3], vtkSmartPointer<vtkPoints> &points, vector<double> &distances, double radius)
{
	VTK_NEW(vtkDoubleArray, pointsArray);
	
	double delta;
	double T; 
	delta = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
	T = 1 / delta + DELTA;
	
	GenerateBalancePoints(pointsArray, normal, center, radius);

	for (vtkIdType i = 0; i < pointsArray->GetNumberOfTuples(); i++)
	{
		// Inputs
		double pointOnCircle[3];
		double pointOverBottom[3];
		double tolerance = 1.e-8;
		pointsArray->GetTuple(i, pointOnCircle);

		//根据参数方程计算沿向量方向的直线上另一个点的坐标
		pointOverBottom[0] = -normal[0] * T + pointOnCircle[0];
		pointOverBottom[1] = -normal[1] * T + pointOnCircle[1];
		pointOverBottom[2] = -normal[2] * T + pointOnCircle[2];

		//找到圆周上的射线与提数据的交点
		double t; // Parametric coordinate of intersection (0 (corresponding to p1) to 1 (corresponding to p2))
		double coords[3]; // The coordinate of the intersection
		double pcoords[3];
		int subId;
		// Note: for a typical use case (ray-triangle intersection), pcoords and subId will not be used
		bspTree->IntersectWithLine(pointOverBottom, pointOnCircle, tolerance, intersectPoints, intersectCells);

		//交点可能不唯一，找到距离最小的角点
		if (intersectPoints->GetNumberOfPoints() != 0)
		{
			double minDistance = 0.0;
			int minIndex = 0;
			double alternativeDistance = vtkMath::Distance2BetweenPoints(pointOnCircle, intersectPoints->GetPoint(0));
			for (int j = 0; j < intersectPoints->GetNumberOfPoints(); j++)
			{
				minDistance = alternativeDistance;
				alternativeDistance = vtkMath::Distance2BetweenPoints(pointOnCircle, intersectPoints->GetPoint(j));
				if (minDistance >= alternativeDistance)
				{
					minDistance = vtkMath::Distance2BetweenPoints(pointOnCircle, intersectPoints->GetPoint(j));
					minIndex = j;
				}
			}
			intersectPoints->GetPoint(minIndex, coords);
		}
		else
		{
			bspTree->IntersectWithLine(pointOverBottom, pointOnCircle, tolerance, t, coords, pcoords, subId);
		}

		//找到交点所在的三角面片的id
		double closestPoint[3];//the coordinates of the closest point will be returned here
		double closestPointDist2; //the squared distance to the closest point will be returned here
		vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
		int subCellId; //this is rarely used (in triangle strips only, I believe)
		cellLocator->FindClosestPoint(coords, closestPoint, cellId, subCellId, closestPointDist2);

		ids.insert(cellId);

		//将圆周上的点和物体上对应的点一次插入到一个容器中，用于绘制直线
		points->InsertNextPoint(coords);
		points->InsertNextPoint(pointOnCircle);

		//根据两点坐标计算距离，并将其插入一个vector中
		double distance = vtkMath::Distance2BetweenPoints(pointOnCircle, coords);
		distances.push_back(distance);
	}
}

double Utils::GetVectorAngle(double *curVector, double *neiVector)
{
	double c = curVector[0] * neiVector[0] + curVector[1] * neiVector[1] + curVector[2] * neiVector[2];
	double d;
	d = sqrt(pow(curVector[0], 2) + pow(curVector[1], 2) + pow(curVector[2], 2)) * 
		sqrt(pow(neiVector[0], 2) + pow(neiVector[1], 2) + pow(neiVector[2], 2));
	double angle = acos(c / d);
	if (angle > 1.57)
		angle = 0;
	return angle;
}

short Utils::TransformVectorToPixel(vector<double> distances)
{
	short pixel = 0;
	int index = distances.size();
	for (vector<double>::iterator distance = distances.begin(); distance != distances.end(); distance++)
	{
		if (*distance > 1)
		{
			pixel += pow(2, index - 1);
		}
		index--;
	}
	return pixel;
}


void Utils::CurvatureThreshold(vtkSmartPointer<vtkPolyData> data, vector<double>& curvatureArray)
{
	VTK_NEW(vtkCurvatures, curvaturesFilter);

	curvaturesFilter->SetInputData(data);
	curvaturesFilter->SetCurvatureTypeToMean();
	curvaturesFilter->Update();

	//获得每个网格数据的曲率
	vtkSmartPointer<vtkDataArray> curvatures =
		static_cast<vtkSmartPointer<vtkDataArray>>(curvaturesFilter->GetOutput()->GetPointData()->GetAttribute(0));
	for (vtkIdType i = 0; i < curvatures->GetNumberOfTuples(); i++)
	{
		curvatureArray.push_back(curvatures->GetTuple(i)[0]);
	}
}


void Utils::ReplacePost(string&s1, const string&s2, const string&s3)
{
	string::size_type pos = 0;
	string::size_type a = s2.size();
	string::size_type b = s3.size();
	while ((pos = s1.find(s2, pos)) != string::npos)
	{
		s1.replace(pos, a, s3);
		pos += b;
	}
}

void Utils::GetViewport(unsigned int size, vector<double*>&  viewports)
{
	int rows = ceil(sqrt(size));
	int *cols = new int[rows];
	int capacityOfRow = floor(sqrt(size));
	int sum = 0;
	for (int i = 0; i < rows - 1; i++)
	{
		cols[i] = capacityOfRow;
		sum += capacityOfRow;
	}
	cols[rows - 1] = 0;
	int residue = size - sum;
	int curIndex = rows - 1;
	for (int j = 1 ; j <= residue; j++)
	{
		if (cols[curIndex] != capacityOfRow + 1)
			cols[curIndex]++;
		else
		{
			curIndex--;
			j--;
		}
	}
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < cols[row]; col++)
		{
			double viewport[4] = {
				col * 1.0 / cols[row],
				row * 1.0 / rows,
				(col + 1) * 1.0 / cols[row],
				(row + 1) * 1.0 / rows
			};
			viewports.push_back(viewport);
		}
	}
}