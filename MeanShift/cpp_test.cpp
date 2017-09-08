#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include "MeanShift.h"
#include <string>
#include <sstream>

using namespace std;

std::vector<std::string> splitWithStl(const std::string &str, const std::string &pattern)
{
	std::vector<std::string> resVec;

	if ("" == str)
	{
		return resVec;
	}
	//方便截取最后一段数据
	std::string strs = str + pattern;

	size_t pos = strs.find(pattern);
	size_t size = strs.size();

	while (pos != std::string::npos)
	{
		std::string x = strs.substr(0, pos);
		resVec.push_back(x);
		strs = strs.substr(pos + 1, size);
		pos = strs.find(pattern);
	}

	return resVec;
}


string readFileIntoString(char * filename)
{
	ifstream ifile(filename);
	//将文件读入到ostringstream对象buf中
	ostringstream buf;
	char ch;
	while (buf&&ifile.get(ch))
		buf.put(ch);
	//返回与流对象buf关联的字符串
	return buf.str();
}

vector<vector<double> > load_points(const char *filename) {
	vector<vector<double> > points;
	ifstream fin(filename);
	string line;
	while (fin >> line)
	{
		vector<string> splitted = splitWithStl(line, ",");
		vector<double> point;
		for (unsigned int i = 4; i < splitted.size(); i++)
		{
			point.push_back(atof(splitted[i].c_str()));
		}
		points.push_back(point);
	}
	fin.close();
	return points;
}

void print_points(vector<vector<double> > points){
    for(int i=0; i<points.size(); i++){
        for(int dim = 0; dim<points[i].size(); dim++) {
            printf("%f ", points[i][dim]);
        }
        printf("\n");
    }
}

double epan_kernel(double distance, double kernel_bandwidth){
	if (distance < kernel_bandwidth)
	{
		return 1 - (distance * distance) / (kernel_bandwidth * kernel_bandwidth);
	}
	return 0;
}

int main(int argc, char **argv)
{
	MeanShift *msp = new MeanShift(epan_kernel);
    double kernel_bandwidth = 0.5;

    vector<vector<double> > points = load_points("1-point-sets.txt");
    vector<Cluster> clusters = msp->cluster(points, kernel_bandwidth);

    FILE *fp = fopen("result.csv", "w");
    if(!fp){
        perror("Couldn't write result.csv");
        exit(0);
    }

    printf("\n====================\n");
    printf("Found %lu clusters\n", clusters.size());
    printf("====================\n\n");
    for(int cluster = 0; cluster < clusters.size(); cluster++) {
      printf("Cluster %i:\n", cluster);
      for(int point = 0; point < clusters[cluster].original_points.size(); point++){
        for(int dim = 0; dim < clusters[cluster].original_points[point].size(); dim++) {
          printf("%f ", clusters[cluster].original_points[point][dim]);
          fprintf(fp, dim?",%f":"%f", clusters[cluster].original_points[point][dim]);
        }
        printf(" -> ");
        for(int dim = 0; dim < clusters[cluster].shifted_points[point].size(); dim++) {
          printf("%f ", clusters[cluster].shifted_points[point][dim]);
        }
        printf("\n");
        fprintf(fp, "\n");
      }
      printf("\n");
    }
    fclose(fp);
	system("pause");
    return 0;
}
