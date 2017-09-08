#include <iostream>
#include <vector>

using namespace std;

int main()
{
	unsigned int size = 7;
	vector<unsigned int> block;
	vector<unsigned int> acc;
	block.push_back(4);
	block.push_back(3);

	unsigned int sum = 0;
	for (unsigned int i = 0; i < block.size(); i++)
	{
		sum += block[i];
		acc.push_back(sum);
	}

	unsigned int** map = new unsigned int*[15];
	for (int i = 0; i < size; i++)
		map[i] = new unsigned int[2];


	for (unsigned int i = 0; i < block.size(); i++)
	{
		for (unsigned int row = 0; row < size; row++)
		{
			for (unsigned int col = 0; col < block[i]; col++)
			{
				if (acc[i] - row < block[i] + 1)
				{
					map[row][0] = row;
					map[row][1] = acc[i];
				}
			}
		}
	}
	
	/*for (unsigned int row = 0; row < 15; row++)
	{
		cout << map[row][0] << "-" << map[row][1];
		cout << endl;
	}*/
	 
	for (int i = 0; i < size - 1; i++)
	{
		for (int j = i + 1 ; j < size; j++)
		{
			if (j > map[i][1] - 1)
			{
				cout << i << "--" << j << endl;
			}
		}
	}
	system("pause");
	return 0;
}