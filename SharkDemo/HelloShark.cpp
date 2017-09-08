#include <shark/Data/Csv.h>  
#include <shark/Algorithms/GradientDescent/CG.h>  
#include <shark/ObjectiveFunctions/ErrorFunction.h>  
#include <shark/ObjectiveFunctions/Loss/SquaredLoss.h>  
#include <shark/Models/LinearModel.h>  
#include <string>  

using namespace std;
using namespace shark;

RegressionDataset loadData(const string& dataFile, const string& labelFile)
{
	Data<RealVector> inputs;
	Data<RealVector> label;
	try
	{
		importCSV(inputs, dataFile, ' ');
		importCSV(label, labelFile, ' ');
	}
	catch (...)
	{
		cerr << "Unable to open file " << dataFile << " and/or " << labelFile << ". Check paths!" << endl;
		exit(EXIT_FAILURE);
	}
	RegressionDataset data(inputs, label);
	return data;
}

void main()
{
	//一定要更改这行代码中的两个文件路径，否则这里会出错！  
	RegressionDataset data = loadData("E:\\SourceCode\\Shark\\examples\\Supervised\\data\\regressionInputs.csv", "E:\\SourceCode\\Shark\\examples\\Supervised\\data\\regressionLabels.csv");

	RegressionDataset test = splitAtElement(data, static_cast<std::size_t>(0.8 * data.numberOfElements()));

	LinearModel<> model(inputDimension(data), labelDimension(data));

	SquaredLoss<> loss;
	ErrorFunction<RealVector, RealVector> errorFunction(data, &model, &loss);

	CG optimazer;
	optimazer.init(errorFunction);
	for (int i = 0; i < 100; i++)
	{
		optimazer.step(errorFunction);
	}

	double trainingError = optimazer.solution().value;

	model.setParameterVector(optimazer.solution().point);
	Data<RealVector> prediction = model(test.inputs());
	double testError = loss.eval(test.labels(), prediction);

	cout << "RESULTS: " << endl;
	cout << "======== \n" << endl;
	cout << "training error " << trainingError << endl;
	cout << "test error: " << testError << endl;

	getchar();
}