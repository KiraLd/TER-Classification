#include "FuzzyClustering.hpp"

FuzzyClustering::FuzzyClustering()
{
	x = 0;
	y = 0;
	membership = NULL;
	centers = NULL;
}

FuzzyClustering::~FuzzyClustering()
{
	if (membership == NULL)
	{
		delete[] membership;
	}
	if (centers == NULL)
	{
		delete[] centers;
	}
}

void FuzzyClustering::setImage(const Mat& img)
{
	this->img = img.clone();
}

void FuzzyClustering::setImage(std::string file, int type)
{
	img = imread(file, type);
	if (!img.data)
	{
		std::cout << "Fichier introuvable" << std::endl;
	}
}

void FuzzyClustering::exportMembership(std::string file)
{

}

void FuzzyClustering::show(int time)
{
	if (time < 1000)
	{
		time = 1000;
	}
	if (membership != NULL)
	{
		for (int i = 0; i < c; i++)
		{
			namedWindow("test");
			imshow("test", membership[i]);
			waitKey(time);
		}
	}
	
}