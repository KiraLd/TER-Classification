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

Mat convertFuzzy(Mat a)
{
	int x = a.rows;
	int y = a.cols;
	Mat b = Mat(x, y, CV_8U);
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			b.at<uchar>(i, j) = (uchar)((fabs(a.at<float>(i, j)) * 255));
		}
	}
	return b;
}

void FuzzyClustering::exportMembership(std::string file)
{
	std::string temp;
	std::ostringstream oss;
	Mat convert_;
	int j = 0;
	for (int i = 0; i < c;i++)
	{
		oss << i;
		temp = file + oss.str() + std::string(".jpg");
		convert_ = convertFuzzy(membership[i]);
		imwrite(temp, convert_);
		oss = std::ostringstream("");
		j++;
	}
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