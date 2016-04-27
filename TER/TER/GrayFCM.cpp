#include "GrayFCM.h"

void GrayFCM::init()
{
	RNG rng;
	if (membership != NULL)
	{
		delete[] membership;
	}
	membership = new Mat[c];
	for (int i = 0; i < c;i++)
	{
		membership[i] = Mat(x, y, CV_32F);
		for (int j = 0;j < x;j++)
		{
			for (int k = 0; k < y;k++)
			{
				membership[i].at<float>(j, k) = rng.uniform(0.f, 1.f);
			}
		}
	}
}

void GrayFCM::update_centers()
{
	float sum;
	Mat* exp = new Mat[c];
	for (int i = 0; i < c;i++)
	{
		centers[i] = 0.f;
		sum = 0.f;
		pow(membership[i], m, exp[i]);
		for (int j = 0; j < x; j++)
		{
			for (int k = 0; k < y;k++)
			{
				centers[i] += exp[i].at<float>(j, k)*img.at<uchar>(j, k);
				sum += exp[i].at<float>(j, k);
			}
		}
		centers[i] /= sum;
		std::cout << centers[i] << std::endl;
	}
	delete[] exp;
}

void GrayFCM::dissimilarity()
{
	for (int i = 0; i < c; i++)
	{
		D[i] = Mat(x, y, CV_32F);
		for (int j = 0; j < x; j++)
		{
			for (int k = 0; k < y; k++)
			{
				D[i].at<float>(j, k) = abs(img.at<uchar>(j, k) - centers[i]);
			}
		}
	}
}

void GrayFCM::update_membership()
{
	float sum;
	float exp;
	if (m < 1.5)
	{
		m = 1.5;
	}
	exp = (1.f / (m - 1));

	for (int i = 0; i < c; i++)
	{
		for (int j = 0; j < x; j++)
		{
			for (int k = 0; k < y; k++)
			{
				sum = 0.f;
				membership[i].at<float>(j, k) = std::powf(1. / D[i].at<float>(j, k), exp);
				for (int l = 0; l < c;l++)
				{
					sum += std::powf(1. / D[l].at<float>(j, k), exp);
				}
				membership[i].at<float>(j, k) /= sum;
			}
		}
	}
}

GrayFCM::GrayFCM(): FuzzyClustering()
{
	D = NULL;
}

GrayFCM::~GrayFCM()
{
	if (D != NULL)
	{
		delete[] D;
	}
}

void GrayFCM::exec(int c, float m, float e, int i_max)
{

	cvtColor(img, img, CV_BGR2GRAY);
	x = img.rows;
	y = img.cols;
	this->c = c;
	this->m = m;
	this->e = e;
	if (centers != NULL)
	{
		delete[] centers;
	}
	centers = new float[c];
	if (D != NULL)
	{
		delete[] D;
	}
	D = new Mat[c];
	init();
	bool arret = false;
	int i = 0;
	float mean = 0.f;
	for (int j = 0; j < c;j++)
	{
		mean += centers[j];
	}
	mean /= (c);
	float mean_1;
	while (!arret && i < i_max)
	{
		update_centers();
		dissimilarity();
		update_membership();
		i++;
		std::cout << i << std::endl;
		mean_1 = mean;
		for (int j = 0; j < c;j++)
		{
			mean += centers[j];
		}
		mean /= (c);
		if (abs(mean - mean_1) < e)
		{
			arret = true;
		}
	}
}