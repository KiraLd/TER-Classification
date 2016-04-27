#include "RgbFCM.h"

void RgbFCM::init()
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

void RgbFCM::update_centers()
{
	float sum;
	Mat* exp = new Mat[c];
	for (int i = 0; i < c;i++)
	{
		centers[3*i] = 0.f;
		centers[3 * i+1] = 0.f;
		centers[3 * i+2] = 0.f;
		sum = 0.f;
		pow(membership[i], m, exp[i]);
		for (int j = 0; j < x; j++)
		{
			for (int k = 0; k < y;k++)
			{
				centers[3*i] += exp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[0];
				centers[3*i+1] += exp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[1];
				centers[3*i+2] += exp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[2];
				sum += exp[i].at<float>(j, k);
			}
		}
		centers[3*i] /= sum;
		centers[3 * i+1] /= sum;
		centers[3 * i+2] /= sum;
	//	std::cout << centers[i] << std::endl;
	}
	delete[] exp;
}

void RgbFCM::dissimilarity()
{
	float r, g, b;
	for (int i = 0; i < c; i++)
	{
		D[i] = Mat(x, y, CV_32F);
		for (int j = 0; j < x; j++)
		{
			r = 0.f;
			g = 0.f;
			b = 0.f;
			for (int k = 0; k < y; k++)
			{
				r = img.at<Vec3b>(j, k)[0] - centers[3 * i];
				g = img.at<Vec3b>(j, k)[1] - centers[3 * i+1];
				b = img.at<Vec3b>(j, k)[2] - centers[3 * i+2];
				r *= r;
				g *= g;
				b *= b;
				D[i].at<float>(j, k) = sqrtf(r+g+b);
			}
		}
	}
}

void RgbFCM::update_membership()
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

RgbFCM::RgbFCM() : FuzzyClustering()
{
	D = NULL;
}

RgbFCM::~RgbFCM()
{
	if (D != NULL)
	{
		delete[] D;
	}
}

void RgbFCM::exec(int c, float m, float e, int i_max)
{

	x = img.rows;
	y = img.cols;
	this->c = c;
	this->m = m;
	this->e = e;
	if (centers != NULL)
	{
		delete[] centers;
	}
	centers = new float[3*c];
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
		mean += centers[j * 3];
		mean += centers[j * 3 + 1];
		mean += centers[j * 3 + 2];
	}
	mean /= (c * 3);
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
			mean += centers[j * 3];
			mean += centers[j * 3 + 1];
			mean += centers[j * 3 + 2];
		}
		mean /= (c * 3);
		if (abs(mean - mean_1) < e)
		{
			arret = true;
		}
	}
}