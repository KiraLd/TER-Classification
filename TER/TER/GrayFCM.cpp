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
		for (int j = x_s;j < x_e;j++)
		{
			for (int k = y_s; k < y_e;k++)
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
		for (int j = x_s; j < x_e; j++)
		{
			for (int k = y_s; k < y_e;k++)
			{
				centers[i] += exp[i].at<float>(j, k)*img.at<uchar>(j, k);
				assert(!isnan(exp[i].at<float>(j, k)));
				sum += exp[i].at<float>(j, k);
			}
		}
		if (sum == 0)
		{
			centers[i] = 0;
		}
		else
		{
			centers[i] /= sum;
		}
	}
	delete[] exp;
}

void GrayFCM::dissimilarity()
{
	for (int i = 0; i < c; i++)
	{
		D[i] = Mat(x, y, CV_32F);
		for (int j = x_s; j < x_e; j++)
		{
			for (int k = y_s; k < y_e; k++)
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
	exp = (2.f / (m - 1));

	for (int i = 0; i < c; i++)
	{
		for (int j = x_s; j < x_e; j++)
		{
			for (int k = y_s; k < y_e; k++)
			{
				sum = 0.f;
				membership[i].at<float>(j, k) = 1.f;
				if (D[i].at<float>(j, k) == 0)
				{
					membership[i].at<float>(j, k) = 0;
				}
				else
				{
					for (int l = 0; l < c; l++)
					{
						if (D[l].at<float>(j, k) != 0)
						{
							sum += pow(D[i].at<float>(j, k) / (D[l].at<float>(j, k)), exp);
						}
						
					}
					if (sum != 0)
					{
						membership[i].at<float>(j, k) /= sum;
					}
					else
					{
						membership[i].at<float>(j, k) =0.f;
					}
				}
			}
		}
	}
	for (int i = x_s; i < x_e; i++)
	{
		for (int j = y_s; j < y_e; j++)
		{
			sum = 0.f;
			for (int k = 0; k < c;k++)
			{
				sum += membership[k].at<float>(i, j);
			}
			if (sum != 0)
			{
				for (int k = 0; k < c; k++)
				{
					membership[k].at<float>(i, j) = membership[k].at<float>(i, j) / sum;
					assert(!isnan(membership[k].at<float>(i, j)));
				}
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

}

void GrayFCM::exec(int c, float m, float e, int i_max,int x_start, int y_start, int x_end, int y_end)
{

	//cvtColor(img, img, CV_BGR2GRAY);
	x = img.rows;
	y = img.cols;
	x_s = x_start;
	y_s = y_start;
	x_e = x_end;
	y_e = y_end;
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
	float mean = c*256.f;
	float mean_1;
	while (!arret && i < i_max)
	{
		update_centers();
		dissimilarity();
		update_membership();
		i++;
		//std::cout << i << std::endl;
		mean_1 = mean;
		mean = 0;
		for (int j = 0; j < c;j++)
		{
			mean += centers[j];
		}
		mean /= (c);
	//	std::cout << mean << std::endl;
		if (abs(mean - mean_1) < e)
		{
			arret = true;
		}
	}
}