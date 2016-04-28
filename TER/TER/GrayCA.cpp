#include "GrayCA.h"

GrayCA::GrayCA(): GrayFCM()
{
	alive = NULL;
	Uexp = NULL;
	card = NULL;
}

GrayCA::~GrayCA()
{
	if (alive != NULL)
	{
		delete[] alive;
	}
	if (Uexp != NULL)
	{
		delete[] Uexp;
	}
}

void GrayCA::init()
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
		alive[i] = true;
		card[i] = 0.0;
		for (int j = 0;j < x;j++)
		{
			for (int k = 0; k < y;k++)
			{
				membership[i].at<float>(j, k) = rng.uniform(0.f, 1.f);
				card[i] += membership[i].at<float>(j, k);
			}
		}
	}
}

void GrayCA::getAlpha(int i)
{
	float n = n0*expf(-(float)i / (float)iter_max);
	float sum1 = 0.0;
	float sum2 = 0.0;
	float sum3 = 0.0;
	for (int j = 0; j < c;j++)
	{
		if (alive[j])
		{
			sum3 = 0.0;
			for (int k = 0; k < x; k++)
			{
				for (int l = 0; l < y; l++)
				{
					sum1 += membership[j].at<float>(k, l)*membership[j].at<float>(k, l)*D[j].at<float>(k, l);
					sum3 += membership[j].at<float>(k, l);
				}
			}
			sum2 += sum3*sum3;
		}
	}
	alpha = n * sum1 / sum2;
}


float GrayCA::Nt(int i,int j)
{
	float sum1 = 0.0;
	float sum2 = 0.0;
	float inv_d = 0.0;
	for (int k = 0; k < c; k++)
	{
		if (alive[k])
		{
			inv_d = 1./D[k].at<float>(i, j);
			sum1 += inv_d * card[k];
			sum2 += inv_d;
		}
	}
	return sum1 / sum2;
}

void GrayCA::update_cluster()
{
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			if (card[i] < seuil)
			{
				alive[i] = false;
			}
		}
	}
}

void GrayCA::exp()
{
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			pow(membership[i], m, Uexp[i]);
		}
	}
}

void GrayCA::update_centers()
{
	float sum;
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			centers[i] = 0.f;
			sum = 0.f;
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y;k++)
				{
					centers[i] += Uexp[i].at<float>(j, k)*img.at<uchar>(j, k);
					sum += Uexp[i].at<float>(j, k);
				}
			}
			centers[i] /= sum;
			std::cout << centers[i] << std::endl;
		}
	}
}

void GrayCA::dissimilarity()
{
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
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
}

void GrayCA::update_membership()
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
		if (alive[i])
		{
			//card[i] = 0.0;
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					sum = 0.f;
					membership[i].at<float>(j, k) = std::powf(1. / D[i].at<float>(j, k), exp);
					for (int l = 0; l < c;l++)
					{
						if (alive[l])
						{
							sum += std::powf(1. / D[l].at<float>(j, k), exp);
						}
					}
					membership[i].at<float>(j, k) /= sum;
					//card[i] += membership[i].at<float>(j, k);
				}
			}
		}
		
	}
}

void GrayCA::update_membership_bias()
{
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y;k++)
				{
					membership[i].at<float>(j, k) += (alpha / D[i].at<float>(j, k))*(card[i] - Nt(j,k));
				}
			}
		}
	}

	for (int i = 0; i < c; i++)
	{
		card[i] = 0.0;
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					card[i] += membership[i].at<float>(j, k);
				}
			}
		}
		
	}
}

void GrayCA::exec(int c, float m, float e, int i_max, float n0,float seuil)
{
	cvtColor(img, img, CV_BGR2GRAY);
	x = img.rows;
	y = img.cols;
	this->c = c;
	this->m = m;
	this->e = e;
	this->n0 = n0;
	this->seuil = seuil;
	this->iter_max = i_max;
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
	if (alive != NULL)
	{
		delete[] alive;
	}
	alive = new bool[c];
	if (Uexp != NULL)
	{
		delete[] Uexp;
	}
	if (card != NULL)
	{
		delete[] card;
	}
	card = new float[c];
	Uexp = new Mat[c];
	init();
	int i = 0;
	while (i < iter_max)
	{
		std::cout << i << std::endl;
		exp();
		update_centers();
		dissimilarity();
		getAlpha(i);
		update_membership();
		update_membership_bias();
		update_cluster();
		i++;
	}
}