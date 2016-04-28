#include "RgbCA.h"

RgbCA::RgbCA() : RgbFCM()
{
	alive = NULL;
	Uexp = NULL;
	card = NULL;
}

RgbCA::~RgbCA()
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

void RgbCA::exportMembership(std::string file)
{
	std::string temp;
	std::ostringstream oss;
	Mat convert_;
	namedWindow("test");
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			oss << i;
			temp = file + oss.str() + std::string(".jpg");
			convert_ = convert(membership[i]);
			imwrite(temp, convert_);
			oss = std::ostringstream("");
		}
	}
}

Mat convert(Mat a)
{
	int x = a.rows;
	int y = a.cols;
	Mat b = Mat(x, y, CV_8U);
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			b.at<uchar>(i, j) = (uchar)a.at<float>(i, j) * 255;
		}
	}
	return b;
}

void RgbCA::init()
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

void RgbCA::getAlpha(int i)
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


float RgbCA::Nt(int i, int j)
{
	float sum1 = 0.0;
	float sum2 = 0.0;
	float inv_d = 0.0;
	for (int k = 0; k < c; k++)
	{
		if (alive[k])
		{
			inv_d = 1. / D[k].at<float>(i, j);
			sum1 += inv_d * card[k];
			sum2 += inv_d;
		}
	}
	return sum1 / sum2;
}

void RgbCA::update_cluster()
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

void RgbCA::exp()
{
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			pow(membership[i], m, Uexp[i]);
		}
	}
}

void RgbCA::update_centers()
{
	float sum;
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			centers[3 * i] = 0.f;
			centers[3 * i + 1] = 0.f;
			centers[3 * i + 2] = 0.f;
			sum = 0.f;
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y;k++)
				{
					centers[3*i] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[0];
					centers[3*i+1] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[1];
					centers[3*i+2] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[2];
					sum += Uexp[i].at<float>(j, k);
				}
			}
			centers[3 * i] /= sum;
			centers[3 * i + 1] /= sum;
			centers[3 * i + 2] /= sum;
			std::cout << centers[i] << std::endl;
		}
	}
}

void RgbCA::dissimilarity()
{
	float r, g, b;
	for (int i = 0; i < c; i++)
	{
		r = 0.f;
		g = 0.f;
		b = 0.f;
		if (alive[i])
		{
			D[i] = Mat(x, y, CV_32F);
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					r = img.at<Vec3b>(j, k)[0] - centers[3 * i];
					g = img.at<Vec3b>(j, k)[1] - centers[3 * i + 1];
					b = img.at<Vec3b>(j, k)[2] - centers[3 * i + 2];
					r *= r;
					g *= g;
					b *= b;
					D[i].at<float>(j, k) = sqrtf(r + g + b);
				}
			}
		}
	}
}

void RgbCA::update_membership()
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

void RgbCA::update_membership_bias()
{
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y;k++)
				{
					membership[i].at<float>(j, k) += (alpha / D[i].at<float>(j, k))*(card[i] - Nt(j, k));
				}
			}
		}
	}

	for (int i = 0; i < c; i++)
	{
		card[i] = 0.0;
		for (int j = 0; j < x; j++)
		{
			for (int k = 0; k < y; k++)
			{
				card[i] += membership[i].at<float>(j, k);
			}
		}
	}
}

void RgbCA::exec(int c, float m, float e, int i_max, float n0, float seuil)
{
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
	centers = new float[3*c];
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
	float mean = 0.f;
	for (int j = 0; j < c;j++)
	{
		mean += centers[j * 3];
		mean += centers[j * 3 + 1];
		mean += centers[j * 3 + 2];
	}
	mean /= (c * 3);
	float mean_1;
	bool arret = false;
	while (!arret && i < iter_max)
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
		mean_1 = mean;
		for (int j = 0; j < c;j++)
		{
			if (alive[j])
			{
				mean += centers[j * 3];
				mean += centers[j * 3 + 1];
				mean += centers[j * 3 + 2];
			}
		}
		mean /= (c * 3);
		if (abs(mean - mean_1) < e)
		{
			arret = true;
		}
	}
}