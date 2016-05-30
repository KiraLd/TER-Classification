#include "CA.hpp"
#include <iostream>
#include <fstream>
#include <ctime>

CA::CA() : RgbFCM()
{
	alive = NULL;
	Uexp = NULL;
	card = NULL;
	LookTable_SQRT = NULL;
	LookTable_D = NULL;
}

CA::~CA()
{
	if (alive != NULL)
	{
		delete[] alive;
	}
	if (Uexp != NULL)
	{
		delete[] Uexp;
	}
	if (LookTable_SQRT != NULL)
	{
		delete[] LookTable_SQRT;
	}
	if (LookTable_D != NULL)
	{
		for (int i = 0; i < 256; i++)
		{
			if (LookTable_D[i] != NULL)
			{
				delete[] LookTable_D[i];
			}
		}
		delete[] LookTable_D;
	}
}

void CA::exportMembership(std::string file)
{
	std::string temp;
	std::ostringstream oss;
	Mat convert_;
	int j = 0;
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			oss << i;
			temp = file + oss.str() + std::string(".jpg");
			convert_ = convert_CA(membership[i]);
			imwrite(temp, convert_);
			oss = std::ostringstream("");
			j++;
		}
	}
}

Mat convert_CA(Mat a)
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

void CA::init()
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
	clock_t tStart = clock();
	/*
	for (int i = 0; i < 10000;i++)
	{
		LookTable_U[i] = expf((float)i*0.01);
	}*/
	if (LookTable_D != NULL)
	{
		for (int i = 0; i < 256; i++)
		{
			delete[] LookTable_D[i];
		}
		delete[] LookTable_D;
	}
	LookTable_D = new float*[256];
	for (int i = 0; i < 256; i++)
	{
		LookTable_D[i] = new float[2560];
		for (int j = 0; j < 2560; j++)
		{
			LookTable_D[i][j] = i - j*0.1;
			LookTable_D[i][j] *= LookTable_D[i][j];
		}
	}
	int n = 256 * 256 * 3*10;
	if (LookTable_SQRT != NULL)
	{
		delete[] LookTable_SQRT;
	}
	LookTable_SQRT = new float[n];
	for (int i = 0; i < n; i++)
	{
		LookTable_SQRT[i] = sqrtf(i*0.1);
	}
	std::cout << LookTable_SQRT[n-1] << std::endl;
	std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;
}

void CA::getAlpha(int i)
{
	float n = n0*expf(-(float)i / t_max);
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
					sum1 += Uexp[j].at<float>(k, l)*D[j].at<float>(k, l);
					sum3 += membership[j].at<float>(k, l);
				}
			}
			sum2 += sum3*sum3;
		}
	}
	sum2 += e;
	alpha = n * sum1 / sum2;

}


float CA::Nt(int i, int j)
{
	float sum1 = 0.0;
	float sum2 = 0.0;
	float inv_d = 0.0;
	for (int k = 0; k < c; k++)
	{
		if (alive[k])
		{
			inv_d = 1. / (D[k].at<float>(i, j)+e);
			sum1 += inv_d * card[k];
			sum2 += inv_d;
		}
	}
	return sum1 / sum2;
}

void CA::update_cluster()
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

void CA::exp()
{
//#pragma omp parallel for
	for (int i = 0; i < c; i++)
	{
		for (int j = 0; j < x; j++)
		{
			for (int k = 0; k < y; k++)
			{
				assert(!isnan(membership[i].at<float>(j,k)));
			}
		}
		
		if (alive[i])
		{
			pow(membership[i], m, Uexp[i]);
		}
	}
}

void CA::update_centers()
{
	float* sum = new float[c];
	float s;
	float a, b, d;
//#pragma omp parallel for
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			centers[3 * i] = 0.f;
			centers[3 * i + 1] = 0.f;
			centers[3 * i + 2] = 0.f;
			sum[i] = 0.f;
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y;k++)
				{
					centers[3 * i] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[0];
					centers[3 * i + 1] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[1];
					centers[3 * i + 2] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[2];
					sum[i] += Uexp[i].at<float>(j, k);
				}
			}
			sum[i] += e;
			
			centers[3 * i] /= sum[i];
			centers[3 * i + 1] /= sum[i];
			centers[3 * i + 2] /= sum[i];
			a = centers[3 * i];
			b = centers[3 * i+1];
			d = centers[3 * i+2];
			s = sum[i];
		}
	}
	delete[] sum;
}

void CA::dissimilarity()
{
	clock_t tStart = clock();
	float* r = new float[c];
	float* g = new float[c];
	float* b = new float[c];
#pragma omp parallel for
	for (int i = 0; i < c; i++)
	{
		r[i] = 0.f;
		g[i] = 0.f;
		b[i] = 0.f;
		if (alive[i])
		{
			D[i] = Mat(x, y, CV_32F);
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					r[i] = LookTable_D[img.at<Vec3b>(j, k)[0]][int(10*centers[3 * i])];
					g[i] = LookTable_D[img.at<Vec3b>(j, k)[1]][int(10*centers[3 * i+1])];
					b[i] = LookTable_D[img.at<Vec3b>(j, k)[2]][int(10*centers[3 * i+2])];
					D[i].at<float>(j, k) = LookTable_SQRT[(int)(10*(r[i] + g[i] + b[i]))];
				}
			}
		}
	}
	delete[] r;
	delete[] g;
	delete[] b;
	//std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;
}

void CA::update_membership()
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
					assert(!isnan(membership[i].at<float>(j, k)));
					for (int l = 0; l < c;l++)
					{
						if (alive[l])
						{
							sum += std::powf(1. / (D[l].at<float>(j, k)+e), exp);
						}
					}
					membership[i].at<float>(j, k) /= sum;
					if (membership[i].at<float>(j, k) < 0.0)
					{
						membership[i].at<float>(j, k) = 0.0f;
					}
					else if (membership[i].at<float>(j, k) > 1.0)
					{
						membership[i].at<float>(j, k) = 1.0f;
					}
					assert(!isnan(membership[i].at<float>(j, k)));
					//card[i] += membership[i].at<float>(j, k);
				}
			}
		}

	}
}

void CA::update_membership_bias()
{
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y;k++)
				{
					assert(!isnan(membership[i].at<float>(j, k)));
					membership[i].at<float>(j, k) += (alpha / (D[i].at<float>(j, k)+e))*(card[i] - Nt(j, k));
					assert(!isnan(alpha));
					assert(!isnan(Nt(j,k)));
					assert(!isnan(membership[i].at<float>(j, k)));
					if (membership[i].at<float>(j, k) < 0.0)
					{
						membership[i].at<float>(j, k) = 0.0f;
					}
					else if (membership[i].at<float>(j, k) > 1.0)
					{
						membership[i].at<float>(j, k) = 1.0f;
					}
					assert(!isnan(membership[i].at<float>(j, k)));
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

void CA::exec(int c, float m, float e, int i_max, float n0, float seuil, float t_max, bool alt)
{
	x = img.rows;
	y = img.cols;
	if (alt)
	{
		this->t_max= sqrt(x*y*n0);
	}
	else
	{
		this->t_max = t_max;
	}
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
	centers = new float[3 * c];
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
	Uexp = new Mat[c];
	if (card != NULL)
	{
		delete[] card;
	}
	card = new float[c];
	
	
	init();


	int i = 0;
	float mean = mean = 255 * 3 * c;
	float mean_1;
	bool arret = false;
	int nb = 0;
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
		mean = 0.f;
		nb = 0;
		
		for (int j = 0; j < c;j++)
		{
			if (alive[j])
			{
				mean += centers[j * 3];
				mean += centers[j * 3 + 1];
				mean += centers[j * 3 + 2];
				nb++;
			}
		}
		std::cout << mean << std::endl;
		if (abs(mean - mean_1) < e)
		{
			arret = true;
		}

	}
	c_final = 0;
	for (int j = 0; j < c; j++)
	{
		if (alive[j])
		{
			c_final++;
		}
	}
	iter = i;
}

void CA::exec()
{
	int c = 10;
	float m = 2.0f;
	float e = 0.01;
	int iter_max = 5000;
	float t_max = 500;
	float n0 = 1;
	float seuil = 5;
	std::string file("test");
	std::ostringstream oss;
	std::string temp;



	std::ofstream ofs;
	ofs.open("alt.txt");
	ofs << "Classe max\tClassique\tNb_iter\tAlt\tNb_iter\n";
	int b = 0;
	for (float i = 1.0; i < 20.0; i += 0.1)
	{
		exec(c, m, e, iter_max, i, seuil, t_max, false);
		ofs << i << "\t\t" << c_final << "\t" << iter;
		exec(c, m, e, iter_max, i, seuil, t_max, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << "Test: " << i - 4 << std::endl;
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();


	ofs.open("alt2.txt");
	ofs << "n0\tcst\tfinal\titer\tAlt\tNb_iter\n";

	for (float i = 1.0; i < 20.0; i += 1.f)
	{
		for (float j = 5.0; j < 200; j++)
		{
			exec(c, m, e, iter_max, i, seuil, j, true);
			ofs << i << "\t" << j << "\t" << c_final << "\t" << iter;
			exec(c, m, e, iter_max, i, seuil, j, true);
			ofs << "\t\t" << c_final << "\t" << iter << "\n";
			std::cout << b << std::endl;
			b++;
		}
	}
	ofs.close();

	//15 
	//Nombre de classe max


	ofs.open("nbclasse.txt");
	ofs << "Classe max\tNbFinal\tNbIter\tAlt\tNb_iter\n";
	std::cout << "Nombre de classe max: 15 tests" << std::endl;
	for (int i = 5; i < 15; i++)
	{
		exec(i, m, e, iter_max, n0, seuil, t_max, false);
		ofs << i << "\t\t" << c_final << "\t" << iter;
		exec(i, m, e, iter_max, n0, seuil, t_max, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();

	ofs.open("Fuzzyness.txt");
	ofs << "Fuzzyness\tNb classe finale\tNb iter\tAlt\tNb_iter\n";
	std::cout << "Coef flou: 15 tests" << std::endl;
	for (float i = 2.0; i < 3.0; i += 0.1)
	{
		exec(c, i, e, iter_max, n0, seuil, t_max, false);
		ofs << i << "\t\t" << c_final << "\t\t\t" << iter;
		exec(c, i, e, iter_max, n0, seuil, t_max, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();

	ofs.open("Amplitude initiale.txt");
	ofs << "N0\tNb classe finale\tNb iter\tAlt\tNb_iter\n";
	std::cout << "Amplitude initiale: 15 tests" << std::endl;
	for (float i = 0.5; i < 2.0; i += 0.1)
	{
		exec(c, m, e, iter_max, i, seuil, t_max, false);
		ofs << i << "\t\t" << c_final << "\t\t\t" << iter;
		exec(c, m, e, iter_max, i, seuil, t_max, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();

	ofs.open("Seuil.txt");
	ofs << "Seuil\tNb classe finale\tNb iter\tAlt\tNb_iter\n";
	std::cout << "Seuil d'agglomeration: 12 tests" << std::endl;
	for (float i = 5.0; i < 600; i += 50.0)
	{
		exec(c, m, e, iter_max, n0, i, t_max, false);
		ofs << i << "\t\t" << c_final << "\t\t\t" << iter;
		exec(c, m, e, iter_max, n0, i, t_max, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();

	ofs.open("Timecst.txt");
	ofs << "Timecst\tNb classe finale\tNb iter\tAlt\tNb_iter\n";
	std::cout << "Constante de temps: 15 tests" << std::endl;
	for (float i = 5.0; i < 60; i += 5.0)
	{
		exec(c, m, e, iter_max, n0, seuil, i, false);
		ofs << i << "\t\t" << c_final << "\t\t\t" << iter;
		exec(c, m, e, iter_max, n0, seuil, i, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();
}