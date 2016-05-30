#include "RgbCA.h"
#include <iostream>
#include <fstream>

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
	int j = 0;
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			oss << i;
			temp = file + oss.str() + std::string(".jpg");
			convert_ = convert(membership[i]);
			imwrite(temp, convert_);
			oss = std::ostringstream("");
			j++;
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
			b.at<uchar>(i, j) = (uchar)((fabs(a.at<float>(i, j)) * 255));
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
	if (i == 0)
	{
		alpha = 0;
	}
	else
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
		if (sum2 == 0)
		{
			alpha = 0.f;
		}
		else
		{
			alpha = n * sum1 / sum2;
		}
	}
	
}

void RgbCA::getAlphaAlt(int i)
{
	if (i != 0)
	{
		float cst = sqrt(x*y*n0);
		float n = n0*expf(-(float)i / cst);
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
		if (sum2 == 0)
		{
			alpha = 0.f;
		}
		else
		{
			alpha = n * sum1 / sum2;
		}
		
	}
	else
	{
		alpha = 0;
	}
	
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
			if (D[k].at<float>(i, j) != 0)
			{
				inv_d = 1. / D[k].at<float>(i, j);
				sum1 += inv_d * card[k];
				sum2 += inv_d;
			}
			
		}
	}
	if (sum2 == 0.f)
	{
		return 0.f;
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
	#pragma omp parallel for
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
	float* sum = new float[c];
	#pragma omp parallel for
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
					centers[3*i] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[0];
					centers[3*i+1] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[1];
					centers[3*i+2] += Uexp[i].at<float>(j, k)*img.at<Vec3b>(j, k)[2];
					sum[i] += Uexp[i].at<float>(j, k);
				}
			}
			if (sum[i] == 0)
			{
				centers[3 * i] = 0;
				centers[3 * i + 1] = 0;
				centers[3 * i + 2] = 0;
			}
			else
			{
				centers[3 * i] /= sum[i];
				centers[3 * i + 1] /= sum[i];
				centers[3 * i + 2] /= sum[i];
			}
			
		}
	}
	delete[] sum;
}

void RgbCA::dissimilarity()
{
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
					r[i] = img.at<Vec3b>(j, k)[0] - centers[3 * i];
					g[i] = img.at<Vec3b>(j, k)[1] - centers[3 * i + 1];
					b[i] = img.at<Vec3b>(j, k)[2] - centers[3 * i + 2];
					r[i] *= r[i];
					g[i] *= g[i];
					b[i] *= b[i];
					D[i].at<float>(j, k) = sqrtf(r[i] + g[i] + b[i]);
				}
			}
		}
	}
	delete[] r;
	delete[] g;
	delete[] b;
}

void RgbCA::update_membership()
{
	float sum;
	float exp;
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
					membership[i].at<float>(j, k) = 1.f;
					if (D[i].at<float>(j, k) == 0)
					{
						membership[i].at<float>(j, k) = 0;
					}
					else
					{
						for (int l = 0; l < c; l++)
						{

							if (alive[l] && D[l].at<float>(j, k) != 0)
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
							membership[i].at<float>(j, k) = 0.f;
						}
					}
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
					if (membership[i].at<float>(j, k) < 0.0)
					{
						membership[i].at<float>(j, k) = 0.0f;
					}
					else if (membership[i].at<float>(j, k) > 1.0)
					{
						membership[i].at<float>(j, k) = 1.0f;
					}
				}
			}
		}
	}
	float sum;
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			sum = 0.f;
			for (int k = 0; k < c;k++)
			{
				if (alive[k])
				{
					sum += membership[k].at<float>(i, j)*membership[k].at<float>(i, j);
				}
			}
			if (sum != 0)
			{
				for (int k = 0; k < c; k++)
				{
					if (alive[k])
					{
						membership[k].at<float>(i, j) = membership[k].at<float>(i, j)*membership[k].at<float>(i, j) / sum;
					}
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

void RgbCA::exec(int c, float m, float e, int i_max, float n0, float seuil,float t_max, bool alt,  RgbFCM* fcm, int it)
{
	x = img.rows;
	y = img.cols;
	this->c = c;
	this->m = m;
	this->e = e;
	this->n0 = n0;
	this->seuil = seuil;
	this->iter_max = i_max;
	this->t_max = t_max;
	if (fcm != NULL)
	{
		fcm->exec(c, m, e, it);
		centers = fcm->centers;
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
		membership = fcm->membership;
		for (int i = 0; i < c;i++)
		{
			alive[i] = true;
			card[i] = 0.0;
			for (int j = 0;j < x;j++)
			{
				for (int k = 0; k < y;k++)
				{
					card[i] += membership[i].at<float>(j, k);
				}
			}
		}
	}
	else
	{
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
		if (card != NULL)
		{
			delete[] card;
		}
		card = new float[c];
		Uexp = new Mat[c];
		init();
	}
	
	
	int i = 0;
	float mean = 255 * 3 * c;
	float mean_1;
	bool arret = false;
	int nb = 0;
	while (!arret && i < iter_max)
	{
		//std::cout << i << std::endl;
		exp();
		update_centers();
		dissimilarity();
		if (alt)
		{
			getAlphaAlt(i);
		}
		else
		{
			getAlpha(i);
		}
		
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
		mean /= (nb * 3);
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

void RgbCA::exec()
{
	int c = 10;
	float m = 2.0f;
	float e = 0.01;
	int iter_max = 50;
	float t_max = 500;
	float n0 = 1;
	float seuil = 5;
	std::string file("test");
	std::ostringstream oss;
	std::string temp;

	std::ofstream ofs;
	int b = 0;
	/*
	
	ofs.open("alt.txt");
	ofs << "Classe max\tClassique\tNb_iter\tAlt\tNb_iter\n";
	
	for (float i = 0.1; i < 20.0; i +=0.1)
	{
		exec(c, m, e, iter_max, i, seuil, t_max, false);
		ofs << i << "\t\t" << c_final << "\t" << iter;
		exec(c, m, e, iter_max, i, seuil, t_max, true);
		ofs<< "\t\t" << c_final << "\t" << iter << "\n";
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
			exec(c, m, e, iter_max, i, seuil, j, false);
			ofs << i <<"\t"<<j<< "\t" << c_final << "\t" << iter;
			exec(c, m, e, iter_max, i, seuil, j, true);
			ofs << "\t\t" << c_final << "\t" << iter << "\n";
			std::cout << b << std::endl;
			b++;
		}
	}
	ofs.close();*/

	//15 
	//Nombre de classe max

	/*
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
	for (float i = 1.1; i < 4.0; i += 0.1)
	{
		exec(c, i, e, iter_max, n0, seuil, t_max,false);
		ofs << i << "\t\t" << c_final << "\t\t\t" << iter;
		exec(c, i, e, iter_max, n0, seuil, t_max, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();
	*/
	/*
	ofs.open("Amplitude initiale.txt");
	ofs << "N0\tNb classe finale\tNb iter\tAlt\tNb_iter\n";
	std::cout << "Amplitude initiale: 15 tests" << std::endl;
	for (float i = 0.5; i < 10.0; i += 0.5)
	{
		exec(c, m, e, iter_max, i, seuil, t_max, false);
		ofs << i << "\t\t" << c_final << "\t\t\t" << iter;
		exec(c, m, e, iter_max, i, seuil, t_max, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();*/
	/*
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
	ofs.close();*/

	ofs.open("Timecst.txt");
	ofs << "Timecst\tNb classe finale\tNb iter\tAlt\tNb_iter\n";
	std::cout << "Constante de temps: 15 tests" << std::endl;
	for (float i = 1.0; i < 20000.0; i += 500.0)
	{
		exec(c, m, e, iter_max, n0,seuil, i, false);
		ofs << i << "\t\t" << c_final << "\t\t\t" << iter;
		exec(c, m, e, iter_max, n0, seuil, i, true);
		ofs << "\t\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		b++;
	}
	ofs.close();
	



	




}