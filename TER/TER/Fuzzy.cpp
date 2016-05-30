#include "Fuzzy.hpp"
#include <iostream>
#include <fstream>

Fuzzy::Fuzzy()
{
	img = nullptr;
	membership = nullptr;
	centers = nullptr;
	distance = nullptr;
	cardinal = nullptr;
	alive = nullptr;
	fcm = true;
}

Fuzzy::~Fuzzy()
{
	
	if (img != nullptr)
	{
		for (int i = 0; i < x; i++)
		{
			delete[] img[i];
		}
		delete[] img;
	}
	if (membership != nullptr)
	{
		for (int i = 0; i < c; i++)
		{
			for (int j = 0; j < x; j++)
			{
				delete[] membership[i][j];
			}
			delete[] membership[i];
		}
		delete[] membership;
	}
	if (centers != nullptr)
	{
		delete[] centers;
	}
	if (distance != nullptr)
	{
		for (int i = 0; i < c; i++)
		{
			for (int j = 0; j < x; j++)
			{
				delete[] distance[i][j];
			}
			delete[] distance[i];
		}
		delete[] distance;
	}
	if (cardinal != nullptr)
	{
		delete[] cardinal;
	}

	if (alive != nullptr)
	{
		delete[] alive;
	}
	
}

void Fuzzy::initFCM()
{
	RNG rng;
	if (membership != nullptr)
	{
		for (int i = 0; i < c; i++)
		{
			for (int j = 0; j < x; j++)
			{
				//delete[] membership[i][j];
			}
			//delete[] membership[i];
		}
		delete[] membership;
	}
	if (alive != nullptr)
	{
		delete[] alive;
	}
	alive = new bool[c];
	membership = new float**[c];
	for (int i = 0; i < c; i++)
	{
		membership[i] = new float*[x];
		alive[i] = true;
		for (int j = 0; j < x; j++)
		{
			membership[i][j] = new float[y];
			for (int k = 0; k < y; k++)
			{
				membership[i][j][k] = rng.uniform(0.f, 1.f);
			}
		}
	}

}

void Fuzzy::initCA()
{
	if (fcm)
	{
		if (cardinal != nullptr)
		{
			delete[] cardinal;
		}
		
		cardinal = new float[c];
		
		for (int i = 0; i < c; i++)
		{
			cardinal[i] = 0.f;
			alive[i] = true;
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					cardinal[i] += membership[i][j][k];
				}
			}
		}
	}
	else
	{
		if (centers != nullptr)
		{
			delete[] centers;
		}
		centers = new float[c];
		if (distance != nullptr)
		{
			for (int i = 0; i < c; i++)
			{
				for (int j = 0; j < x; j++)
				{
					//delete[] distance[i][j];
				}
				//delete[] distance[i];
			}
			delete[] distance;
		}
		distance = new float**[c];
		for (int i = 0; i < c; i++)
		{
			distance[i] = new float*[x];
			for (int j = 0; j < x; j++)
			{
				distance[i][j] = new float[y];
			}
		}
		RNG rng;
		if (membership != nullptr)
		{
			for (int i = 0; i < c; i++)
			{
				for (int j = 0; j < x; j++)
				{
					//delete[] membership[i][j];
				}
				//delete[] membership[i];
			}
			delete[] membership;
		}
		if (alive != nullptr)
		{
			delete[] alive;
		}
		alive = new bool[c];
		if (cardinal != nullptr)
		{
			delete[] cardinal;
		}

		cardinal = new float[c];

		membership = new float**[c];
		for (int i = 0; i < c; i++)
		{
			membership[i] = new float*[x];
			alive[i] = true;
			centers[i] = rng.uniform(0.f, 255.f);
			for (int j = 0; j < x; j++)
			{
				membership[i][j] = new float[y];
			}
		}
	}
}

void Fuzzy::centersFCM()
{
	float sum;
	float exp;
	for (int i = 0; i < c; i++)
	{
		centers[i] = 0.f;
		sum = 0.f;
		for (int j = 0; j < x; j++)
		{
			for (int k = 0; k < y; k++)
			{
				exp = membership[i][j][k];
				exp *= exp;
				centers[i] += exp * img[j][k];
				sum += exp;
			}
		}
		if (sum == 0)
		{
			centers[i] = 0.f;
		}
		else
		{
			centers[i] /= sum;
		}
		assert(!isnan(centers[i]));
	}
	
}

void Fuzzy::centersCA()
{
	float sum;
	float Us;
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			centers[i] = 0.f;
			sum = 0.f;
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					Us = membership[i][j][k];
					Us *= Us;
					centers[i] += Us * img[j][k];
					sum += Us;
				}
			}
			if (sum == 0.f)
			{
				centers[i] = 0.f;
			}
			else
			{
				centers[i] /= sum;
			}
		}
	}
}

void Fuzzy::distance_euclidienne()
{
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					distance[i][j][k] = abs(img[j][k] - centers[i]);
				}
			}
		}
	}
}

void Fuzzy::updateU_FCM()
{
	float sum;
	float exp = 2.f / (m - 1.f);
	float D1, D2;
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					sum = 0.f;
					membership[i][j][k] = 1.f;
					D1 = distance[i][j][k];
					if (D1 == 0.f)
					{
						membership[i][j][k] = 0.f;
					}
					else
					{
						for (int l = 0; l < c; l++)
						{
							D2 = distance[l][j][k];
							if (D2 != 0)
							{
								sum += pow(D1 / D2, exp);
							}
						}
						if (sum != 0.f)
						{
							membership[i][j][k] /= sum;
							if (membership[i][j][k] > 1.f)
							{
								membership[i][j][k] = 1.f;
							}
							else if (membership[i][j][k] < 0.f)
							{
								membership[i][j][k] = 0.f;
							}
						}
						else
						{
							membership[i][j][k] = 0.f;
						}
					}
				}
			}
		}
	}

	if (!fcm)
	{
		for (int i = 0; i < x; i++)
		{
			for (int j = 0; j < y; j++)
			{
				sum = 0.f;
				for (int k = 0; k < c; k++)
				{
					sum += membership[k][i][j];
				}
				if (sum != 0)
				{
					for (int k = 0; k < c; k++)
					{
						membership[k][i][j] = membership[k][i][j] / sum;
					}
				}
			}
		}
	}
}

void Fuzzy::updateU_CA()
{
	float d;
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					d = distance[i][j][k];
					if (d != 0.f)
					{
						membership[i][j][k] += (alpha / d)*(cardinal[i] - Nt(j, k));
					}
					if (membership[i][j][k] > 1.f)
					{
						membership[i][j][k] = 1.f;
					}
					else if (membership[i][j][k] < 0.f)
					{
						membership[i][j][k] = 0.f;
					}
				}
			}
		}
	}
	/*
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			d = 0.f;
			for (int k = 0; k < c; k++)
			{
				if (alive[k])
				{
					d += membership[k][i][j];
				}
			}
			if (d != 0.f)
			{
				for (int k = 0; k < c; k++)
				{
					if (alive[k])
					{
						membership[k][i][j] = membership[k][i][j] / d;
					}
				}
			}
		}
	}*/

	for (int i = 0; i < c; i++)
	{
		cardinal[i] = 0.f;
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					cardinal[i] += membership[i][j][k];
				}
			}
		}
	}
}

float Fuzzy::Nt(int i, int j)
{
	float sum1 = 0.f;
	float sum2 = 0.f;
	float inv_d = 0.f;
	float d;
	for (int k = 0; k < c; k++)
	{
		if (alive[k])
		{
			d = distance[k][i][j];
			if (d != 0)
			{
				inv_d = 1.f / d;
				sum1 += inv_d * cardinal[k];
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

void Fuzzy::setImage(const Mat& img)
{
	Mat img_gray;
	cvtColor(img, img_gray, CV_BGR2GRAY);
	
	if (this->img != nullptr)
	{
		for (int i = 0; i < x;i++)
		{
			delete[] this->img[i];
		}
		delete[] this->img;
	}
	x = img_gray.rows;
	y = img_gray.cols;
	this->img = new float*[x];
	for (int i = 0; i < x; i++)
	{
		this->img[i] = new float[y];
		for (int j = 0; j < y; j++)
		{
			this->img[i][j] = (float)img_gray.at<uchar>(i, j);
		}
	}
}

void Fuzzy::exportMembership(std::string file)
{
	std::string temp;
	std::ostringstream oss;
	Mat convert = Mat(x, y, CV_8U);
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			oss << i;
			temp = file + oss.str() + std::string(".jpg");
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					convert.at<uchar>(j, k) = (uchar)(fabs(membership[i][j][k]) * 255);
				}
			}
			imwrite(temp, convert);
			oss = std::ostringstream("");
		}
		
	}
}

void Fuzzy::FCM(int c, float m, float e, int iter)
{
	this->c = c;
	this->m = m;
	this->e = e;
	this->iter_max = iter;
	if (centers != nullptr)
	{
		delete[] centers;
	}
	centers = new float[c];
	if (distance != nullptr)
	{
		for (int i = 0; i < c; i++)
		{
			for (int j = 0; j < x; j++)
			{
				//delete[] distance[i][j];
			}
			//delete[] distance[i];
		}
		delete[] distance;
	}
	distance = new float**[c];
	for (int i = 0; i < c; i++)
	{
		distance[i] = new float*[x];
		for (int j = 0; j < x; j++)
		{
			distance[i][j] = new float[y];
		}
	}

	initFCM();
	bool arret = false;
	float C, Cprev;
	C = c*256.f;
	int i = 0;
	while (!arret && i < iter_max)
	{
		centersFCM();
		distance_euclidienne();
		updateU_FCM();
		Cprev = C;
		C = 0.f;
		for (int j = 0; j < c; j++)
		{
			C += centers[j];
		}
		C /= c;
		if (abs(C - Cprev) < e)
		{
			arret = true;
		}
		i++;
	}
	this->iter = i;

}

void Fuzzy::getAlpha(int i)
{
	if (i == 0)
	{
		alpha = 0;
	}
	else
	{
		float n = n0 * expf(-(float)i / time);
		float sum1 = 0.f;
		float sum2 = 0.f;
		float sum3 = 0.f;
		float U;
		for (int j = 0; j < c; j++)
		{
			if (alive[j])
			{
				sum3 = 0.f;
				for (int k = 0; k < x; k++)
				{
					for (int l = 0; l < y; l++)
					{
						U = membership[j][k][l];
						sum1 += U*U*distance[j][k][l];
						sum3 += U;
					}
				}
				sum2 += sum3*sum3;
			}
		}
		if (sum2 == 0.f)
		{
			alpha = 0.f;
		}
		else
		{
			alpha = n*sum1 / sum2;
		}
	}
}

void Fuzzy::CA(int c, float m, float e, int iter, float n0, float seuil, float cst, bool fcm, int fcm_iter)
{
	this->c = c;
	this->m = m;
	this->e = e;
	this->iter_max = iter;
	this->n0 = n0;
	this->seuil = seuil;
	this->time = cst;
	this->fcm = fcm;
	if (fcm)
	{
		this->fcm = false;
		FCM(c, m, e, fcm_iter);
		this->fcm = true;
	}
	initCA();
	int i = 0;
	float C = 255.f * c;
	float Cprev;
	int c_alt = c;
	bool arret = false;
	bool stop;
	while (i < iter_max && !arret)
	{
		distance_euclidienne();
		getAlpha(i);
		updateU_FCM();
		updateU_CA();
		centersCA();
		
		stop = true;
		for (int j = 0; j < c; j++)
		{
			if (alive[j])
			{
				if (cardinal[j] < seuil)
				{
					alive[j] = false;
				}
			}
		}
		Cprev = C;
		C = 0;
		c_alt = 0;
		for (int j = 0; j < c; j++)
		{
			if (alive[j])
			{
				C += centers[j];
				c_alt++;
			}
		}

		if (c_alt < 1)
		{
			arret = true;
		}
		C /= (float)c_alt;
		if (fabs(C - Cprev) < e)
		{
			arret = true;
		}
		i++;
	}
	this->iter = i;
	c_final = c_alt;
}

void Fuzzy::CA(int c, float e, int iter, bool fcm, int fcm_iter)
{
	CA(c, 2.0f, e, iter, 1.f, sqrt(x*y), sqrt(x*y), fcm, fcm_iter);
}

void Fuzzy::exec()
{
	int c = 8;
	float m = 2.0f;
	float e = 0.01f;
	int iter_max = 50;
	float t_max = 10.f;
	float seuil = (2.f / 3.f)*sqrt(x*y);
	float n0 = 5.f;
	std::ostringstream oss;
	std::ostringstream name;
	std::string temp;

	std::ofstream ofs;
	int b = 0;
	
	ofs.open("ca/Classemax.txt");
	ofs << "C\tF\tIter\n";

	int size = x*y / 100;

	//~400 test

	
	for (int i = 5; i < 20; i++)
	{
		name << b;
		temp = std::string("ca/c_") + name.str() + std::string("_");
		CA(i, m, e, iter_max, n0, seuil,t_max, true, 5);
		ofs << i << "\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();


	ofs.open("ca/n0.txt");
	ofs << "n0\tF\tIter\n";

	for (float i = 0.1f; i < 10.f; i += 0.1f)
	{
		name << b;
		temp = std::string("ca/n0_") + name.str() + std::string("_");
		CA(c, m, e, iter_max, i, seuil, t_max, true, 5);
		ofs << i << "\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();


	ofs.open("ca/seuil.txt");
	ofs << "seuil\tF\tIter\n";
	for (float i = 1.f; i < 100.f; i += 1.f)
	{
		name << b;
		temp = std::string("ca/seuil_") + name.str() + std::string("_");
		CA(c, m, e, iter_max, n0, i, t_max, true, 5);
		ofs << i << "\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();

	ofs.open("ca/cst.txt");
	ofs << "cst\tF\tIter\n";
	for (float i = 1.f; i < 100.f; i += 1.f)
	{
		name << b;
		temp = std::string("ca/cst_") + name.str() + std::string("_");
		CA(c, m, e, iter_max, n0, seuil, i, true, 5);
		ofs << i << "\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();

	std::cout << "Nombre d'essai: " << std::endl;
	int n = (20 - 5)* 9*9;
	std::cout << n << std::endl;
	ofs.open("ca/test.txt");
	ofs << "C\t"<< "N0\t" << "S\t" << "T\t" << "I\t" << "F\t\n";
	//c
	for (int i = 6; i < 20; i++)
	{
		//n0
		for (float j = 1.f; j < 1.01f; j += 0.1f)
		{
			//seuil
			for (float k = 1.f; k < 10.f; k += 1.f)
			{
				//cst
				for (float l = 1.f; l < 10.f; l += 1.f)
				{
					CA(i, 2.0, e, iter_max, j, k, l, false, 5);
					if (c_final == 5)
					{
						ofs << i << "\t" << j << "\t" << k << "\t" << l << "\t" << c_final << "\t" << iter << "\n";
					}
					std::cout<<b<<std::endl;
					b++;
				}
			}
		}
	}

	ofs.close();



}