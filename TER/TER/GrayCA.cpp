#include "GrayCA.h"
#include <iostream>
#include <fstream>

GrayCA::GrayCA(): GrayFCM()
{
	alive = NULL;
	Uexp = NULL;
	card = NULL;
}

Mat convertCA(Mat a)
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

void GrayCA::exportMembership(std::string file)
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
			convert_ = convertCA(membership[i]);
			imwrite(temp, convert_);
			oss = std::ostringstream("");
			j++;
		}
	}
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
		centers[i] = rng.uniform(0.f, 255.f);
		for (int j = 0;j < x;j++)
		{
			for (int k = 0; k < y;k++)
			{
				//centers[i] = rng.gaussian()
				membership[i].at<float>(j, k) = 0.f;
				card[i] += membership[i].at<float>(j, k);
			}
		}
	}
}

void GrayCA::initWindows()
{
	GrayFCM fcm;
	fcm.setImage(img);
	int c_x = x / 16;
	int c_y = y / 16;
	int bordure_x = x % 16;
	int bordure_y = y % 16;
	int x_16 = x - bordure_x - 16;
	int y_16 = y - bordure_y - 16;
	if (bordure_x != 0)
	{
		c_x++;
	}
	if (bordure_y != 0)
	{
		c_y++;
	}
	int c = 2 * c_x*c_y;
	if (membership != NULL)
	{
		delete[] membership;
	}
	membership = new Mat[c];
	if (centers != NULL)
	{
		delete[] centers;
	}
	centers = new float[c];
	int nb = 0;
	for (int i = 0; i < x_16; i +=16)
	{
		for (int j = 0; j < y_16; j += 16)
		{
			fcm.exec(2, 2.0, 0.01, 500, i, j, i + 16, j + 16);
			membership[nb] = Mat(fcm.membership[0]);
			centers[nb] = fcm.centers[0];
			nb++;
			membership[nb] = Mat(fcm.membership[1]);
			centers[nb] = fcm.centers[1];
		}
	}
}

void GrayCA::getAlpha(int i)
{
	if (i == 0)
	{
		alpha = 0;
	}
	else
	{
		
		//float time = sqrt(x)*sqrt(y)*n0;
		//std::cout << time << std::endl;
		//float n = (c*sqrt(x)*sqrt(y)/(x*y))*expf(-(float)i / time);
		
		float n = n0*expf(-(float)i / cst_time);
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


float GrayCA::Nt(int i,int j)
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
				inv_d = 1.f / D[k].at<float>(i, j);
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
			if (sum == 0)
			{
				centers[i] = 0;
			}
			else
			{
				centers[i] /= sum;
			}
					
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
	exp = (2.f / (m - 1));

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
					if (D[i].at<float>(j, k) != 0)
					{
						membership[i].at<float>(j, k) += (alpha / D[i].at<float>(j, k))*(card[i] - Nt(j, k));
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
					sum += membership[k].at<float>(i, j);
				}
			}
			if (sum != 0)
			{
				for (int k = 0; k < c; k++)
				{
					if (alive[k])
					{
						membership[k].at<float>(i, j) = membership[k].at<float>(i, j) / sum;
					}
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

void GrayCA::exec(int c, float m, float e, int i_max, float n0,float seuil,float cst, bool fcm, int iter)
{
	//cvtColor(img, img, CV_BGR2GRAY);
	x = img.rows;
	y = img.cols;
	this->c = c;
	this->m = m;
	this->e = e;
	this->n0 = n0;
	this->seuil = seuil;
	this->iter_max = i_max;
	cst_time = cst;
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
	float mean = c * 256.f;
	float mean_1;
	bool arret = false;
	int k;
	if (fcm == true)
	{
		GrayFCM gfcm;
		gfcm.setImage(img);
		gfcm.exec(c, m, e, iter,0,0,x,y);
		if (membership != NULL)
		{
			delete[] membership;
		}
		membership = new Mat[c];
		for (int i = 0; i < c; i++)
		{
			membership[i] = Mat(gfcm.membership[i]);
			centers[i] = gfcm.centers[i];
			alive[i] = true;
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
	else
	{
		init();
	}
	
	int i = 0;
	while (i < iter_max && !arret)
	{
		//std::cout << i << std::endl;
		
		dissimilarity();
		getAlpha(i);
		update_membership();
		update_membership_bias();
		update_cluster();
		exp();
		update_centers();
		i++;
		mean_1 = mean;
		mean = 0;
		k = 0;
		for (int j = 0; j < c; j++)
		{
			if (alive[j])
			{
				k++;
				mean += centers[i];
			}
		}
		if (k <= 2)
		{
			arret = true;
		}
		else
		{
			mean /= (float)k;
		}
		
		if (abs(mean - mean_1) < e)
		{
			arret = true;
		}
	}
	this->iter = i;
	c_final = 0;
	

	for (int j = 0; j < c; j++)
	{
		if (alive[j])
		{
			c_final++;
		}
	}
	
}

void GrayCA::exec()
{
	int c = 8;
	float m = 2.0f;
	float e = 0.01f;
	int iter_max = 50;
	float t_max = 10.f;
	float n0 = 5.f;
	float seuil = 5.f;
	std::ostringstream oss;
	std::ostringstream name;
	std::string temp;

	std::ofstream ofs;
	int b = 0;
	ofs.open("ca2/Classemax.txt");
	ofs << "C\tF\tIter\n";
	int x = img.rows;
	int y = img.cols;
	int size = x*y / 100;

	//~400 test

	for (int i = 5; i < 50; i++)
	{
		name << b;
		temp = std::string("ca2/c_") + name.str() + std::string("_");
		exec(i, m, e, iter_max, n0, seuil, t_max,true,5);
		ofs << i <<"\t"<< c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();


	ofs.open("ca2/fuzzyness.txt");
	ofs << "m\tF\tIter\n";

	for (float i = 1.1f; i < 5.f; i += 0.1f)
	{
		name << b;
		temp = std::string("ca2/m_") + name.str() + std::string("_");
		exec(c, i, e, iter_max, n0, seuil, t_max,true,5);
		ofs << i << "\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();


	ofs.open("ca2/n0.txt");
	ofs << "n0\tF\tIter\n";

	for (float i = 0.1f; i < 10.f; i += 0.1f)
	{
		name << b;
		temp = std::string("ca2/n0_") + name.str() + std::string("_");
		exec(c, m, e, iter_max, i, seuil, t_max,true,5);
		ofs << i << "\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();


	ofs.open("ca2/seuil.txt");
	ofs << "seuil\tF\tIter\n";
	for (float i = 1.f; i < 100.f; i += 1.f)
	{
		name << b;
		temp = std::string("ca2/seuil_") + name.str() + std::string("_");
		exec(c, m, e, iter_max, n0, i, t_max,true,5);
		ofs << i << "\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();
	
	ofs.open("ca2/cst.txt");
	ofs << "cst\tF\tIter\n";
	for (float i = 1.f; i < 100.f; i += 1.f)
	{
		name << b;
		temp = std::string("ca2/cst_") + name.str() + std::string("_");
		exec(c, m, e, iter_max, n0,seuil, i,true, 5);
		ofs << i << "\t" << c_final << "\t" << iter << "\n";
		std::cout << b << std::endl;
		exportMembership(temp);
		name = std::ostringstream("");
		b++;
	}
	ofs.close();
	
	ofs.open("ca2/test.txt");
	ofs << "C\t" << "M\t" << "N0\t" << "S\t" << "T\t" << "I\t" << "F\t\n";


	ofs.close();

	
	
}