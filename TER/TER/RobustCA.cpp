#include "RobustCA.hpp"
#include <iostream>
#include <fstream>

RobustCA::RobustCA() : RgbFCM()
{
	alive = NULL;
	Uexp = NULL;
	card = NULL;
	robust_card = NULL;
	robust_D = NULL;
	loss_D = NULL;
	T = nullptr;
	S = NULL;
}

RobustCA::~RobustCA()
{
	if (alive != NULL)
	{
		delete[] alive;
	}
	if (Uexp != NULL)
	{
		delete[] Uexp;
	}
	if (T != NULL)
	{
		delete[] T;
	}
	if (S != NULL)
	{
		delete[] S;
	}
	if (loss_D != NULL)
	{
		delete[] loss_D;
	}
	if (robust_card != NULL)
	{
		delete[] robust_card;
	}
	if (robust_D != NULL)
	{
		delete[] robust_D;
	}
}

void RobustCA::exportMembership(std::string file)
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
			convert_ = convertRobust(membership[i]);
			imwrite(temp, convert_);
			oss = std::ostringstream("");
			j++;
		}
	}
}

Mat convertRobust(Mat a)
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

void RobustCA::init()
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
		robust_card[i] = 0.0;
		robust_D[i] = Mat(x, y, CV_32F);
		for (int j = 0;j < x;j++)
		{
			for (int k = 0; k < y;k++)
			{
				membership[i].at<float>(j, k) = rng.uniform(0.f, 1.f);
				card[i] += membership[i].at<float>(j, k);
				robust_D[i].at<float>(j, k) = 1.f;
				robust_card[i] += membership[i].at<float>(j, k);
				
			}
		}
	}
}

void RobustCA::getAlpha(int i)
{
	if (i == 0)
	{
		alpha = 0.f;
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
						sum1 += membership[j].at<float>(k, l)*membership[j].at<float>(k, l)*loss_D[j].at<float>(k, l);
						sum3 += membership[j].at<float>(k, l)*robust_D[j].at<float>(k, l);
						assert(!isnan(membership[j].at<float>(k, l)));
						assert(!isnan(loss_D[j].at<float>(k, l)));
						assert(!isnan(robust_D[j].at<float>(k, l)));
					}
				}
				sum2 += sum3*sum3;
			}
		}
		assert(!isnan(sum1));
		assert(!isnan(sum2));
		if (sum2 > e)
		{
			alpha = n * sum1 / sum2;
			assert(!isnan(alpha));
		}
		else
		{
			alpha = 0.f;
		}
		
	}
	

}

void RobustCA::getAlphaAlt(int i)
{
	if (i == 0)
	{
		alpha = 0;
	}
	else
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
						sum1 += Uexp[j].at<float>(k, l)*loss_D[j].at<float>(k, l);
						sum3 += membership[j].at<float>(k, l)*robust_D[j].at<float>(k, l);;
						assert(!isnan(membership[j].at<float>(k, l)));
						assert(!isnan(loss_D[j].at<float>(k, l)));
						assert(!isnan(robust_D[j].at<float>(k, l)));
					}
				}
				sum2 += sum3*sum3;
			}
		}
		assert(!isnan(sum2));
		assert(!isnan(sum1));
		if (sum2 > e)
		{
			alpha = n * sum1 / sum2;
			assert(!isnan(alpha));
		}
		else
		{
			alpha = 0.f;
		}
	}
}


float RobustCA::Nt(int i, int j)
{
	float sum1 = 0.0;
	float sum2 = 0.0;
	float inv_d = 0.0;
	for (int k = 0; k < c; k++)
	{
		if (alive[k])
		{
			if (loss_D[k].at<float>(i, j) > e)
			{
				inv_d = 1.f / loss_D[k].at<float>(i, j);
				assert(!isnan(loss_D[k].at<float>(i, j)));
				assert(!isinf(loss_D[k].at<float>(i, j)));
				assert(!isnan(robust_card[k]));
				sum1 += inv_d * robust_card[k];
				sum2 += inv_d;
			}
		}
	}
	if (sum2 > e)
	{
		assert(!isnan(sum1));
		assert(!isnan(sum2));
		return sum1 / sum2;
	}
	else
	{
		return 0.f;
	}
	
}

void RobustCA::update_cluster()
{
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			if (robust_card[i] < seuil)
			{
				alive[i] = false;
			}
		}
	}
}

void RobustCA::exp()
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

void RobustCA::update_centers()
{
	float* sum = new float[c];
	RNG rng;
#pragma omp parallel for
	for (int i = 0; i < c;i++)
	{
		if (alive[i])
		{
			centers[i] = 0.f;
			sum[i] = 0.f;
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y;k++)
				{
					centers[i] += Uexp[i].at<float>(j, k)*(float)img.at<uchar>(j, k)*robust_D[i].at<float>(j,k);
					sum[i] += Uexp[i].at<float>(j, k)*robust_D[i].at<float>(j, k);
					assert(!isnan(Uexp[i].at<float>(j, k)));
					assert(!isnan(robust_D[i].at<float>(j, k)));
				}
			}
			if (sum[i] > e)
			{
				centers[i] /= sum[i];
				assert(!isnan(sum[i]));
				assert(!isnan(centers[i]));
				assert(!isinf(centers[i]));
			}
			else
			{
				centers[i] = rng.uniform(0.0, 255.0);
			}
		}
	}
	delete[] sum;
}

void RobustCA::dissimilarity()
{
	float* r = new float[c];
#pragma omp parallel for
	for (int i = 0; i < c; i++)
	{
		r[i] = 0.f;
		if (alive[i])
		{
			D[i] = Mat(x, y, CV_32F);
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y; k++)
				{
					r[i] = fabs((float)img.at<uchar>(j, k) - centers[i]);
					D[i].at<float>(j, k) = r[i];
				}
			}
		}
		
	}
	delete[] r;
}

void	RobustCA::wp()
{
	delta = delta - 1.f;
	delta = max(delta, 4.f);
	float d;
	float* K = new float[c];
	for (int k = 0; k < c; k++)
	{
		if (alive[k])
		{
			
			for (int i = 0; i < x; i++)
			{
				for (int j = 0; j < y; j++)
				{
					d = D[k].at<float>(i, j);
					if (d >= 0.f && d <= T[k])
					{
						robust_D[k].at<float>(i, j) = 0.f;
						if(T[k]>e)
						{ 
							robust_D[k].at<float>(i, j) = 1.f - (d*d) / (2.f * T[k] * T[k]);
							assert(!isnan(robust_D[k].at<float>(i, j)));
							assert(!isinf(robust_D[k].at<float>(i, j)));
							//std::cout <<1.f-((d*d)/(2.f*T[k]*T[k]))<< std::endl;
							//std::cout << robust_D[k].at<float>(i, j) << std::endl;
							assert(robust_D[k].at<float>(i, j) >= 0);
						}
						
					}
					else if (d <= T[k] + delta*S[k])
					{
						robust_D[k].at<float>(i, j) = 0.f;
						if (delta > e && S[k] > e)
						{
							robust_D[k].at<float>(i, j) = d - (T[k] + delta*S[k]);
							robust_D[k].at<float>(i, j) *= robust_D[k].at<float>(i, j);
							robust_D[k].at<float>(i, j) /= 2 * delta*delta*S[k];
							assert(robust_D[k].at<float>(j, k) >= 0);
							assert(!isnan(robust_D[k].at<float>(i, j)));
							assert(!isinf(robust_D[k].at<float>(i, j)));
						}
					}
					else
					{
						robust_D[k].at<float>(i, j) = 0.f;
					}
					assert(!isnan(robust_D[k].at<float>(i, j)));
					assert(!isinf(robust_D[k].at<float>(i, j)));
				}
			}
		}
	}
	float max = 0.f;
	float K_temp;
	for (int i = 0; i < c; i++)
	{
		max = 0.f;
		K_temp = 0.f;
		if (alive[i])
		{
			for (int j = 0; j < c; j++)
			{
				if (alive[j])
				{
					K_temp = (5.f * T[j] + delta*S[j]) / 6.f;
					if (K_temp > max)
					{
						max = K_temp;
					}
				}
			}
		}
		K[i] = max - ((5.f * T[i] + delta*S[i]) / 6.f);
		assert(!isnan(K[i]));
		assert(!isinf(K[i]));
	}

	for (int k = 0; k < c; k++)
	{
		if (alive[k])
		{
			loss_D[k] = Mat(x, y, CV_32F);
			for (int i = 0; i < x; i++)
			{
				for (int j = 0; j < y; j++)
				{
					d = D[k].at<float>(i, j);
					if (d >= 0.f && d <= T[k])
					{
						loss_D[k].at<float>(i, j) = 0;
						if (T[k] > e)
						{
							loss_D[k].at<float>(i, j) = d - ((d*d*d) / (6 * T[k] * T[k]));
						}
						assert(!isnan(loss_D[k].at<float>(i, j)));
						assert(!isinf(loss_D[k].at<float>(i, j)));
						
					}
					else if (d <= T[k] + delta*S[k])
					{
						loss_D[k].at<float>(i, j) = 0.f;
						if (delta > e && S[k] > e)
						{
							loss_D[k].at<float>(i, j) = d - (T[k] + delta*S[k]);
							loss_D[k].at<float>(i, j) *= loss_D[k].at<float>(i, j) * loss_D[k].at<float>(i, j);
							loss_D[k].at<float>(i, j) /= (6.f * delta*delta*S[k] * S[k]);
							loss_D[k].at<float>(i, j) += (5.f * T[k] + delta*S[k]) / 6.f;
							assert(!isnan(loss_D[k].at<float>(i, j)));
							assert(!isinf(loss_D[k].at<float>(i, j)));
						}
					}
					else
					{
						loss_D[k].at<float>(i, j) = ((5.f * T[k] + delta*S[k]) / 6.f) + K[k];
						assert(!isnan(loss_D[k].at<float>(i, j)));
						assert(!isinf(loss_D[k].at<float>(i, j)));
					}
					assert(!isnan(loss_D[k].at<float>(i, j)));
					assert(!isinf(loss_D[k].at<float>(i, j)));
				}
			}
		}
	}
	delete[] K;
}

void RobustCA::TS()
{
	Mat s(x, y, CV_32F);
	Mat* tri = new Mat[c];
	std::vector<float> vecFromMat;
	int size = x*y;
	float moyenne;
	int l;
	bool arret;
	float color;
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{ 
			vecFromMat = vector<float>();
			moyenne = 0.f;
			tri[i] = Mat(D[i]);
			//median 
			for (int j = 0; j < x;j++)
			{
				for (int k = 0; k < y; k++)
				{
					arret = false;
					l = 0;
					while (!arret && l < c)
					{
						if (alive[l])
						{
							if (l != i && D[i].at<float>(j, k) > D[l].at<float>(j, k))
							{
								arret = true;
							}
						}
						l++;
					}	
					if (!arret)
					{
						color = (float)img.at<uchar>(j,k);
						vecFromMat.push_back(color);
					}
				}
				

			}
			if (vecFromMat.size() <= 1)
			{
				T[i] = 0.f;
				S[i] = 0.f;
			}
			else
			{
				std::nth_element(vecFromMat.begin(), vecFromMat.begin() + vecFromMat.size() / 2, vecFromMat.end());
				T[i] = vecFromMat[vecFromMat.size() / 2];
				size = vecFromMat.size();
				for (int j = 0; j < size; j++)
				{
					vecFromMat[j] = fabs(vecFromMat[j] - T[i]);
				}
				std::nth_element(vecFromMat.begin(), vecFromMat.begin() + vecFromMat.size() / 2, vecFromMat.end());
				S[i] = vecFromMat[vecFromMat.size() / 2];
				assert(!isnan(T[i]));
				assert(!isinf(T[i]));
				assert(!isnan(S[i]));
				assert(!isinf(S[i]));
				vecFromMat.~vector();
			}

		}
	}
	delete[] tri;
}

void RobustCA::update_membership()
{
	float sum;

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
					
					if (loss_D[i].at<float>(j, k) > e)
					{
						membership[i].at<float>(j, k) = 1.f / loss_D[i].at<float>(j, k);
						for (int l = 0; l < c;l++)
						{
							if (alive[l] && loss_D[l].at<float>(j, k) > e)
							{
								sum += 1.f / loss_D[l].at<float>(j, k);
								assert(!isnan(loss_D[l].at<float>(j, k)));
								assert(!isinf(loss_D[l].at<float>(j, k)));
								assert(!isnan(sum));
								assert(!isinf(sum));
							}
						}
						assert(!isnan(sum));
						assert(!isinf(sum));
						membership[i].at<float>(j, k) /= sum;
						if (membership[i].at<float>(j, k) < 0.0)
						{
							membership[i].at<float>(j, k) = 0.0f;
						}
						else if (membership[i].at<float>(j, k) > 1.0)
						{
							membership[i].at<float>(j, k) = 1.0f;
						}
					}
					else
					{
						membership[i].at<float>(j, k) = 0.f;
					}
					
					//card[i] += membership[i].at<float>(j, k);
				}
			}
		}

	}
}

void RobustCA::update_membership_bias()
{
	for (int i = 0; i < c; i++)
	{
		if (alive[i])
		{
			for (int j = 0; j < x; j++)
			{
				for (int k = 0; k < y;k++)
				{
					if (loss_D[i].at<float>(j, k) > e)
					{
						membership[i].at<float>(j, k) += (alpha / loss_D[i].at<float>(j, k))*(robust_card[i] - Nt(j, k));
						assert(!isnan(alpha));
						assert(!isinf(alpha));
						assert(!isnan(Nt(j,k)));
						assert(!isinf(Nt(j, k)));
						assert(!isnan(loss_D[i].at<float>(j, k)));
						assert(!isinf(loss_D[i].at<float>(j, k)));
						assert(!isnan(robust_card[i]));
						assert(!isinf(robust_card[i]));
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
		robust_card[i] = 0.0f;
		for (int j = 0; j < x; j++)
		{
			for (int k = 0; k < y; k++)
			{
				card[i] += membership[i].at<float>(j, k);
				robust_card[i] += membership[i].at<float>(j, k)*robust_D[i].at<float>(j, k);
				assert(robust_D[i].at<float>(j, k) >= 0.f);
			}
		}
		//std::cout << "Card: " << i << ", " << card[i] << std::endl;
	}
}

void RobustCA::exec(int c, float m, float e, int i_max, float n0, float seuil, float t_max, bool alt, float delta)
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
	this->delta = delta;
	if(centers != NULL)
	{
		delete[] centers;
	}
	centers = new float[c];
	if (D != NULL)
	{
		delete[] D;
	}
	D = new Mat[c];
	if (robust_D != NULL)
	{
		delete[] robust_D;
	}
	robust_D = new Mat[c];
	if (robust_card != NULL)
	{
		delete[] robust_card;
	}
	robust_card = new float[c];
	if (loss_D != NULL)
	{
		delete[] loss_D;
	}
	loss_D = new Mat[c];
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
	if (T != NULL)
	{
		delete[] T;
	}
	T = new float[c];
	if (S != NULL)
	{
		delete[] S;
	}
	S = new float[c];
	init();


	int i = 0;
	float mean = c * 256.f;
	int classe = c;
	float mean_1;
	bool arret = false;
	while (!arret && i < iter_max)
	{
		
		exp();
		update_centers();
		dissimilarity();
		TS();
		wp();
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
		classe = 0;
		for (int j = 0; j < c;j++)
		{
			if (alive[j])
			{
				mean += centers[j];
				classe++;
			}
		}
		if (classe <= 2)
		{
			arret = true;
		}
		mean /= (float)(classe);
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
	//std::cout << "Classes finales: " << c_final << std::endl;
	iter = i;
}

void RobustCA::exec()
{
	int c = 10;
	float m = 2.0f;
	float e = 0.01f;
	int iter_max = 50;
	float t_max = 500.f;
	float n0 = 1.f;
	float seuil = 5.f;
	std::string file("test");
	std::ostringstream oss;
	std::ostringstream name;
	std::string temp;

	std::ofstream ofs;
	int b = 0;
	ofs.open("Robust.txt");
	ofs << "C\tn0\tseuil\tClasse\tI\n";
	int min, min_id;
	min_id = 0;
	min = 9999;
	for (int l = 8; l < 16; l++)
	{
		for (float h = 0.5f; h < 1.f;h += 0.1f)
		{
			for (float g = 0.1f; g < 1.f; g += 0.1f)
			{
				name << b;
				file = std::string("robust/test_") + name.str() + std::string("_");
				exec(l, m, e, iter_max, h, g, t_max, false, 12.f);
				ofs <<l<<"\t"<<h<<"\t"<<g<<"\t"<< "\t" << "\t" << c_final << "\t" << iter << "\n";
				if (c_final < min)
				{
					min = c_final;
					min_id = b;
				}
				
				std::cout << "Nombre de classe: " << c_final << std::endl;
				exportMembership(file);
				name = std::ostringstream("");
				b++;
				std::cout << b << std::endl;
			}
		}
	}
	ofs << "Classe min: " << min << " indice: " << min_id;
	ofs.close();
}
