#include <opencv2\opencv.hpp>
#include "FuzzyClustering.hpp"
#include "GrayFCM.h"
#include "RgbFCM.h"
#include "GrayCA.h"
#include "RgbCA.h"
#include <ctime>
using namespace cv;


Mat* fcm_init(int x, int y, int c)
{
	Mat* membership = new Mat[c];
	RNG rng;
	for (int i = 0; i < c;i++)
	{
		membership[i] = Mat(x, y, CV_32F);
		for (int j = 0;j < x;j++)
		{
			for (int k = 0; k < y;k++)
			{
				membership[i].at<float>(j, k) = rng.uniform(0.f,1.f);
			}
		}
	}
	return membership;
}

void update_centers(const Mat& img, Mat* membership,float* centers, int c, float m)
{
	int x = img.rows;
	int y = img.cols;
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

void dissimilarity(const Mat& img, float* centers, int c, Mat* D)
{
	int x = img.rows;
	int y = img.cols;
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

void update_membership(Mat* D, Mat* U, int x, int y, int c, float m)
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
				U[i].at<float>(j,k) = std::powf(1. / D[i].at<float>(j, k), exp);
				for (int l = 0; l < c;l++)
				{
					sum += std::powf(1. / D[l].at<float>(j, k), exp);
				}
				U[i].at<float>(j, k) /= sum;
			}
		}
	}
}

Mat* fcm_NDG(const Mat& img, int c, float m, float e, int i_max)
{
	int x = img.rows;
	int y = img.cols;
	float* centers = new float[c];
	bool arret = false;
	Mat* D = new Mat[c];
	Mat* membership = fcm_init(x, y, c);
	int i = 0;
	while (!arret && i<i_max)
	{
		std::cout << i << std::endl;
		update_centers(img, membership, centers, c, m);
		dissimilarity(img, centers, c, D);
		update_membership(D, membership, x, y, c, m);
		i++;
	}
	delete[] D;
	delete[] centers;
	return membership;
}





int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		return -1;
	}
	std::string in = std::string(argv[1]);
	Mat img;
	img = imread(in, CV_LOAD_IMAGE_COLOR);
	if (!img.data)
	{
		std::cout << "Fichier introuvable" << std::endl;
		return -1;
	}

	clock_t tStart = clock();
	RgbCA ca;
	ca.setImage(img);
	int nb_c = 10;
	ca.exec(nb_c, 2.0, 0.01, 500, 1,50);
	std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;
	int s = 0;
	for (int i = 0; i < nb_c; i++)
	{
		namedWindow("test");
		if (ca.alive[i])
		{
			//imshow("test", ca.membership[i]);
			//waitKey(3000);
			s++;
		}
	}
	std::cout << s << std::endl;
	std::string file("test");
	ca.exportMembership(file);
	
	/*
	clock_t tStart = clock();
	GrayFCM fcm;
	fcm.setImage(img);
	fcm.exec(3, 2.0, 0.01, 10);
	std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;
	namedWindow("test");
	imshow("test", fcm.membership[0]);
	waitKey(3000);
	imshow("test", fcm.membership[1]);
	waitKey(3000);
	imshow("test", fcm.membership[2]);
	waitKey(3000);
	*/
	

	/*
	clock_t tStart = clock();
	RgbFCM fcm;
	fcm.setImage(img);
	fcm.exec(3, 2.0, 0.01, 25);
	std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;
	namedWindow("test");
	imshow("test", fcm.membership[0]);
	waitKey(3000);
	imshow("test", fcm.membership[1]);
	waitKey(3000);
	imshow("test", fcm.membership[2]);
	waitKey(3000);
	*/
	return 0;
}