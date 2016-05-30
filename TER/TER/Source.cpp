#include <opencv2\opencv.hpp>
#include "FuzzyClustering.hpp"
#include "GrayFCM.h"
#include "RgbFCM.h"
#include "GrayCA.h"
#include "RobustCA.hpp"
#include "RgbCA.h"
#include <ctime>
#include "CA.hpp"
#include "Fuzzy.hpp"
using namespace cv;




int main(int argc, char* argv[])
{	
	/*
	Mat a = Mat(1024, 1024, CV_32F);
	float tmp = 0;
	clock_t tStart = clock();
	for (int k = 0; k < 200; k++)
	{
		for (int i = 0; i < 1024; i++)
		{
			const float* Mi = a.ptr<float>(i);
			for (int j = 0; j < 1024; j++)
			{
				tmp = Mi[j];
			}
		}
	}
	std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;

	tStart = clock();
	float** b = new float*[1024];
	for (int i = 0; i < 1024; i++)
	{
		b[i] = new float[1024];
	}
	for (int k = 0; k < 200; k++)
	{
		for (int i = 0; i < 1024; i++)
		{
			for (int j = 0; j < 1024; j++)
			{
				tmp = b[i][j];
			}
		}
	}
	
	std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;*/
	

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
	Fuzzy ca;
	ca.setImage(img);
	ca.exec();


	GrayCA gca;
	cvtColor(img, img, CV_BGR2GRAY);
	gca.setImage(img);
	gca.exec();
	//ca.CA(8, 2.f, 0.01, 50, 1.0, 1.0, 10.0, true, 5);
	//ca.exportMembership(std::string("ca/ca"));

	std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;

	/*
	tStart = clock();
	GrayCA gca;
	cvtColor(img, img, CV_BGR2GRAY);
	gca.setImage(img);

	gca.exec(8, 2.0, 0.01, 50, 1.1, 1.0, 10.f, true,5);
	std::cout << (double)(clock() - tStart) / CLOCKS_PER_SEC << std::endl;
	*/

	/*
	GrayFCM fcm;
	cvtColor(img, img, CV_BGR2GRAY);
	fcm.setImage(img);
	fcm.exec(2, 2.0, 0.01, 500, 0, 0, 16, 16);
	fcm.exportMembership(std::string("ca/ca"));
	*/
	
	/*
	GrayCA gca;
	cvtColor(img, img, CV_BGR2GRAY);
	gca.setImage(img);

	gca.exec(8, 2.0, 0.01, 50, 1.1, 1.0, 10.f,true);
	gca.exportMembership(std::string("ca/ca"));
	gca.exec();*/
	
	return 0;
}