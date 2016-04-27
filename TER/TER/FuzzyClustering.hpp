#pragma once
#ifndef _FUZZY_C_
#define _FUZZY_C_
#include <opencv2\opencv.hpp>
using namespace cv;

struct FuzzyClustering
{
	Mat img;
	int x;
	int y;
	int c;
	float m;
	float e;
	Mat* membership;
	float* centers;
	
	FuzzyClustering();
	~FuzzyClustering();
	void setImage(const Mat& img);
	void setImage(std::string file, int type);
	void exportMembership(std::string file);
	void show(int time);
	virtual void exec(int c, float m, float e, int i_max) = 0;
};

#endif
