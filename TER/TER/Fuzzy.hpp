#pragma once
#ifndef _FUZZY_
#define _FUZZY_

#include <opencv2\opencv.hpp>
using namespace cv;

struct Fuzzy
{
	float** img;
	float*** membership;
	float* centers;
	float*** distance;
	float* cardinal;
	bool* alive;
	float m, e, n0, seuil, time,alpha;
	int x, y, c,iter,c_final,iter_max;
	bool fcm;
	

	Fuzzy();
	~Fuzzy();
	virtual void initFCM();
	virtual void initCA();
	virtual void centersFCM();
	virtual void centersCA();
	virtual void distance_euclidienne();
	virtual void updateU_FCM();
	virtual void updateU_CA();
	virtual void getAlpha(int i);
	virtual float Nt(int i, int j);
	virtual void setImage(const Mat& img);
	virtual void exportMembership(std::string file);
	virtual void FCM(int c, float m, float e, int iter);
	virtual void CA(int c, float m, float e, int iter, float n0, float seuil, float cst, bool fcm, int fcm_iter);
	virtual void CA(int c, float e, int iter, bool fcm, int fcm_iter);
	virtual void exec();


};


#endif
