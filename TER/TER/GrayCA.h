#pragma once
#ifndef _GRAY_CA_
#define _GRAY_CA_
#include "GrayFCM.h"

struct GrayCA: public GrayFCM
{
	bool* alive;
	Mat* Uexp;
	float* card;
	int iter_max;
	float cst_time;
	int iter;
	int c_final;
	float n0;
	float alpha;
	float seuil;
	GrayCA();
	~GrayCA();
	virtual void init();
	virtual void initWindows();
	virtual void update_centers();
	virtual void dissimilarity();
	virtual void update_membership();
	virtual void update_membership_bias();
	float Nt(int i, int j);
	void getAlpha(int i);
	void update_cluster();
	void exp();
	void exec(int c, float m, float e, int i_max, float n0, float seuil,float cst,bool fcm,int iter);
	void exec();
	void exportMembership(std::string file);
	
};
Mat convertCA(Mat a);
#endif