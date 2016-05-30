#pragma once
#pragma once
#ifndef _CA_
#define _CA_
#include "RgbFCM.h"
Mat convert_CA(Mat a);
struct CA : public RgbFCM
{
	bool* alive;
	Mat* Uexp;
	float* card;
	int iter_max;
	int c_final;
	float n0;
	float alpha;
	float seuil;
	float t_max;
	float LookTable_U[10000];
	float** LookTable_D;
	float* LookTable_SQRT;
	int iter;
	CA();
	~CA();
	virtual void init();
	virtual void update_centers();
	virtual void dissimilarity();
	virtual void update_membership();
	virtual void update_membership_bias();
	void exportMembership(std::string file);
	float Nt(int i, int j);
	void getAlpha(int i);
	void update_cluster();
	void exp();
	void exec(int c, float m, float e, int i_max, float n0, float seuil, float t_max, bool alt);
	void exec();
};

#endif