#pragma once
#ifndef _RGB_CA_
#define _RGB_CA_
#include "RgbFCM.h"
Mat convert(Mat a);
struct RgbCA : public RgbFCM
{
	bool* alive;
	Mat* Uexp;
	float* card;
	int iter_max;
	float n0;
	float alpha;
	float seuil;
	RgbCA();
	~RgbCA();
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
	void exec(int c, float m, float e, int i_max, float n0, float seuil);
};

#endif#pragma once
