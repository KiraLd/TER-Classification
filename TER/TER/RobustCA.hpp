#pragma once
#ifndef _ROBUST_CA_
#define _ROBUST_CA_
#include "RgbFCM.h"
Mat convertRobust(Mat a);
struct RobustCA : public RgbFCM
{
	bool* alive;
	Mat* Uexp;
	Mat* robust_D;
	Mat* loss_D;
	float* T;
	float* S;
	float* robust_card;
	float* card;
	int iter_max;
	int c_final;
	float n0;
	float alpha;
	float seuil;
	float t_max;
	float delta;
	int iter;
	RobustCA();
	~RobustCA();
	virtual void init();
	virtual void update_centers();
	virtual void dissimilarity();
	void wp();
	void TS();
	virtual void update_membership();
	virtual void update_membership_bias();
	void exportMembership(std::string file);
	float Nt(int i, int j);
	void getAlpha(int i);
	void getAlphaAlt(int i);
	void update_cluster();
	void exp();
	void exec(int c, float m, float e, int i_max, float n0, float seuil, float t_max, bool alt,float delta);
	void exec();
};

#endif

