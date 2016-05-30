#pragma once
#ifndef _NDG_FCM_
#define _NDG_FCM_
#include "FuzzyClustering.hpp"

struct GrayFCM: public FuzzyClustering
{
	Mat* D;
	int x_s, y_s, x_e, y_e;
	virtual void init();
	virtual void update_centers();
	virtual void dissimilarity();
	virtual void update_membership();
	GrayFCM();
	~GrayFCM();
	virtual void exec(int c, float m, float e, int i_max, int x_start = 0, int y_start = 0, int x_end = 16, int y_end = 16);
	virtual void exec(int c, float m, float e, int i_max);
};

#endif