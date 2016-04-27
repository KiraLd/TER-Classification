#pragma once
#ifndef _NDG_FCM_
#define _NDG_FCM_
#include "FuzzyClustering.hpp"

struct GrayFCM: public FuzzyClustering
{
	Mat* D;
	void init();
	void update_centers();
	void dissimilarity();
	void update_membership();
	GrayFCM();
	~GrayFCM();
	void exec(int c, float m, float e, int i_max);
};

#endif