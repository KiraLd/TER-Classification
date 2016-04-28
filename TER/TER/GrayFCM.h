#pragma once
#ifndef _NDG_FCM_
#define _NDG_FCM_
#include "FuzzyClustering.hpp"

struct GrayFCM: public FuzzyClustering
{
	Mat* D;
	virtual void init();
	virtual void update_centers();
	virtual void dissimilarity();
	virtual void update_membership();
	GrayFCM();
	~GrayFCM();
	virtual void exec(int c, float m, float e, int i_max);
};

#endif