#pragma once
#pragma once
#ifndef _RGB_FCM_
#define _RGB_FCM_
#include "FuzzyClustering.hpp"

struct RgbFCM : public FuzzyClustering
{
	Mat* D;
	void init();
	void update_centers();
	void dissimilarity();
	void update_membership();
	RgbFCM();
	~RgbFCM();
	void exec(int c, float m, float e, int i_max);
};

#endif
