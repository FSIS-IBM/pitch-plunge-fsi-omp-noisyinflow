#include<cmath>
#include "nearestIBpoint.h"

int nearestIBpoint(int N_marker_temp, float* x_s_temp, float* y_s_temp, float xtemp, float ytemp)
{
    float min_dis_temp, distance_temp;
    int k1_temp;

    min_dis_temp = 10.0;
    
    for(int ctr=0;ctr<N_marker_temp;ctr++)
    {
        distance_temp = sqrt((*(x_s_temp+ctr)-xtemp)*(*(x_s_temp+ctr)-xtemp) + (*(y_s_temp+ctr)-ytemp)*(*(y_s_temp+ctr)-ytemp));
        if (distance_temp < min_dis_temp)
        {
            min_dis_temp = distance_temp;
            k1_temp = ctr;
        }
    }

    return k1_temp;
}
