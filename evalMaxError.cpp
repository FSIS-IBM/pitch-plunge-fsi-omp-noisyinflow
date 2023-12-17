#include<cmath>
#include "evalMaxError.h"
float evalMaxError(int m_temp, int n_temp, float* Q, float* QA_old)
{
    float big, error;

    big=0.0;
	#pragma omp parallel for default(shared) private(error) reduction(max:big) schedule(dynamic)
    for(int i=0;i<n_temp;i++)
    {
        for(int j=0;j<m_temp;j++)
        {
            error=abs(*(Q+i*m_temp+j)-*(QA_old+i*m_temp+j));
            if (error>big)
                big = error;
        }
    }

    return big;
}
