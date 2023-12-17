#include<cmath>
#include "assign_flag.h"
#include "nearestIBpoint.h"

// Function for tagging the flag to the different meshes

void assign_flag(int m_temp, int n_temp, int N_marker_temp, int c1_temp, int c2_temp, int c3_temp, int c4_temp, float* x_s_temp, float* y_s_temp, float* n1_temp, float* n2_temp, float* xtemp, float* ytemp, int* flagtemp)
{
    float temp;
    int k1;
    
	#pragma omp parallel for default(shared) private(k1,temp) schedule(dynamic)
    for(int i=0;i<n_temp;i++)
    {
        for (int j=0;j<m_temp;j++)
        {
            if(j>c1_temp && j<c2_temp && i>c3_temp && i<c4_temp)
            {
                k1 = nearestIBpoint(N_marker_temp, (float*)x_s_temp, (float*)y_s_temp, *(xtemp+i*m_temp+j), *(ytemp+i*m_temp+j));
                
                temp = ((*(xtemp+i*m_temp+j)) - (*(x_s_temp+k1))) * (*(n1_temp+k1)) + ((*(ytemp+i*m_temp+j)) - (*(y_s_temp+k1))) * (*(n2_temp+k1));

                if (temp>=0)
                    *(flagtemp+i*m_temp+j) = 1;
                else
                    *(flagtemp+i*m_temp+j) = 0;
            }
            else
                *(flagtemp+i*m_temp+j) = 1;
        }
    }
}
