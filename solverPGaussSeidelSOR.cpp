#include "solverPGaussSeidelSOR.h"
#include "evalMaxError.h"

void solverPGaussSeidelSOR(int m, int n, float* a_p, float* b_p, float* c_p, float* d_p, float* e_p, float* f_p, float* p_c, float* pA_old, float* p_c_red, float* p_c_black, float* p_c_red_old, float* p_c_black_old, float* max_error3)
{
	// Gauss Seidel iterations with red-black SOR for pressure correction equation
	
	for (int p=0; p<GSite2; p++)
	{
		#pragma omp parallel for collapse(2)
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<m;j++)
			{
				*(pA_old +i*m +j) = *(p_c+i*m+j);
				*(p_c_red_old +m*i+j) = *(p_c_red +m*i+j);
				*(p_c_black_old +m*i+j) = *(p_c_black +m*i+j);
			}
		}
		
		// For red and black pressure correction points
		
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=0; i<n; i++)
		{
			int j1;
			j1 = (i%2==0) ? 0 : 1;
			
			for (int j=j1; j<m; j+=2)
			{
				if (i==0)//bottom row
				{
					*(p_c_red+m*i+j) = *(p_c_black+m*(i+1)+j);
				}
				else if (i>0 && i<n-1 && j==0)//left column
				{
					*(p_c_red+m*i+j) = *(p_c_black+m*i+j+1);
				}
				else if (i>0 && i<n-1 && j==m-1)//right column
				{
					*(p_c_red+m*i+j) = *(p_c_black+m*i+j-1);
				}
				else if (i==n-1)//top row
				{
					*(p_c_red+m*i+j) = *(p_c_black+m*(i-1)+j);
				}
				else
				{
					*(p_c_red+m*i+j) = alpha_SOR*((*(b_p+i*m+j)*(*(p_c_black+m*(i+1)+j)) + *(c_p+i*m+j)*(*(p_c_black+m*(i-1)+j)) + *(d_p+i*m+j)*(*(p_c_black+m*i+j+1)) + *(e_p+i*m+j)*(*(p_c_black+m*i+j-1)) + *(f_p+i*m+j)) / (*(a_p+i*m+j))) + (1.0-alpha_SOR)*(*(p_c_red_old +m*i+j));
				}
			}
		}
		
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=0; i<n; i++)
		{
			int j1;
			j1 = (i%2==0) ? 1 : 0;
			
			for (int j=j1; j<m; j+=2)
			{
				if (i==0)//bottom row
				{
				   *(p_c_black+m*i+j) = *(p_c_red+m*(i+1)+j);
				}
				else if (i>0 && i<n-1 && j==0)//left column
				{
				   *(p_c_black+m*i+j) = *(p_c_red+m*i+j+1);
				}
				else if (i>0 && i<n-1 && j==m-1)//right column
				{
					*(p_c_black+m*i+j) = *(p_c_red+m*i+j-1);
				}
				else if (i==n-1)//top row
				{
					*(p_c_black+m*i+j) = *(p_c_red+m*(i-1)+j);
				}
				else
				{
					*(p_c_black+m*i+j) = alpha_SOR*((*(b_p+i*m+j)*(*(p_c_red+m*(i+1)+j))+ *(c_p+i*m+j)*(*(p_c_red+m*(i-1)+j)) + *(d_p+i*m+j)*(*(p_c_red+m*i+j+1)) + *(e_p+i*m+j)*(*(p_c_red+m*i+j-1)) + *(f_p+i*m+j)) / (*(a_p+i*m+j))) + (1.0-alpha_SOR)*(*(p_c_black_old +m*i+j));
				}
			}
		}
		
		// generating the actual pressure correction matrix
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=0;i<n;i++)
		{
			int j1;
			
			j1 = (i%2==0) ? 0 : 1;
			for (int j=j1; j<m; j+=2)
				*(p_c+i*m+j)=*(p_c_red+m*i+j);

			j1 = (i%2==0) ? 1 : 0;
			for (int j=j1; j<m; j+=2)
				*(p_c+i*m+j)=*(p_c_black +m*i+j);
		}


		// Calculate error
		*(max_error3+p) = evalMaxError(m, n, (float*)p_c, (float*)pA_old);
		
		
		// Checking if the error is below tolerance
		if (*(max_error3+p) <= tol)
			break;
	}
}