#include "solverUVGaussSeidelSOR.h"
#include "evalMaxError.h"

// Solving momentum equations using Gauss Seidel Successive Over relaxation
void solverUVGaussSeidelSOR(int m, int m_u, int n_u, int m_v, int n_v, float* a_u, float* b_u, float* c_u, float* d_u, float* e_u, float* f_u, float* Bf_u, float* a_v, float* b_v, float* c_v, float* d_v, float* e_v, float* f_v, float* Bf_v, float* x1, float* y1, float* y, float* pressure, float mass_in, float domain_width, float* uA_old, float* vA_old, float* u, float* v, float* u_red, float* v_red, float* u_black, float* v_black, float* u_red_old, float* v_red_old, float* u_black_old, float* v_black_old, float* max_error1, float* max_error2)
{
    float mass_out;
        
		
	#pragma omp parallel for default(shared) schedule(dynamic)
	for(int i=0;i<n_u;i++)
	{
		int j1;
		
		j1 = (i%2==0) ? 1 : 2;
		for(int j=j1;j<m_u;j+=2)
		{
		   *(u_red+i*m_u+j)=*(u+i*m_u+j);
			*(u_black+i*m_u+j)=0.0;
		}

		j1 = (i%2==0) ? 2 : 1;
		for(int j=j1;j<m_u;j+=2)
		{
			*(u_red+i*m_u+j)=0.0;
			*(u_black+i*m_u+j)=*(u+i*m_u+j);
		}
	}

	#pragma omp parallel for default(shared) schedule(dynamic)
	for(int i=0;i<n_v;i++)
	{
		int j1;
		
		j1 = (i%2==0) ? 1 : 2;
		for(int j=j1;j<m_v;j+=2)
		{
			*(v_red+i*m_v+j)=*(v+i*m_v+j);
			*(v_black+i*m_v+j)=0.0;
		}

		j1 = (i%2==0) ? 2 : 1;
		for(int j=j1;j<m_v;j+=2)
		{
			*(v_red+i*m_v+j)=0.0;
			*(v_black+i*m_v+j)=*(v+i*m_v+j);
		}
	}
	
	
	// Solution of the x and y momentum equation using Gauss seidel SOR with Red-Black tagging
	
	for (int p=0; p<GSite; p++) // GSite - Gauss Seidel Iterations
	{
		#pragma omp parallel for collapse(2)
		for(int i=0;i<n_u;i++)
		{
			for(int j=0;j<m_u;j++)
			{
				*(uA_old+i*m_u+j)=*(u+i*m_u+j);
				*(u_red_old+i*m_u+j)=*(u_red+i*m_u+j);
				*(u_black_old+i*m_u+j)=*(u_black+i*m_u+j);
			}
		}

		// Calculate u-red and u-black velocities
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=0; i<n_u; i++)
		{
			int j1;
			j1 = (i%2==0) ? 1 : 2;
			
			for (int j=j1; j<m_u; j+=2)
			{
				if(i==0)
				{
					*(u_red+i*m_u+j) = *(u_black+(i+1)*m_u+j);
				}
				else if(i==n_u-1)
				{
				   *(u_red+i*m_u+j) = *(u_black+(i-1)*m_u+j);
				}
				else
				{
					if(j==m_u-1)
						*(u_red+i*m_u+j)=*(u_black+i*m_u+j-1);
					else
						*(u_red+i*m_u+j) = alpha_SOR*((*(b_u+i*m_u+j)*(*(u_black+(i+1)*m_u+j)) + *(c_u+i*m_u+j)*(*(u_black+(i-1)*m_u+j)) + *(d_u+i*m_u+j)*(*(u_black+i*m_u+j+1))+ *(e_u+i*m_u+j)*(*(u_black+i*m_u+j-1)) + *(f_u+i*m_u+j) + (dt/(*(x1+i*m+(j+1)) - *(x1+i*m+j)))*(*(pressure+i*m+j) - *(pressure+i*m+(j+1))) + dt*(*(Bf_u+i*m_u+j))) / (*(a_u+i*m_u+j))) + (1.0-alpha_SOR)*(*(u_red_old+i*m_u+j));
				}
			}
		}
		
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=0; i<n_u; i++)
		{
			int j1;
			j1 = (i%2==0) ? 2 : 1;
			
			for (int j=j1; j<m_u; j+=2)
			{
				if(i==0)
				{
					*(u_black+i*m_u+j) = *(u_red+(i+1)*m_u+j);
				}
				else if(i==n_u-1)
				{
					*(u_black+i*m_u+j) = *(u_red+(i-1)*m_u+j);
				}
				else
				{
					if(j==m_u-1)
						*(u_black+i*m_u+j)=*(u_red+i*m_u+j-1);
					else
						*(u_black+i*m_u+j) = alpha_SOR*((*(b_u+i*m_u+j)*(*(u_red+(i+1)*m_u+j)) + *(c_u+i*m_u+j)*(*(u_red+(i-1)*m_u+j)) + *(d_u+i*m_u+j)*(*(u_red+i*m_u+j+1)) + *(e_u+i*m_u+j)*(*(u_red+i*m_u+j-1)) + *(f_u+i*m_u+j) + (dt/(*(x1+i*m+(j+1)) - *(x1+i*m+j)))*(*(pressure+i*m+j) - *(pressure+i*m+(j+1))) + dt*(*(Bf_u+i*m_u+j))) / (*(a_u+i*m_u+j))) + (1.0-alpha_SOR)*( *(u_black_old+i*m_u+j));
				}
			}
		}
		
		// Generating the actual u using red and black u velocities
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=0;i<n_u;i++)
		{
			int j1;
			
			j1 = (i%2==0) ? 1 : 2;
			for(int j=j1;j<m_u;j+=2)
				*(u+i*m_u+j)=*(u_red+i*m_u+j);

			j1 = (i%2==0) ? 2 : 1;
			for(int j=j1;j<m_u;j+=2)
				*(u+i*m_u+j)=*(u_black+i*m_u+j);
		}
		

		// For every grid point and iteration calculation of error
		max_error1[p] = evalMaxError(m_u, n_u, (float*)u, (float*)uA_old); // For u velocity


		// Assignment of current values to the old values
		#pragma omp parallel for collapse(2)
		for(int i=0;i<n_v;i++)
		{
			for(int j=0;j<m_v;j++)
			{
				*(vA_old+i*m_v+j)=*(v+i*m_v+j);
				*(v_red_old+i*m_v+j)=*(v_red+i*m_v+j);
				*(v_black_old+i*m_v+j)=*(v_black+i*m_v+j);
			}
		}
		
		// Calculate v-red and v-black velocities
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=1; i<n_v-1; i++)
		{
			int j1;
			j1 = (i%2==0) ? 1 : 2;
			
			for (int j=j1; j<=m_v-1; j+=2)
			{
				if(j==m_v-1)
					*(v_red+i*m_v+j)=*(v_black+m_v*i+j-1);
				else
					*(v_red+i*m_v+j) = alpha_SOR*((*(b_v+i*m_v+j)*(*(v_black+m_v*(i+1)+j)) + *(c_v+i*m_v+j)*(*(v_black+m_v*(i-1)+j)) + *(d_v+i*m_v+j)*(*(v_black+m_v*i+j+1)) + *(e_v+i*m_v+j)*(*(v_black+m_v*i+j-1)) + *(f_v+i*m_v+j) + (dt/(*(y1+(i+1)*m+j) - *(y1+i*m+j)))*(*(pressure+i*m+j) - *(pressure+(i+1)*m+j)) + dt*(*(Bf_v+i*m_v+j))) / (*(a_v+i*m_v+j))) + (1.0-alpha_SOR)*(*(v_red_old+i*m_v+j));
			}
		}
		
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=1; i<n_v-1; i++)
		{
			int j1;
			j1 = (i%2==0) ? 2 : 1;
			
			for (int j=j1; j<=m_v-1; j+=2)
			{
				if(j==m_v-1)
					*(v_black+i*m_v+j)=*(v_red+m_v*i+j-1);
				else
					*(v_black+i*m_v+j) = alpha_SOR*((*(b_v+i*m_v+j)*(*(v_red+(i+1)*m_v+j)) + *(c_v+i*m_v+j)*(*(v_red+(i-1)*m_v+j)) + *(d_v+i*m_v+j)*(*(v_red+i*m_v+j+1)) + *(e_v+i*m_v+j)*(*(v_red+i*m_v+j-1)) + *(f_v+i*m_v+j) + (dt/(*(y1+(i+1)*m+j) - *(y1+i*m+j)))*(*(pressure+i*m+j) - *(pressure+(i+1)*m+j)) + dt*(*(Bf_v+i*m_v+j))) / (*(a_v+i*m_v+j)))+(1.0-alpha_SOR)*(*(v_black_old+i*m_v+j));
			}
		}

		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=0;i<n_v;i++)
		{
			int j1;
			
			j1 = (i%2==0) ? 1 : 2;
			for(int j=j1;j<m_v;j+=2)
				*(v+i*m_v+j)=(*(v_red+i*m_v+j));

			j1 = (i%2==0) ? 2 : 1;
			for(int j=j1;j<m_v;j+=2)
				*(v+i*m_v+j)=*(v_black+i*m_v+j);
		}
		
		
		// For every grid point and iteration calculation of error
		max_error2[p] = evalMaxError(m_v, n_v, (float*)v, (float*)vA_old);
		
		
		if (max_error1[p] <= tol && max_error2[p] <= tol) // Checking if the maximum error is less than tolerance
		{
			break;
		}
	}
	

	// Output mass flow correction 
	for(int ctr=1;ctr<=40;ctr++)
	{
		mass_out = 0.0;
		#pragma omp parallel for default(shared) reduction(+:mass_out) schedule(dynamic)
		for(int i=1;i<n_u-1;i++)
		{
			mass_out = mass_out + (*(u+i*m_u+(m_u-1)) * (*(y+i*(m-1)+(m_u-1)) - *(y+(i-1)*(m-1)+(m_u-1))));
		}
		
		#pragma omp parallel for default(shared) schedule(dynamic)
		for(int i=1; i<n_u-1; i++)
		{
			*(u+i*m_u+(m_u-1)) = *(u+i*m_u+(m_u-1)) + ((mass_in-mass_out)/domain_width);
		}
	}
	
	//cout << mass_in << "  " << (mass_in-mass_out) << "  " << mass_out << endl << endl;
}
