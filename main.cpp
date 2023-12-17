#include<iostream>
#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<omp.h>

#include "nearestIBpoint.h"
#include "assign_flag.h"
#include "evalMaxError.h"
#include "solverUVGaussSeidelSOR.h"
#include "solverPGaussSeidelSOR.h"

using namespace std;

const float pi = 3.141592654;
float Re = 300.0;
int GSite = 100, GSite2 = 500;
float alpha_SOR=1.50, tol=0.000001;
float dt=0.0002;
int start_timestep = 0, end_timestep = 750000, writeInterval=654;
float amp_x = 0.0, omega_x = 0.0, amp_y = 0.50, omega_y = 4.0; //body kinematics
float pitch, pitch_m = 0.0, pitch_a = 0.0, omega_pitch = 0.0, pitch_dot, phase = 0.0, pitch_axis = -0.50;
float AR = 0.12;
float lift3 = 0.0, drag3 = 0.0;
float damping = 0.0, freq = 16.0/3.0, cubic_stiff = 16.0, dens_ratio = 1.0;

float fung1(float var1, float var2, float var3) // var1 is time, var2 is position, var3 is velocity
{
     float fun1;
     fun1 = var3;
     return fun1;
}
float fung2(float var1, float var2, float var3) // var1 is time, var2 is position, var3 is velocity
{
     float fun2;
     fun2 = - damping*var3 - freq*freq*(-var2 + cubic_stiff*pow(var2,3)) - (pitch_axis/(((1.0+AR*AR)/16.0)+pitch_axis*pitch_axis))*cos(var2)*(-amp_y*omega_y*omega_y*cos(omega_y*var1)) + ((32.0*dens_ratio)/pi)*(pitch_axis/(AR*(1.0+AR*AR+16.0*pitch_axis*pitch_axis)))*(lift3*cos(var2) + drag3*sin(var2));
     return fun2;
}

int main()
{
    float x_lower_left = -7.50, x_lower_right = 30.00, y_lower_left = -12.50, y_upper_left = 12.50;
    float Hx2 = 2.00, Hx1 = -0.5, Hy2 = 1.75, Hy1 = -1.75;
    int m_div1/* = 15*/, m_div2 = 625, m_div3/* = 75*/, m, n_div1/* = 5*/, n_div2 = 875, n_div3/* = 15*/, n,m_u,n_u,m_v,n_v;
    float radius_a = 0.5, radius_b = AR*radius_a;
    int N_marker = 1500; //Please take a number divisible by 4
    float rate_xu = 1.05, rate_xd = 1.008, rate_yu = 1.025, rate_yd = 1.025;
    float q_x1,q_x2,q_y1,q_y2,d_x1,d_x2,d_y1,d_y2,dx1, dx2, dy1, dy2;
    //float probe_location = 2.50; /*with respect to origin along x-axis */
    float time1 = start_timestep*dt, z0 = 0.0, z1 = 0.0, rk1[2],rk2[2],rk3[2],rk4[2];
	int c1,c2,c3,c4,c5,probe_i,probe_j;

    dx2 = (Hx2-Hx1)/m_div2;
    m_div1 = ceil((log((((Hx1-x_lower_left)/dx2)*(rate_xu-1))+1))/(log(rate_xu)));
    m_div3 = ceil((log((((x_lower_right-Hx2)/dx2)*(rate_xd-1))+1))/(log(rate_xd)));
    d_x1 = (dx2*((pow(rate_xu,m_div1))-1))/(rate_xu-1);

    dy2 = (Hy2-Hy1)/n_div2;
    n_div1 = ceil((log((((Hy1-y_lower_left)/dy2)*(rate_yu-1))+1))/(log(rate_yu)));
    n_div3 = ceil((log((((y_upper_left-Hy2)/dy2)*(rate_yd-1))+1))/(log(rate_yd)));
    d_y1 = (dy2*((pow(rate_yu,n_div1))-1))/(rate_yu-1);

    m=m_div1+m_div2+m_div3+2;
    n=n_div1+n_div2+n_div3+2;
    m_u = m-1;
    n_u = n;
    m_v = m;
    n_v = n-1;

	// ##### Defining the dense zone #####

	/*c1 = m_div1 + (((xc-Hx1)-(radius_a+10*dx2))/dx2);
	c2 = m_div1+m_div2-(((Hx2-xc)-(radius_a+10*dx2))/dx2);
	c3 = n_div1 + (((yc-Hy1)-(radius_b+10*dy2))/dy2);
	c4 = n_div1+n_div2-(((Hy2-yc)-(radius_b+10*dy2))/dy2);*/
	c1 = m_div1;			// Left boundary
	c2 = m_div1+m_div2;		// Right boundary
	c3 = n_div1;			// Bottom
	c4 = n_div1+n_div2;		// Top

    /*c5 = ceil((log((((probe_location-Hx2)/dx2)*(rate_xd-1))+1))/(log(rate_xd)));
    probe_i = n_div1+ceil((n_div2/2))+1;
    probe_j = m_div1+m_div2+c5;*/

    float x[n-1][m-1], y[n-1][m-1], x1[n][m], y1[n][m], x_u[n_u][m_u], y_u[n_u][m_u], x_v[n_v][m_v], y_v[n_v][m_v];
    int flag[n-1][m-1], flag_tu[n_u][m_u], flag_u[n_u][m_u], flag_tv[n_v][m_v], flag_v[n_v][m_v], flag1[n][m];
    float x_s[N_marker], y_s[N_marker], n1[N_marker], n2[N_marker], t1[N_marker], t2[N_marker],v1,v2, theta, xc, yc, x_pivot, y_pivot, u_pivot, v_pivot;
    float domain_width, distance, min_dis, phi, temp;
	int k1, sig, n1_i, n1_j;
	const float dist = sqrt(dx2*dx2+dy2*dy2);


	//  ##########  Generation of the structured non-uniform Cartesian mesh for the flow solution   ##########

    // ##### Mesh Generation: Calculation of the coordinates of the grid points #####
    for(int i=0;i<n-1;i++)
    {
        dx1 = dx2*(pow(rate_xu,(m_div1-1)));
        for(int j=0;j<m-1;j++)
        {
            if (j<=m_div1)
            {
                if (j==0)
                    q_x1=Hx1-d_x1;
                else
                {
                    q_x1 = q_x1 + dx1;
                    dx1 = dx1/rate_xu;
                }
                x[i][j] = q_x1;
            }
            else if(j<=m_div1+m_div2)
            {
                x[i][j] = q_x1 + (j-m_div1)*dx2;
                q_x2 = x[i][j];
                d_x2 = dx2;
            }
            else if(j<=m_div1+m_div2+m_div3)
            {
                q_x2 = q_x2 + d_x2;
                x[i][j] = q_x2;
                d_x2 = d_x2 * rate_xd;
            }
        }
    }
    for(int j=0;j<m-1;j++)
    {
        dy1 = dy2*(pow(rate_yu,(n_div1-1)));
        for(int i=0;i<n-1;i++)
        {
            if (i<=n_div1)
            {
                if (i==0)
                    q_y1=Hy1-d_y1;
                else
                {
                    q_y1 = q_y1 + dy1;
                    dy1 = dy1/rate_yu;
                }
                y[i][j] = q_y1;
            }
            else if(i<=n_div1+n_div2)
            {
                y[i][j] = q_y1 + (i-n_div1)*dy2;
                q_y2 = y[i][j];
                d_y2 = dy2;
            }
            else if(i<=n_div1+n_div2+n_div3)
            {
                q_y2 = q_y2 + d_y2;
                y[i][j] = q_y2;
                d_y2 = d_y2 * rate_yd;
            }
        }
    }

    // ##### Calculation of the coordinate of the center of each grid cell #####
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            if (j==0)
            {
                if(i<n-1)
                    x1[i][j] = x[i][j];
                else
                    x1[i][j] = x[i-1][j];
            }
            else if (j==m-1)
            {
                if(i<n-1)
                    x1[i][j] = x[i][j-1];
                else
                    x1[i][j] = x[i-1][j-1];
            }
            else
            {
                if(i<n-1)
                    x1[i][j] = (x[i][j]+x[i][j-1])/2.0;
                else
                    x1[i][j] = (x[i-1][j]+x[i-1][j-1])/2.0;
            }

            if (i==0)
            {
                if(j<m-1)
                    y1[i][j] = y[i][j];
                else
                    y1[i][j] = y[i][j-1];
            }
            else if (i==n-1)
            {
                if(j<m-1)
                    y1[i][j] = y[i-1][j];
                else
                    y1[i][j] = y[i-1][j-1];
            }
            else
            {
                if(j<m-1)
                    y1[i][j] = (y[i][j]+y[i-1][j])/2.0;
                else
                    y1[i][j] = (y[i][j-1]+y[i-1][j-1])/2.0;
            }
        }
    }

    // ##### Calculation of the coordinate of the u and v velocity locations #####
    for (int i=0;i<n_u;i++)
    {
        for(int j=0;j<m_u;j++)
        {
            x_u[i][j] = x[0][j];
            y_u[i][j] = y1[i][0];
        }
    }
    for (int i=0;i<n_v;i++)
    {
        for(int j=0;j<m_v;j++)
        {
            x_v[i][j] = x1[0][j];
            y_v[i][j] = y[i][0];
        }
    }

    domain_width = 0.0;
    for(int i=1;i<n_u-1;i++)
    {
        domain_width = domain_width + (y[i][m_u-1]-y[i-1][m_u-1]);
    }


	// Prints the grid details and iteration details
	cout << "m_div1 = " << m_div1 << "  m_div2 = " << m_div2 << "  m_div3 = " << m_div3 << endl;
	cout << "n_div1 = " << n_div1 << "  n_div2 = " << n_div2 << "  n_div3 = " << n_div3 << endl;
	cout << "m = " << m << "  n = " << n << endl;
	cout << "m_u = " << m_u << "   n_u = " << n_u << "   m_v = " << m_v << "   n_v = " << n_v << endl;
	cout << "dx1 = " << dx1 << "   dx2 = " << dx2 << endl;
	cout << "dy1 = " << dy1 << "   dy2 = " << dy2 << endl;
	/*cout << probe_i << " " << probe_j << endl;
    cout << x[probe_i][probe_j] << " " << y[probe_i][probe_j] << endl;
    cout << x1[probe_i][probe_j] << " " << y1[probe_i][probe_j] << endl;
    cout << x_u[probe_i][probe_j] << " " << y_u[probe_i][probe_j] << endl;
    cout << x_v[probe_i][probe_j] << " " << y_v[probe_i][probe_j] << endl;*/

    theta = (2.0*pi)/N_marker;

    if(start_timestep == 0)
    {
        // ##### Solid body kinematics #####

		x_pivot = amp_x*sin(omega_x*0*dt);
        y_pivot = amp_y*cos(omega_y*0*dt);
        //pitch = pitch_m + pitch_a*sin(omega_pitch*0*dt + phase);
        pitch = z0;
        xc = -pitch_axis*cos(pitch) + x_pivot;
        yc = pitch_axis*sin(pitch) + y_pivot;


		// ##### Defining the Solid body #####

        for(int i=0;i<N_marker;i++)
        {
            x_s[i]=radius_a*cos(pi-theta*i)*cos((pitch)) + radius_b*sin(pi-theta*i)*sin((pitch)) + xc;
            y_s[i]=-radius_a*cos(pi-theta*i)*sin((pitch)) + radius_b*sin(pi-theta*i)*cos((pitch)) + yc;
        }
        for (int i = 0;i<N_marker;i++)
        {
            if(i==0)
            {
                v1 = x_s[i+1] - x_s[N_marker-1];
                v2 = y_s[i+1] - y_s[N_marker-1];
            }
            else if (i==N_marker-1)
            {
                v1 = x_s[0] - x_s[i-1];
                v2 = y_s[0] - y_s[i-1];
            }
            else
            {
                v1 = x_s[i+1] - x_s[i-1];
                v2 = y_s[i+1] - y_s[i-1];
            }
            n1[i] = -v2;
            n2[i] = v1;
            n1[i] = (-v2)/sqrt((-v2)*(-v2)+v1*v1);
            n2[i] = v1/sqrt((-v2)*(-v2)+v1*v1);
            //t1[i] = v1/sqrt(v1*v1+v2*v2);
            //t2[i] = v2/sqrt(v1*v1+v2*v2);
        }
        ofstream file("solidBoundary.dat");
        if (file.is_open())
        {
            for(int i=0; i<N_marker; i++)
            {
                file << x_s[i] << " " << y_s[i] << " " << n1[i] << " " << n2[i] << endl;
            }
            file << x_s[0] << " " << y_s[0] << " " << n1[0] << " " << n2[0];
            file.close();
        }


		//#############    Identifying the fluid points, solid points and the forcing points     #############

		//Identifying the fluid points, solid points and the forcing points for the Eulerian mesh
		assign_flag(m-1, n-1, N_marker, c1, c2, c3, c4, (float*)x_s, (float*)y_s, (float*)n1, (float*)n2, (float*)x, (float*)y, (int*)flag);

		//Identifying the fluid points, solid points and the forcing points for the cell centers
		assign_flag(m, n, N_marker, c1, c2, c3, c4, (float*)x_s, (float*)y_s, (float*)n1, (float*)n2, (float*)x1, (float*)y1, (int*)flag1);

		//Identifying the fluid points, solid points and the forcing points for the u-mesh grid points
		assign_flag(m_u, n_u, N_marker, c1, c2, c3, c4, (float*)x_s, (float*)y_s, (float*)n1, (float*)n2, (float*)x_u, (float*)y_u, (int*)flag_tu);

		#pragma omp parallel for collapse(2)
		for(int i=0;i<n_u;i++)
		{
			for (int j=0;j<m_u;j++)
			{
				flag_u[i][j] = (flag_tu[i][j]==0) ? 3 : 1;
			}
		}

		//Identifying the fluid points, solid points and the forcing points for the v-mesh grid points
		assign_flag(m_v, n_v, N_marker, c1, c2, c3, c4, (float*)x_s, (float*)y_s, (float*)n1, (float*)n2, (float*)x_v, (float*)y_v, (int*)flag_tv);

		#pragma omp parallel for collapse(2)
		for(int i=0;i<n_v;i++)
		{
			for (int j=0;j<m_v;j++)
			{
				flag_v[i][j] = (flag_tv[i][j]==0) ? 3 : 1;
			}
		}


        ofstream file1("grid.dat");
        if (file1.is_open())
        {
            file1 << "ZONE T=DATA I=" << m-1 << " " << "J=" << n-1 << endl;
            for(int i=0;i<n-1;i++)
            {
                for(int j=0;j<m-1;j++)
                {
                file1 << x[i][j] << " " << y[i][j] << " " << flag[i][j] << endl;
                }
            }
            file1.close();
        }
        ofstream file2("grid1.dat");
        if (file2.is_open())
        {
            file2 << "ZONE T=DATA I=" << m << " " << "J=" << n << endl;
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<m;j++)
                {
                    file2 << x1[i][j] << " " << y1[i][j] << " " << flag1[i][j] << endl;
                }
            }
            file2.close();
        }
        ofstream file3("grid_u.dat");
        if (file3.is_open())
        {
            file3 << "ZONE T=DATA I=" << m_u << " " << "J=" << n_u << endl;
            for(int i=0;i<n_u;i++)
            {
                for(int j=0;j<m_u;j++)
                {
                    file3 << x_u[i][j] << " " << y_u[i][j] << " " << flag_u[i][j] << endl;
                }
            }
            file3.close();
        }
        ofstream file4("grid_v.dat");
        if (file4.is_open())
        {
            file4 << "ZONE T=DATA I=" << m_v << " " << "J=" << n_v << endl;
            for(int i=0;i<n_v;i++)
            {
                for(int j=0;j<m_v;j++)
                {
                    file4 << x_v[i][j] << " " << y_v[i][j] << " " << flag_v[i][j] << endl;
                }
            }
            file4.close();
        }
    }


	float u[n_u][m_u], v[n_v][m_v], u_1[n_u][m_u], v_1[n_v][m_v], u_old[n_u][m_u], v_old[n_v][m_v], pressure[n][m], uf1[n][m], vf1[n][m], u_inlet[750001];
	float ue,uw,un,us,ve,vw,vn,vs,ue_old,uw_old,un_old,us_old,ve_old,vw_old,vn_old,vs_old;
	float a_u[n_u][m_u], b_u[n_u][m_u], c_u[n_u][m_u], d_u[n_u][m_u], e_u[n_u][m_u], f_u[n_u][m_u], Bf_u[n_u][m_u];
	float a_v[n_v][m_v], b_v[n_v][m_v], c_v[n_v][m_v], d_v[n_v][m_v], e_v[n_v][m_v], f_v[n_v][m_v], Bf_v[n_v][m_v];
	float p_c[n][m], a_p[n][m], b_p[n][m], c_p[n][m], d_p[n][m], e_p[n][m], f_p[n][m];
	float big,residual_flux,error,A1,B1,C1,mass_in,mass_out;
	float x_intersect, y_intersect, x_mirror, y_mirror, x_mirror1, y_mirror1, slope, dist_a, dist_b;
	float U_intp, V_intp, u_mirror, v_mirror, u_body, v_body;
    float uA_old[n_u][m_u], vA_old[n_v][m_v],pA_old[n][m],max_error1[GSite],max_error2[GSite],max_error_Bf_u,max_error_Bf_v,max_error3[GSite2];
    //float surface_pr[N_marker],lift3,drag3;
    float uf, uf_n, uf_s, uf_old, vf, vf_e, vf_w, vf_old, temp1, temp_D, temp_L;
    float u_red[n_u][m_u],u_black[n_u][m_u],v_red[n_v][m_v],v_black[n_v][m_v],p_c_red[n][m],p_c_black[n][m];
    float u_red_old[n_u][m_u],u_black_old[n_u][m_u],v_red_old[n_v][m_v],v_black_old[n_v][m_v],p_c_red_old[n][m],p_c_black_old[n][m];
    string fileName,fileName1,fileName4,fileName7,fileName8,fileName9,fileName10;
    //float x_1,y_1,x_2,y_2,x_3,y_3,u_2,u_3,v_2,v_3,p_2,p_3,b2_u_x,b3_u_y,b2_v_x,b3_v_y,b1,b2,b3,Det,dpdn;


    ifstream fileread("u_noisy_delw5_q0pt0655_2.txt");
    if (fileread.is_open())
    {
        for(int i=0; i<750001; i++)
        {
            fileread >> u_inlet[i];
        }
        fileread.close();
    }


	//  ##########  Initialisation of the fields including boundary conditions   ##########

	// When starting from a different time step
    if (start_timestep != 0)
    {
        stringstream ss;
        ss << start_timestep;
        fileName = ss.str();

        fileName1 = "Cavity_u"+fileName+".txt";
        ifstream fileread1(fileName1.c_str());
        if (fileread1.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    fileread1 >> u[i][j];
                }
            }
            fileread1.close();
        }
        else
        {
            cout << "Files not present. \nPlease restart from time step = 0 or the time step for which files exist" << endl;
            return 1;
        }

        fileName1 = "Cavity_u_old"+fileName+".txt";
        ifstream fileread1a(fileName1.c_str());
        if (fileread1a.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    fileread1a >> u_old[i][j];
                }
            }
            fileread1a.close();
        }

        fileName1 = "Cavity_v"+fileName+".txt";
        ifstream fileread2(fileName1.c_str());
        if (fileread2.is_open())
        {
            for(int i=0; i<n_v; i++)
            {
                for(int j=0; j<m_v; j++)
                {
                    fileread2 >> v[i][j];
                }
            }
            fileread2.close();
        }

        fileName1 = "Cavity_v_old"+fileName+".txt";
        ifstream fileread2a(fileName1.c_str());
        if (fileread2a.is_open())
        {
            for(int i=0; i<n_v; i++)
            {
                for(int j=0; j<m_v; j++)
                {
                    fileread2a >> v_old[i][j];
                }
            }
            fileread2a.close();
        }

        fileName1 = "pressure"+fileName+".txt";
        ifstream fileread3(fileName1.c_str());
        if (fileread3.is_open())
        {
            for(int i=0; i<n; i++)
            {
                for(int j=0; j<m; j++)
                {
                    fileread3 >> pressure[i][j];
                }
            }
            fileread3.close();
        }

        fileName1 = "p_c"+fileName+".txt";
        ifstream fileread3a(fileName1.c_str());
        if (fileread3a.is_open())
        {
            for(int i=0; i<n; i++)
            {
                for(int j=0; j<m; j++)
                {
                    fileread3a >> p_c[i][j];
                }
            }
            fileread3a.close();
        }

		for(int i=0;i<n;i++)
        {
            int j1;

			j1 = (i%2==0) ? 0 : 1;
			for (int j=j1; j<m; j+=2)
            {
                p_c_red[i][j]=p_c[i][j];
                p_c_black[i][j]=0.0;
            }

			j1 = (i%2==0) ? 1 : 0;
            for (int j=j1; j<m; j+=2)
            {
                p_c_red[i][j]=0.0;
                p_c_black[i][j]=p_c[i][j];
            }
        }
    }
    else // When starting from time1 = 0
    {
        //declaring the initial condition
        for(int i=0;i<n_u;i++)
        {
            for(int j=0;j<m_u;j++)
				u[i][j] = (flag_u[i][j]==1) ? 1.0 : 0.0;
        }
        for(int i=0;i<n_u;i++)
        {
            for(int j=0;j<m_u;j++)
            {
                u_red[i][j]=u[i][j];
                u_black[i][j]=u[i][j];
            }
        }
        for(int i=0;i<n_v;i++)
        {
            for(int j=0;j<m_v;j++)
            {
                v[i][j]=0.0f;
                v_red[i][j]=0.0f;
                v_black[i][j]=0.0f;
            }
        }
        for (int i=0; i<n; i++)
        {
            for(int j=0; j<m; j++)
            {
                pressure[i][j] = 0.0;
                p_c[i][j]= 0.0;
                p_c_red[i][j]=0.0;
                p_c_black[i][j]=0.0;
            }
        }
        for(int i=0;i<n_u;i++)
        {
            for(int j=0;j<m_u;j++)
            {
                Bf_u[i][j]=0.0;
            }
        }
        for(int i=0;i<n_v;i++)
        {
            for(int j=0;j<m_v;j++)
            {
                Bf_v[i][j]=0.0f;
            }
        }
        /*
        // Guess pressure field
        ifstream fileread("pressure1.txt");
		if (fileread.is_open())
		{
			for(i=0; i<n; i++)
			{
				for(j=0; j<m; j++)
				{
					fileread >> pressure[i][j];
				}
			}
			fileread.close();
		}*/

        // boundary conditions for u
        for(int i=0; i<n_u; i++)
        {
            for (int j=0; j<m_u; j++)
            {
                /*if ( i==0 ) // bottom
                {
                    u[i][j]= 0.0;
                }
                else if ( i== n_u-1 ) // top
                {
                    u[i][j]= 0.0;
                }*/
                if ( j==0 && i!=0 && i != n_u-1) //left
                {
                    u[i][j]= u_inlet[0];
                }
                /*else if ( j==m_u-1 && i!=0 && i != n_u-1 ) //right
                {
                    u[i][j]= 0.5;
                }*/
            }
        }

        //Boundary condition for v
        for(int i=0; i<n_v; i++)
        {
            for (int j=0; j<m_v; j++)
            {
                if ( i==0 ) // bottom
                {
                    v[i][j]= 0.0f;
                }
                else if ( i== n_v-1 ) // top
                {
                    v[i][j]= 0.0f;
                }
                if ( j==0 && i!=0 && i != n_v-1) //left
                {
                    v[i][j]= 0.0f;
                }
                /*else if ( j==m_v-1 && i!=0 && i != n_v-1 ) //right
                {
                    v[i][j]= 0.0;
                }*/
            }
        }

        //storing the values of u and v
        for (int i=0; i<n_u; i++)
        {
            for(int j=0; j<m_u; j++)
            {
                u_old[i][j]= u[i][j];
            }
        }
        for (int i=0; i<n_v; i++)
        {
            for(int j=0; j<m_v; j++)
            {
                v_old[i][j]= v[i][j];
            }
        }

        // Estimation of u from x-momentum
        for(int i=0; i<n_u; i++)
        {
            for (int j=1; j<m_u; j++)
            {
                if (j==m_u-1 && i!=0 && i!=n_u-1) //right boundary
                {
                    u[i][j] = u_old[i][j-1];
                }
                if (i==0)  //bottom boundary
                {
                    u[i][j] = u_old[i+1][j];
                }
                else if (i==1 && j!=m_u-1)//bottom row
                {
                    ue = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u_old[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u_old[i][j+1];
                    uw = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j];
                    un = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                    us = u_old[i-1][j];
                    vn = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                    vs = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j+1];
                    u[i][j] = u_old[i][j] - (dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw) - (dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us) + (dt/((x1[i][j+1]-x1[i][j])*Re))*(((u_old[i][j+1]-u_old[i][j])/(x[i][j+1]-x[i][j]))-((u_old[i][j]-u_old[i][j-1])/(x[i][j]-x[i][j-1]))) + (dt/((y[i][j]-y[i-1][j])*Re))*(((u_old[i+1][j]-u_old[i][j])/(y1[i+1][j]-y1[i][j]))-((u_old[i][j]-u_old[i-1][j])/(y1[i][j]-y1[i-1][j]))) + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]);
                }
                else if (i==n_u-2 && j!=m_u-1)//top row
                {
                    ue = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u_old[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u_old[i][j+1];
                    uw = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j];
                    un = u_old[i+1][j];
                    us = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i][j]+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i-1][j];
                    vn = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                    vs = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j+1];
                    u[i][j] = u_old[i][j] - (dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw) - (dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us) + (dt/((x1[i][j+1]-x1[i][j])*Re))*(((u_old[i][j+1]-u_old[i][j])/(x[i][j+1]-x[i][j]))-((u_old[i][j]-u_old[i][j-1])/(x[i][j]-x[i][j-1]))) + (dt/((y[i][j]-y[i-1][j])*Re))*(((u_old[i+1][j]-u_old[i][j])/(y1[i+1][j]-y1[i][j]))-((u_old[i][j]-u_old[i-1][j])/(y1[i][j]-y1[i-1][j]))) + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]);
                }
                else if (i==n_u-1) //top boundary
                {
                    u[i][j] = u_old[i-1][j];
                }
                else if (j!=m_u-1)
                {
                    ue = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u_old[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u_old[i][j+1];
                    uw = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j];
                    un = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                    us = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i][j]+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i-1][j];
                    vn = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                    vs = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j+1];
                    u[i][j] = u_old[i][j] - (dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw) - (dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us) + (dt/((x1[i][j+1]-x1[i][j])*Re))*(((u_old[i][j+1]-u_old[i][j])/(x[i][j+1]-x[i][j]))-((u_old[i][j]-u_old[i][j-1])/(x[i][j]-x[i][j-1]))) + (dt/((y[i][j]-y[i-1][j])*Re))*(((u_old[i+1][j]-u_old[i][j])/(y1[i+1][j]-y1[i][j]))-((u_old[i][j]-u_old[i-1][j])/(y1[i][j]-y1[i-1][j]))) + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]);
                }
            }
        }

        //Estimation of v from y-momentum
        for(int i=1; i<n_v-1; i++)
        {
            for (int j=1; j<=m_v-1; j++)
            {
                if (j==1)//left column
                {
                    ue = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                    uw = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i][j-1];
                    ve = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                    vw = v_old[i][j-1];
                    vn = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v_old[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v_old[i][j];
                    vs = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
                    v[i][j] = v_old[i][j] - (dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw) - (dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs) + (dt/((x[i][j]-x[i][j-1])*Re))*(((v_old[i][j+1]-v_old[i][j])/(x1[i][j+1]-x1[i][j]))-((v_old[i][j]-v_old[i][j-1])/(x1[i][j]-x1[i][j-1]))) + (dt/((y1[i+1][j]-y1[i][j])*Re))*(((v_old[i+1][j]-v_old[i][j])/(y[i+1][j]-y[i][j]))-((v_old[i][j]-v_old[i-1][j])/(y[i][j]-y[i-1][j]))) + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]);
                }
                else if (j==m_v-2)//right column
                {
                    ue = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                    uw = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i][j-1];
                    ve = v_old[i][j+1];
                    vw = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j-1]+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j];
                    vn = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v_old[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v_old[i][j];
                    vs = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
                    v[i][j] = v_old[i][j] - (dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw) - (dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs) + (dt/((x[i][j]-x[i][j-1])*Re))*(((v_old[i][j+1]-v_old[i][j])/(x1[i][j+1]-x1[i][j]))-((v_old[i][j]-v_old[i][j-1])/(x1[i][j]-x1[i][j-1]))) + (dt/((y1[i+1][j]-y1[i][j])*Re))*(((v_old[i+1][j]-v_old[i][j])/(y[i+1][j]-y[i][j]))-((v_old[i][j]-v_old[i-1][j])/(y[i][j]-y[i-1][j]))) + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]);
                }
                else if (j==m_v-1) //right boundary
                {
                    v[i][j] = v_old[i][j-1];
                }
                else
                {
                    ue = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                    uw = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i][j-1];
                    ve = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                    vw = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j-1]+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j];
                    vn = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v_old[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v_old[i][j];
                    vs = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
                    v[i][j] = v_old[i][j] - (dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw) - (dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs) + (dt/((x[i][j]-x[i][j-1])*Re))*(((v_old[i][j+1]-v_old[i][j])/(x1[i][j+1]-x1[i][j]))-((v_old[i][j]-v_old[i][j-1])/(x1[i][j]-x1[i][j-1]))) + (dt/((y1[i+1][j]-y1[i][j])*Re))*(((v_old[i+1][j]-v_old[i][j])/(y[i+1][j]-y[i][j]))-((v_old[i][j]-v_old[i-1][j])/(y[i][j]-y[i-1][j]))) + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]);
                }
            }
        }
    }
            /*ofstream filew1("Cavity_u.txt");
            if (filew1.is_open())
            {
                for(int i=0; i<n_u; i++)
                {
                    for(int j=0; j<m_u; j++)
                    {
                        filew1 << u[i][j] << " ";
                    }
                    filew1 << endl;
                }
                filew1.close();
            }

            ofstream filew2("Cavity_v.txt");
            if (filew2.is_open())
            {
                for(int i=0; i<n_v; i++)
                {
                    for(int j=0; j<m_v; j++)
                    {
                        filew2 << v[i][j] << " ";
                    }
                    filew2 << endl;
                }
                filew2.close();
            }

            ofstream filew3("pressure.txt");
            if (filew3.is_open())
            {
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<m; j++)
                    {
                        filew3 << pressure[i][j] << endl;
                    }
                }
                filew3.close();
            }

            ofstream filew4("Cavity_uold.txt");
            if (filew4.is_open())
            {
                for(int i=0; i<n_u; i++)
                {
                    for(int j=0; j<m_u; j++)
                    {
                        filew4 << u_old[i][j] << " ";
                    }
                    filew4 << endl;
                }
                filew4.close();
            }

            ofstream filew5("Cavity_vold.txt");
            if (filew5.is_open())
            {
                for(int i=0; i<n_v; i++)
                {
                    for(int j=0; j<m_v; j++)
                    {
                        filew5 << v_old[i][j] << " ";
                    }
                    filew5 << endl;
                }
                filew5.close();
            }*/

    // Calculation of the coefficients of x-momentum equation
    for(int i=1; i<n_u-1; i++)
    {
        for (int j=1; j<m_u-1; j++)
        {
            if (i==1)//bottom row
            {
                b_u[i][j] = (dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(1/(y1[i+1][j]-y1[i][j]));
                c_u[i][j] = 0.0;
                d_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1/(x[i][j+1]-x[i][j]));
                e_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1/(x[i][j]-x[i][j-1]));
            }
            else if (i==n_u-2)//top row
            {
                b_u[i][j] = 0.0;
                c_u[i][j] = (dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(1/(y1[i][j]-y1[i-1][j]));
                d_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1/(x[i][j+1]-x[i][j]));
                e_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1/(x[i][j]-x[i][j-1]));
            }
            else
            {
                b_u[i][j] = (dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(1/(y1[i+1][j]-y1[i][j]));
                c_u[i][j] = (dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(1/(y1[i][j]-y1[i-1][j]));
                d_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1/(x[i][j+1]-x[i][j]));
                e_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1/(x[i][j]-x[i][j-1]));
            }
            a_u[i][j] = 1+b_u[i][j]+c_u[i][j]+d_u[i][j]+e_u[i][j];
        }
    }

    // Calculation of the co-efficients of the y-momentum equation
    for(int i=1; i<n_v-1; i++)
    {
        for (int j=1; j<m_v-1; j++)
        {
            b_v[i][j] = (dt/(2.0*Re*(y1[i+1][j]-y1[i][j])))*(1/(y[i+1][j]-y[i][j]));
            c_v[i][j] = (dt/(2.0*Re*(y1[i+1][j]-y1[i][j])))*(1/(y[i][j]-y[i-1][j]));
            d_v[i][j] = (dt/(2.0*Re*(x[i][j]-x[i][j-1])))*(1/(x1[i][j+1]-x1[i][j]));
            e_v[i][j] = (dt/(2.0*Re*(x[i][j]-x[i][j-1])))*(1/(x1[i][j]-x1[i][j-1]));
            a_v[i][j] = 1.0+b_v[i][j]+c_v[i][j]+d_v[i][j]+e_v[i][j];
        }
    }

    //Calculation of the coefficients of the pressure correction equation
    for(int i=1; i<n-1; i++)
    {
        for (int j=1; j<m-1; j++)
        {
            if (i==1 && j==1)//bottom left corner
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 0.0;
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 0.0;
            }
            else if (i==1 && j>1 && j<m-2)//bottom row
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 0.0;
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else if (i==1 && j==m-2)//bottom right corner
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 0.0;
                d_p[i][j] = 0.0;
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else if (i>1 && i<n-2 && j==1)//left column
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 0.0;
            }
            else if (i>1 && i<n-2 && j==m-2)//right column
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 0.0;
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else if (i==n-2 && j==1)//top left corner
            {
                b_p[i][j] = 0.0;
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 0.0;
            }
            else if (i==n-2 && j>1 && j<m-2)//top row
            {
                b_p[i][j] = 0.0;
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else if (i ==n-2 && j==m-2)//top right corner
            {
                b_p[i][j] = 0.0;
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 0.0;
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            a_p[i][j] = b_p[i][j] + c_p[i][j] + d_p[i][j] + e_p[i][j];
        }
    }


    ofstream file5c("forcec.txt");
    //ofstream file6("residual_flux.txt");
    //ofstream file6a("velatprobe.txt");


	//   ###############     Time Marching Begins      ##############//

    for(int k=start_timestep+1;k<=end_timestep;k++)   //time marching
    {
        // ##### Solid body motion #####

		rk1[0] = fung1(time1,z0,z1);
        rk1[1] = fung2(time1,z0,z1);
        rk2[0] = fung1(time1+(dt/2.0),z0+dt*(rk1[0]/2.0),z1+dt*(rk1[1]/2.0));
        rk2[1] = fung2(time1+(dt/2.0),z0+dt*(rk1[0]/2.0),z1+dt*(rk1[1]/2.0));
        rk3[0] = fung1(time1+(dt/2.0),z0+dt*(rk2[0]/2.0),z1+dt*(rk2[1]/2.0));
        rk3[1] = fung2(time1+(dt/2.0),z0+dt*(rk2[0]/2.0),z1+dt*(rk2[1]/2.0));
        rk4[0] = fung1(time1+dt,z0+dt*rk3[0],z1+dt*rk3[1]);
        rk4[1] = fung2(time1+dt,z0+dt*rk3[0],z1+dt*rk3[1]);
        z0 = z0 + (dt/6.0)*(rk1[0]+(2.0*rk2[0])+(2.0*rk3[0])+rk4[0]);
        z1 = z1 + (dt/6.0)*(rk1[1]+(2.0*rk2[1])+(2.0*rk3[1])+rk4[1]);

        time1 = k*dt;		// Time instance calculation

        x_pivot = amp_x*sin(omega_x*k*dt);
        y_pivot = amp_y*cos(omega_y*k*dt);
        //pitch = pitch_m + pitch_a*sin(omega_pitch*k*dt + phase);
        pitch = z0;
        xc = -pitch_axis*cos(pitch) + x_pivot;
        yc = pitch_axis*sin(pitch) + y_pivot;
        u_pivot = amp_x*omega_x*cos(omega_x*dt*k);
        v_pivot = -amp_y*omega_y*sin(omega_y*dt*k);
        //pitch_dot = pitch_a*omega_pitch*cos(omega_pitch*k*dt + phase);
        pitch_dot = z1;


		// ##### Defining the Solid body #####

		for(int i=0;i<N_marker;i++)
        {
            x_s[i]=radius_a*cos(pi-theta*i)*cos((pitch)) + radius_b*sin(pi-theta*i)*sin((pitch)) + xc;
            y_s[i]=-radius_a*cos(pi-theta*i)*sin((pitch)) + radius_b*sin(pi-theta*i)*cos((pitch)) + yc;
        }
        for (int i = 0;i<N_marker;i++)
        {
            if(i==0)
            {
                v1 = x_s[i+1] - x_s[N_marker-1];
                v2 = y_s[i+1] - y_s[N_marker-1];
            }
            else if (i==N_marker-1)
            {
                v1 = x_s[0] - x_s[i-1];
                v2 = y_s[0] - y_s[i-1];
            }
            else
            {
                v1 = x_s[i+1] - x_s[i-1];
                v2 = y_s[i+1] - y_s[i-1];
            }
            n1[i] = -v2;
            n2[i] = v1;
            n1[i] = (-v2)/sqrt((-v2)*(-v2)+v1*v1);
            n2[i] = v1/sqrt((-v2)*(-v2)+v1*v1);
            //t1[i] = v1/sqrt(v1*v1+v2*v2);
            //t2[i] = v2/sqrt(v1*v1+v2*v2);
        }

		// Writing the solid boundary after every user specified write interval
		if ( (k%writeInterval) == 0)
        {
            stringstream ss;
            ss << k;
            fileName = ss.str();
            fileName = "solidBoundary"+fileName+".dat";
            ofstream file7(fileName.c_str());
            if (file7.is_open())
            {
                for(int i=0; i<N_marker; i++)
                {
                    file7 << x_s[i] << " " << y_s[i] << endl;
                }
                file7 << x_s[0] << " " << y_s[0];
                file7.close();
            }
        }


		//#############    Identifying the fluid points, solid points and the forcing points     #############

		//Identifying the fluid points, solid points and the forcing points for the cell centers
		assign_flag(m, n, N_marker, c1, c2, c3, c4, (float*)x_s, (float*)y_s, (float*)n1, (float*)n2, (float*)x1, (float*)y1, (int*)flag1);

		//Identifying the fluid points, solid points and the forcing points for the u-mesh grid points
		assign_flag(m_u, n_u, N_marker, c1, c2, c3, c4, (float*)x_s, (float*)y_s, (float*)n1, (float*)n2, (float*)x_u, (float*)y_u, (int*)flag_tu);

		#pragma omp parallel for collapse(2)
		for(int i=0;i<n_u;i++)
		{
			for (int j=0;j<m_u;j++)
			{
				flag_u[i][j] = (flag_tu[i][j]==0) ? 3 : 1;
			}
		}

		//Identifying the fluid points, solid points and the forcing points for the v-mesh grid points
		assign_flag(m_v, n_v, N_marker, c1, c2, c3, c4, (float*)x_s, (float*)y_s, (float*)n1, (float*)n2, (float*)x_v, (float*)y_v, (int*)flag_tv);

		#pragma omp parallel for collapse(2)
		for(int i=0;i<n_v;i++)
		{
			for (int j=0;j<m_v;j++)
			{
				flag_v[i][j] = (flag_tv[i][j]==0) ? 3 : 1;
			}
		}

		// Alpha_i flags to the data points for the q - forcing terms (tu and tv for u and v velocities)
		#pragma omp parallel for collapse(2)
		for(int i=0;i<n_u;i++)
		{
			for(int j=0;j<m_u;j++)
			{
				flag_tu[i][j] = (flag_u[i][j]==3) ? 1 : 0;
			}
		}
		#pragma omp parallel for collapse(2)
		for(int i=0;i<n_v;i++)
		{
			for(int j=0;j<m_v;j++)
			{
				flag_tv[i][j] = (flag_v[i][j]==3) ? 1 : 0;
			}
		}


		//  ##########   Estimating the intermediate u velocity i.e. generating u_tilde and v_tilde velocity field considering Bf_u = Bf_v = 0

		//  Calculates intermediate u velocity (u_tilde)
		#pragma omp parallel for default(shared) private(ue_old,uw_old,un_old,us_old,vn_old,vs_old,ue,uw,un,us,vn,vs) schedule(dynamic)
		for(int i=c3; i<c4; i++)
		{
			for (int j=c1; j<c2; j++)
			{
				ue_old = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u_old[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u_old[i][j+1];
				uw_old = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j];
				un_old = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
				us_old = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i][j]+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i-1][j];
				vn_old = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
				vs_old = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j+1];
				ue = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u[i][j+1];
				uw = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u[i][j];
				un = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u[i][j];
				us = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*u[i][j]+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*u[i-1][j];
				vn = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v[i][j+1];
				vs = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v[i-1][j+1];
				u_1[i][j] = u[i][j] - 1.5*((dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw)+(dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us)) + 0.5*((dt/(x1[i][j+1]-x1[i][j]))*(ue_old*ue_old-uw_old*uw_old)+(dt/(y[i][j]-y[i-1][j]))*(vn_old*un_old-vs_old*us_old)) + (dt/(Re*(x1[i][j+1]-x1[i][j])))*(((u[i][j+1]-u[i][j])/(x[i][j+1]-x[i][j]))-((u[i][j]-u[i][j-1])/(x[i][j]-x[i][j-1])))+(dt/(Re*(y[i][j]-y[i-1][j])))*(((u[i+1][j]-u[i][j])/(y1[i+1][j]-y1[i][j]))-((u[i][j]-u[i-1][j])/(y1[i][j]-y1[i-1][j]))) + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]);
			}
		}

		// Calculates intermediate v velocity (v_tilde)
		#pragma omp parallel for default(shared) private(ue_old,uw_old,ve_old,vw_old,vn_old,vs_old,ue,uw,ve,vw,vn,vs) schedule(dynamic)
		for(int i=c3; i<c4; i++)
		{
			for (int j=c1; j<c2; j++)
			{
				ue_old = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
				uw_old = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i][j-1];
				ve_old = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
				vw_old = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j-1]+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j];
				vn_old = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v_old[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v_old[i][j];
				vs_old = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
				ue = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u[i][j];
				uw = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u[i][j-1];
				ve = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v[i][j+1];
				vw = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*v[i][j-1]+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*v[i][j];
				vn = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v[i][j];
				vs = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v[i-1][j];
				v_1[i][j] = v[i][j] - 1.5*((dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw)+(dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs)) + 0.5*((dt/(x[i][j]-x[i][j-1]))*(ue_old*ve_old-uw_old*vw_old)+(dt/(y1[i+1][j]-y1[i][j]))*(vn_old*vn_old-vs_old*vs_old)) + (dt/(Re*(x[i][j]-x[i][j-1])))*(((v[i][j+1]-v[i][j])/(x1[i][j+1]-x1[i][j]))-((v[i][j]-v[i][j-1])/(x1[i][j]-x1[i][j-1])))+(dt/(Re*(y1[i+1][j]-y1[i][j])))*(((v[i+1][j]-v[i][j])/(y[i+1][j]-y[i][j]))-((v[i][j]-v[i-1][j])/(y[i][j]-y[i-1][j]))) + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]);
			}
		}


        /*if ( (k%writeInterval) == 0)
        {
            stringstream ss;
            ss << k;
            fileName = ss.str();
            fileName = "grid1_"+fileName+".dat";
            ofstream file8(fileName.c_str());
            if (file8.is_open())
            {
                file8 << "ZONE T=DATA I=" << m << " " << "J=" << n << endl;
                for(int i=0;i<n;i++)
                {
                    for(int j=0;j<m;j++)
                    {
                        file8 << x1[i][j] << " " << y1[i][j] << " " << flag1[i][j] << endl;
                    }
                }
                file8.close();
            }
            fileName = ss.str();
            fileName = "grid_"+fileName+".dat";
            ofstream file8a(fileName.c_str());
            if (file8a.is_open())
            {
                file8a << "ZONE T=DATA I=" << m-1 << " " << "J=" << n-1 << endl;
                for(int i=0;i<n-1;i++)
                {
                    for(int j=0;j<m-1;j++)
                    {
                        file8a << x[i][j] << " " << y[i][j] << " " << flag[i][j] << endl;
                    }
                }
                file8a.close();
            }
            fileName = ss.str();
            fileName = "grid_u_"+fileName+".dat";
            ofstream file9(fileName.c_str());
            if (file9.is_open())
            {
                file9 << "ZONE T=DATA I=" << m_u << " " << "J=" << n_u << endl;
                for(int i=0;i<n_u;i++)
                {
                    for(int j=0;j<m_u;j++)
                    {
                        file9 << x_u[i][j] << " " << y_u[i][j] << " " << flag_u[i][j] << endl;
                    }
                }
                file9.close();
            }
            fileName = ss.str();
            fileName = "grid_v_"+fileName+".dat";
            ofstream file10(fileName.c_str());
            if (file10.is_open())
            {
                file10 << "ZONE T=DATA I=" << m_v << " " << "J=" << n_v << endl;
                for(int i=0;i<n_v;i++)
                {
                    for(int j=0;j<m_v;j++)
                    {
                        file10 << x_v[i][j] << " " << y_v[i][j] << " " << flag_v[i][j] << endl;
                    }
                }
                file10.close();
            }
        }*/


		//  ##########   Computation of body forces Bf_u and Bf_v    ##########

        // interpolating the body force at the u forcing point (Bf_u)
        for(int i=1; i<n_u-1; i++)
        {
            for (int j=1; j<m_u-1; j++)
            {
                if (flag_u[i][j]!=3)	//Body forces at the fluid points are zero
                {
                    Bf_u[i][j] = 0.0;
                }
                else	// For the solid points i.e. the forcing points
                {
                    min_dis = 10.0;
                    for(int ctr=0;ctr<N_marker;ctr++)
                    {
                        distance = sqrt((x_s[ctr]-x_u[i][j])*(x_s[ctr]-x_u[i][j]) + (y_s[ctr]-y_u[i][j])*(y_s[ctr]-y_u[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = ctr;
                        }
                    }
                    if(min_dis<=0.000001)	 // For cases when the points overlap with the solid surface
                    {
                        u_body = u_pivot + (pitch_dot*(y_s[k1]-y_pivot));
                        U_intp = u_body;
                        Bf_u[i][j] = (U_intp-u_1[i][j])/dt;
                    }
                    else	// Caclulate the intersecting points SA
                    {
                        if (k1==0)
                        {
                            if((x_s[k1+1]-x_s[N_marker-1])==0)
                            {
                                x_intersect = x_s[N_marker-1];
                                y_intersect = y_u[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                            {
                                x_intersect = x_u[i][j];
                                y_intersect = y_s[N_marker-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1]);
                                x_intersect = (y_u[i][j]-y_s[N_marker-1]+slope*x_s[N_marker-1]+(1/slope)*x_u[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[N_marker-1] - slope*x_s[N_marker-1];
                            }
                        }
                        else if (k1==N_marker-1)
                        {
                            if((x_s[0]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y_u[i][j];
                            }
                            else if ((y_s[0]-y_s[k1-1])==0)
                            {
                                x_intersect = x_u[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1]);
                                x_intersect = (y_u[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x_u[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        else
                        {
                            if((x_s[k1+1]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y_u[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[k1-1])==0)
                            {
                                x_intersect = x_u[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1]);
                                x_intersect = (y_u[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x_u[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        u_body = u_pivot + (pitch_dot*(y_intersect-y_pivot));
                        x_mirror = x_intersect + (x_intersect-x_u[i][j]);	 // Mirror points x -cordinate
                        y_mirror = y_intersect + (y_intersect-y_u[i][j]); 	 // Mirror points y - coordinate
                        n1_i = c3+floor((y_mirror-y_u[c3][c1])/dy2);	 // Index location for lower left neighbour of mirror point - x index
                        n1_j = c1+floor((x_mirror-x_u[c3][c1])/dx2);	 // Index location for lower left neighbour of mirror Point - y index

						// To check if all four neighouring points are fluid points or not
						if(flag_u[n1_i][n1_j]==1 && flag_u[n1_i][n1_j+1]==1 && flag_u[n1_i+1][n1_j+1]==1 && flag_u[n1_i+1][n1_j]==1)
                        {
                            // Bilinear interpolation from surrounding points to the IA point
							u_mirror = (((((x_u[n1_i+1][n1_j+1]-x_mirror)*u_1[n1_i+1][n1_j]+(x_mirror-x_u[n1_i+1][n1_j])*u_1[n1_i+1][n1_j+1])/(x_u[n1_i+1][n1_j+1]-x_u[n1_i+1][n1_j]))*(y_mirror-y_u[n1_i][n1_j]))+((((x_u[n1_i][n1_j+1]-x_mirror)*u_1[n1_i][n1_j]+(x_mirror-x_u[n1_i][n1_j])*u_1[n1_i][n1_j+1])/(x_u[n1_i][n1_j+1]-x_u[n1_i][n1_j]))*(y_u[n1_i+1][n1_j]-y_mirror)))/(y_u[n1_i+1][n1_j]-y_u[n1_i][n1_j]);
                            dist_a = sqrt(pow((y_mirror-y_intersect),2)+pow((x_mirror-x_intersect),2));
                            dist_b = sqrt(pow((y_intersect-y_u[i][j]),2)+pow((x_intersect-x_u[i][j]),2));
                            U_intp = (1+(dist_b/dist_a))*u_body - (dist_b/dist_a)*u_mirror; 	// Velocity extrapolation
                            Bf_u[i][j] = (U_intp-u_1[i][j])/dt;		// Momentum forcing calculation
                        }
                        else	// change to a different IA point along the normal
                        {
                            if (k1==0)
                            {
                                if((x_s[k1+1]-x_s[N_marker-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            else if (k1==N_marker-1)
                            {
                                if((x_s[0]-x_s[k1-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[0]-y_s[k1-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            else
                            {
                                if((x_s[k1+1]-x_s[k1-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[k1+1]-y_s[k1-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            n1_i = c3+floor((y_mirror1-y_u[c3][c1])/dy2); 	 //neighbour index - x
                            n1_j = c1+floor((x_mirror1-x_u[c3][c1])/dx2); 	 //neighbour index - y
                            if(flag_u[n1_i][n1_j]==1 && flag_u[n1_i][n1_j+1]==1 && flag_u[n1_i+1][n1_j+1]==1 && flag_u[n1_i+1][n1_j]==1) 	// Check if all of them are inside fluid
                            {
                                u_mirror = (((((x_u[n1_i+1][n1_j+1]-x_mirror1)*u_1[n1_i+1][n1_j]+(x_mirror1-x_u[n1_i+1][n1_j])*u_1[n1_i+1][n1_j+1])/(x_u[n1_i+1][n1_j+1]-x_u[n1_i+1][n1_j]))*(y_mirror1-y_u[n1_i][n1_j]))+((((x_u[n1_i][n1_j+1]-x_mirror1)*u_1[n1_i][n1_j]+(x_mirror1-x_u[n1_i][n1_j])*u_1[n1_i][n1_j+1])/(x_u[n1_i][n1_j+1]-x_u[n1_i][n1_j]))*(y_u[n1_i+1][n1_j]-y_mirror1)))/(y_u[n1_i+1][n1_j]-y_u[n1_i][n1_j]);
                                dist_a = sqrt(pow((y_mirror1-y_intersect),2)+pow((x_mirror1-x_intersect),2));
                                dist_b = sqrt(pow((y_intersect-y_u[i][j]),2)+pow((x_intersect-x_u[i][j]),2));
                                U_intp = (1+(dist_b/dist_a))*u_body - (dist_b/dist_a)*u_mirror;
                                Bf_u[i][j] = (U_intp-u_1[i][j])/dt;
                            }
                            else
                            {
                                cout << "what will happen to me --> u" << endl;
                                /*cout << i << " --> " << j  << " --> " << x_u[i][j] << " --> " << y_u[i][j] << " --> " << "what will happen to me" << endl;
                                cout << x_intersect << " --> " << y_intersect  << " --> " << x_mirror << " --> " << y_mirror << " --> " << x_mirror1 << " -- >" << y_mirror1 << endl;
                                cout << n1_i << " --> " << n1_j << endl;
                                cout << flag_u[n1_i][n1_j] << " --> " << flag_u[n1_i][n1_j+1]  << " --> " << flag_u[n1_i+1][n1_j+1] << " --> " << flag_u[n1_i+1][n1_j] << endl;*/
                                flag_u[i][j] = 1;
                                flag_tu[i][j] = 0;
                                Bf_u[i][j] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        // interpolating the body force at the v forcing point (Bf_v)
        for(int i=1; i<n_v-1; i++)
        {
            for (int j=1; j<m_v-1; j++)
            {
                if (flag_v[i][j]!=3)
                {
                    Bf_v[i][j] = 0.0;
                }
                else
                {
                    min_dis = 10.0;
                    for(int ctr=0;ctr<N_marker;ctr++)
                    {
                        distance = sqrt((x_s[ctr]-x_v[i][j])*(x_s[ctr]-x_v[i][j]) + (y_s[ctr]-y_v[i][j])*(y_s[ctr]-y_v[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = ctr;
                        }
                    }
                    if(min_dis<=0.000001)
                    {
                        v_body = v_pivot + (-pitch_dot*(x_s[k1]-x_pivot));
                        V_intp = v_body;
                        Bf_v[i][j] = (V_intp-v_1[i][j])/dt;
                    }
                    else
                    {
                        if (k1==0)
                        {
                            if((x_s[k1+1]-x_s[N_marker-1])==0)
                            {
                                x_intersect = x_s[N_marker-1];
                                y_intersect = y_v[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                            {
                                x_intersect = x_v[i][j];
                                y_intersect = y_s[N_marker-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1]);
                                x_intersect = (y_v[i][j]-y_s[N_marker-1]+slope*x_s[N_marker-1]+(1/slope)*x_v[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[N_marker-1] - slope*x_s[N_marker-1];
                            }
                        }
                        else if (k1==N_marker-1)
                        {
                            if((x_s[0]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y_v[i][j];
                            }
                            else if ((y_s[0]-y_s[k1-1])==0)
                            {
                                x_intersect = x_v[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1]);
                                x_intersect = (y_v[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x_v[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        else
                        {
                            if((x_s[k1+1]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y_v[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[k1-1])==0)
                            {
                                x_intersect = x_v[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1]);
                                x_intersect = (y_v[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x_v[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        v_body = v_pivot + (-pitch_dot*(x_intersect-x_pivot));
                        x_mirror = x_intersect + (x_intersect-x_v[i][j]);
                        y_mirror = y_intersect + (y_intersect-y_v[i][j]);
                        n1_i = c3+floor((y_mirror-y_v[c3][c1])/dy2);
                        n1_j = c1+floor((x_mirror-x_v[c3][c1])/dx2);
                        if(flag_v[n1_i][n1_j]==1 && flag_v[n1_i][n1_j+1]==1 && flag_v[n1_i+1][n1_j+1]==1 && flag_v[n1_i+1][n1_j]==1)
                        {
                            v_mirror = (((((x_v[n1_i+1][n1_j+1]-x_mirror)*v_1[n1_i+1][n1_j]+(x_mirror-x_v[n1_i+1][n1_j])*v_1[n1_i+1][n1_j+1])/(x_v[n1_i+1][n1_j+1]-x_v[n1_i+1][n1_j]))*(y_mirror-y_v[n1_i][n1_j]))+((((x_v[n1_i][n1_j+1]-x_mirror)*v_1[n1_i][n1_j]+(x_mirror-x_v[n1_i][n1_j])*v_1[n1_i][n1_j+1])/(x_v[n1_i][n1_j+1]-x_v[n1_i][n1_j]))*(y_v[n1_i+1][n1_j]-y_mirror)))/(y_v[n1_i+1][n1_j]-y_v[n1_i][n1_j]);
                            dist_a = sqrt(pow((y_mirror-y_intersect),2)+pow((x_mirror-x_intersect),2));
                            dist_b = sqrt(pow((y_intersect-y_v[i][j]),2)+pow((x_intersect-x_v[i][j]),2));
                            V_intp = (1+(dist_b/dist_a))*v_body - (dist_b/dist_a)*v_mirror;
                            Bf_v[i][j] = (V_intp-v_1[i][j])/dt;
                        }
                        else
                        {
                            if (k1==0)
                            {
                                if((x_s[k1+1]-x_s[N_marker-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            else if (k1==N_marker-1)
                            {
                                if((x_s[0]-x_s[k1-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[0]-y_s[k1-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            else
                            {
                                if((x_s[k1+1]-x_s[k1-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[k1+1]-y_s[k1-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            n1_i = c3+floor((y_mirror1-y_v[c3][c1])/dy2);
                            n1_j = c1+floor((x_mirror1-x_v[c3][c1])/dx2);
                            if(flag_v[n1_i][n1_j]==1 && flag_v[n1_i][n1_j+1]==1 && flag_v[n1_i+1][n1_j+1]==1 && flag_v[n1_i+1][n1_j]==1)
                            {
                                v_mirror = (((((x_v[n1_i+1][n1_j+1]-x_mirror1)*v_1[n1_i+1][n1_j]+(x_mirror1-x_v[n1_i+1][n1_j])*v_1[n1_i+1][n1_j+1])/(x_v[n1_i+1][n1_j+1]-x_v[n1_i+1][n1_j]))*(y_mirror1-y_v[n1_i][n1_j]))+((((x_v[n1_i][n1_j+1]-x_mirror1)*v_1[n1_i][n1_j]+(x_mirror1-x_v[n1_i][n1_j])*v_1[n1_i][n1_j+1])/(x_v[n1_i][n1_j+1]-x_v[n1_i][n1_j]))*(y_v[n1_i+1][n1_j]-y_mirror1)))/(y_v[n1_i+1][n1_j]-y_v[n1_i][n1_j]);
                                dist_a = sqrt(pow((y_mirror1-y_intersect),2)+pow((x_mirror1-x_intersect),2));
                                dist_b = sqrt(pow((y_intersect-y_v[i][j]),2)+pow((x_intersect-x_v[i][j]),2));
                                V_intp = (1+(dist_b/dist_a))*v_body - (dist_b/dist_a)*v_mirror;
                                Bf_v[i][j] = (V_intp-v_1[i][j])/dt;
                            }
                            else
                            {
                                cout << "what will happen to me --> v" << endl;
                                flag_v[i][j] = 1;
                                flag_tv[i][j] = 0;
                                v[i][j] = v_1[i][j];
                                Bf_v[i][j] = 0.0;
                            }
                        }
                    }
                }
            }
        }


		//  ###########   Calculation of the RHS of x and y-momentum equation (Different from momentum forcing Bf_u)   ##########

		// Calculation of the RHS of x-momentum equation (Different from momentum forcing Bf_u)
		#pragma omp parallel for default(shared) private(ue_old,uw_old,un_old,us_old,vn_old,vs_old,ue,uw,un,us,vn,vs,A1,B1,C1) schedule(dynamic)
		for(int i=1; i<n_u-1; i++)
		{
			for (int j=1; j<m_u-1; j++)
			{
				ue_old = ((x[i][j+1]-x1[i][j+1])*u_old[i][j]+(x1[i][j+1]-x[i][j])*u_old[i][j+1])/(x[i][j+1]-x[i][j]);
				uw_old = ((x[i][j]-x1[i][j])*u_old[i][j-1]+(x1[i][j]-x[i][j-1])*u_old[i][j])/(x[i][j]-x[i][j-1]);
				un_old = ((y[i][j]-y1[i][j])*u_old[i+1][j]+(y1[i+1][j]-y[i][j])*u_old[i][j])/(y1[i+1][j]-y1[i][j]);
				us_old = ((y[i-1][j]-y1[i-1][j])*u_old[i][j]+(y1[i][j]-y[i-1][j])*u_old[i-1][j])/(y1[i][j]-y1[i-1][j]);
				vn_old = ((x1[i][j+1]-x[i][j])*v_old[i][j]+(x[i][j]-x1[i][j])*v_old[i][j+1])/(x1[i][j+1]-x1[i][j]);
				vs_old = ((x1[i][j+1]-x[i][j])*v_old[i-1][j]+(x[i][j]-x1[i][j])*v_old[i-1][j+1])/(x1[i][j+1]-x1[i][j]);
				ue = ((x[i][j+1]-x1[i][j+1])*u[i][j]+(x1[i][j+1]-x[i][j])*u[i][j+1])/(x[i][j+1]-x[i][j]);
				uw = ((x[i][j]-x1[i][j])*u[i][j-1]+(x1[i][j]-x[i][j-1])*u[i][j])/(x[i][j]-x[i][j-1]);
				un = ((y[i][j]-y1[i][j])*u[i+1][j]+(y1[i+1][j]-y[i][j])*u[i][j])/(y1[i+1][j]-y1[i][j]);
				us = ((y[i-1][j]-y1[i-1][j])*u[i][j]+(y1[i][j]-y[i-1][j])*u[i-1][j])/(y1[i][j]-y1[i-1][j]);
				vn = ((x1[i][j+1]-x[i][j])*v[i][j]+(x[i][j]-x1[i][j])*v[i][j+1])/(x1[i][j+1]-x1[i][j]);
				vs = ((x1[i][j+1]-x[i][j])*v[i-1][j]+(x[i][j]-x1[i][j])*v[i-1][j+1])/(x1[i][j+1]-x1[i][j]);
				A1 = 1.5*((dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw)+(dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us));
				B1 = 0.5*((dt/(x1[i][j+1]-x1[i][j]))*(ue_old*ue_old-uw_old*uw_old)+(dt/(y[i][j]-y[i-1][j]))*(vn_old*un_old-vs_old*us_old));
				C1 = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(((u[i][j+1]-u[i][j])/(x[i][j+1]-x[i][j]))-((u[i][j]-u[i][j-1])/(x[i][j]-x[i][j-1])))+(dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(((u[i+1][j]-u[i][j])/(y1[i+1][j]-y1[i][j]))-((u[i][j]-u[i-1][j])/(y1[i][j]-y1[i-1][j])));
				f_u[i][j] = u[i][j]-A1+B1+C1;
			}
		}

		// Calculation of the RHS term of the y-momentum equation (Different from momemntum forcing Bf_v)
		#pragma omp parallel for default(shared) private(ue_old,uw_old,ve_old,vw_old,vn_old,vs_old,ue,uw,ve,vw,vn,vs,A1,B1,C1) schedule(dynamic)
		for(int i=1; i<n_v-1; i++)
		{
			for (int j=1; j<m_v-1; j++)
			{
				ue_old = ((y[i][j]-y1[i][j])*u_old[i+1][j]+(y1[i+1][j]-y[i][j])*u_old[i][j])/(y1[i+1][j]-y1[i][j]);
				uw_old = ((y[i][j-1]-y1[i][j-1])*u_old[i+1][j-1]+(y1[i+1][j-1]-y[i][j-1])*u_old[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]);
				ve_old = ((x1[i][j+1]-x[i][j])*v_old[i][j]+(x[i][j]-x1[i][j])*v_old[i][j+1])/(x1[i][j+1]-x1[i][j]);
				vw_old = ((x1[i][j]-x[i][j-1])*v_old[i][j-1]+(x[i][j-1]-x1[i][j-1])*v_old[i][j])/(x1[i][j]-x1[i][j-1]);
				vn_old = ((y1[i+1][j]-y[i][j])*v_old[i+1][j]+(y[i+1][j]-y1[i+1][j])*v_old[i][j])/(y[i+1][j]-y[i][j]);
				vs_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
				ue = ((y[i][j]-y1[i][j])*u[i+1][j]+(y1[i+1][j]-y[i][j])*u[i][j])/(y1[i+1][j]-y1[i][j]);
				uw = ((y[i][j-1]-y1[i][j-1])*u[i+1][j-1]+(y1[i+1][j-1]-y[i][j-1])*u[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]);
				ve = ((x1[i][j+1]-x[i][j])*v[i][j]+(x[i][j]-x1[i][j])*v[i][j+1])/(x1[i][j+1]-x1[i][j]);
				vw = ((x1[i][j]-x[i][j-1])*v[i][j-1]+(x[i][j-1]-x1[i][j-1])*v[i][j])/(x1[i][j]-x1[i][j-1]);
				vn = ((y1[i+1][j]-y[i][j])*v[i+1][j]+(y[i+1][j]-y1[i+1][j])*v[i][j])/(y[i+1][j]-y[i][j]);
				vs = ((y1[i][j]-y[i-1][j])*v[i][j]+(y[i][j]-y1[i][j])*v[i-1][j])/(y[i][j]-y[i-1][j]);
				A1 = 1.5*((dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw)+(dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs));
				B1 = 0.5*((dt/(x[i][j]-x[i][j-1]))*(ue_old*ve_old-uw_old*vw_old)+(dt/(y1[i+1][j]-y1[i][j]))*(vn_old*vn_old-vs_old*vs_old));
				C1 = (dt/(2.0*Re*(x[i][j]-x[i][j-1])))*(((v[i][j+1]-v[i][j])/(x1[i][j+1]-x1[i][j]))-((v[i][j]-v[i][j-1])/(x1[i][j]-x1[i][j-1])))+(dt/(2.0*Re*(y1[i+1][j]-y1[i][j])))*(((v[i+1][j]-v[i][j])/(y[i+1][j]-y[i][j]))-((v[i][j]-v[i-1][j])/(y[i][j]-y[i-1][j])));
				f_v[i][j] = v[i][j]-A1+B1+C1;
			}
		}


		// Terms from the previous time step needed for force calculation

		// For drag (A volume integration term temp_D is calculated using reduction)
		temp_D=0.0;
        #pragma omp parallel for default(shared) private(uf_old,uf_n,uf_s,un,us) reduction(+:temp_D) schedule(dynamic)
        for(int i=c3;i<c4;i++)
        {
            for (int j=c1;j<c2;j++)
            {
                if(flag1[i][j]==0)
                {
                    uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
                    uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
                    uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
                    un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
                    us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
                    temp_D += ((-0.5)*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
                }
            }
        }

		// For Lift (A volume integration term temp_L is calculated using reduction)
		temp_L=0.0;
        #pragma omp parallel for default(shared) private(vf_old,vf_e,vf_w,ve,vw) reduction(+:temp_L) schedule(dynamic)
        for(int i=c3;i<c4;i++)
        {
            for (int j=c1;j<c2;j++)
            {
                if(flag1[i][j]==0)
                {
                    vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
                    vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
                    vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
                    ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
                    vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
                    temp_L += ((-0.5)*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
                }
            }
        }


        //storing the values of u and v as old values. u and v will now be modified
		for (int i=0; i<n_u; i++)
        {
            for(int j=0; j<m_u; j++)
            {
                u_old[i][j]= u[i][j];
            }
        }
		for (int i=0; i<n_v; i++)
        {
            for(int j=0; j<m_v; j++)
            {
                v_old[i][j]= v[i][j];
            }
        }


        /*ofstream filew6("a_u.txt");
        if (filew6.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    filew6 << a_u[i][j] << " ";
                }
                filew6 << endl;
            }
            filew6.close();
        }
        ofstream filew7("b_u.txt");
        if (filew7.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    filew7 << b_u[i][j] << " ";
                }
                filew7 << endl;
            }
            filew7.close();
        }
        ofstream filew8("c_u.txt");
        if (filew8.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    filew8 << c_u[i][j] << " ";
                }
                filew8 << endl;
            }
            filew8.close();
        }
        ofstream filew9("d_u.txt");
        if (filew9.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    filew9 << d_u[i][j] << " ";
                }
                filew9 << endl;
            }
            filew9.close();
        }
        ofstream filew10("e_u.txt");
        if (filew10.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    filew10 << e_u[i][j] << " ";
                }
                filew10 << endl;
            }
            filew10.close();
        }
        ofstream filew11("f_u.txt");
        if (filew11.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    filew11 << f_u[i][j] << " ";
                }
                filew11 << endl;
            }
            filew11.close();
        }
        ofstream filew17("f_v.txt");
        if (filew17.is_open())
        {
            for(int i=0; i<n_v; i++)
            {
                for(int j=0; j<m_v; j++)
                {
                    filew17 << f_v[i][j] << " ";
                }
                filew17 << endl;
            }
            filew17.close();
        }*/


		for(int i=0; i<n_u; i++)
        {
            u[i][0]= u_inlet[k];
        }


		// Calculates net mass flow rate at the inlet
		mass_in = 0.0;
		#pragma omp parallel for default(shared) reduction(+:mass_in) schedule(dynamic)
        for(int i=1;i<n_u-1;i++)
        {
            mass_in += (u[i][0]*(y[i][m_u-1]-y[i-1][m_u-1]));
        }


		// SOlving momentum equations using Gauss Seidel Successive Over relaxation with red black tagging
		solverUVGaussSeidelSOR(m, m_u, n_u, m_v, n_v, (float*)a_u, (float*)b_u, (float*)c_u, (float*)d_u, (float*)e_u, (float*)f_u, (float*)Bf_u, (float*)a_v, (float*)b_v, (float*)c_v, (float*)d_v, (float*)e_v, (float*)f_v, (float*)Bf_v, (float*)x1, (float*)y1, (float*)y, (float*)pressure, mass_in, domain_width, (float*)uA_old, (float*)vA_old, (float*)u, (float*)v, (float*)u_red, (float*)v_red, (float*)u_black, (float*)v_black, (float*)u_red_old, (float*)v_red_old, (float*)u_black_old, (float*)v_black_old, (float*)max_error1, (float*)max_error2);


        //   #############      Calculation of f_p term of the pressure correction equation     ############
        for(int i=1; i<n-1; i++)
        {
            for (int j=1; j<m-1; j++)
            {
                if(flag_u[i][j]==1 && flag_u[i][j-1]==1 && flag_v[i][j]==1 && flag_v[i-1][j]==1)
                    f_p[i][j] = -(1.0/dt)*(((u[i][j]-u[i][j-1])/(x[i][j]-x[i][j-1]))+((v[i][j]-v[i-1][j])/(y[i][j]-y[i-1][j])));
                else if (flag_u[i][j]==0 && flag_u[i][j-1]==0 && flag_v[i][j]==0 && flag_v[i-1][j]==0)
                    f_p[i][j] = 0.0;
                else
                {
                    min_dis = 10;
                    for(int ctr=0;ctr<N_marker;ctr++)
                    {
                        distance = sqrt((x_s[ctr]-x1[i][j])*(x_s[ctr]-x1[i][j]) + (y_s[ctr]-y1[i][j])*(y_s[ctr]-y1[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = ctr;
                        }
                    }
                    if(min_dis<=0.000001)
                    {
                        u_body = u_pivot + (pitch_dot*(y_s[k1]-y_pivot));
                        v_body = v_pivot + (-pitch_dot*(x_s[k1]-x_pivot));
                    }
                    else
                    {
                        if (k1==0)
                        {
                            if((x_s[k1+1]-x_s[N_marker-1])==0)
                            {
                                x_intersect = x_s[N_marker-1];
                                y_intersect = y1[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                            {
                                x_intersect = x1[i][j];
                                y_intersect = y_s[N_marker-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1]);
                                x_intersect = (y1[i][j]-y_s[N_marker-1]+slope*x_s[N_marker-1]+(1/slope)*x1[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[N_marker-1] - slope*x_s[N_marker-1];
                            }
                        }
                        else if (k1==N_marker-1)
                        {
                            if((x_s[0]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y1[i][j];
                            }
                            else if ((y_s[0]-y_s[k1-1])==0)
                            {
                                x_intersect = x1[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1]);
                                x_intersect = (y1[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x1[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        else
                        {
                            if((x_s[k1+1]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y1[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[k1-1])==0)
                            {
                                x_intersect = x1[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1]);
                                x_intersect = (y1[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x1[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        u_body = u_pivot + (pitch_dot*(y_intersect-y_pivot));
                        v_body = v_pivot + (-pitch_dot*(x_intersect-x_pivot));
                    }
                    f_p[i][j] = -(1.0/dt)*((((u[i][j]-u[i][j-1])/(x[i][j]-x[i][j-1]))+((v[i][j]-v[i-1][j])/(y[i][j]-y[i-1][j]))) - (((flag_tu[i][j]*(u[i][j]-u_body))/(x[i][j]-x[i][j-1]))+(-(flag_tu[i][j-1]*(u[i][j-1]-u_body))/(x[i][j]-x[i][j-1]))+((flag_tv[i][j]*(v[i][j]-v_body))/(y[i][j]-y[i-1][j]))+(-(flag_tv[i-1][j]*(v[i-1][j]-v_body))/(y[i][j]-y[i-1][j]))));
                }
            }
        }


         /*   ofstream filew23("f_p.txt");
        if (filew23.is_open())
        {
            for(int i=0; i<n; i++)
            {
                for(int j=0; j<m; j++)
                {
                    filew23 << f_p[i][j]*dt << " ";
                }
                filew23 << endl;
            }
            filew23.close();
        }
        ofstream filew1("Cavity_u.txt");
        if (filew1.is_open())
        {
            for(int i=0; i<n_u; i++)
            {
                for(int j=0; j<m_u; j++)
                {
                    filew1 << u[i][j] << " ";
                }
                filew1 << endl;
            }
            filew1.close();
        }
        ofstream filew2("Cavity_v.txt");
        if (filew2.is_open())
        {
            for(int i=0; i<n_v; i++)
            {
                for(int j=0; j<m_v; j++)
                {
                    filew2 << v[i][j] << " ";
                }
                filew2 << endl;
            }
            filew2.close();
        }*/
        /*big = 0;
        for (int i=1;i<n-1;i++)
        {
            for (int j=1;j<m-1;j++)
            {
                error=abs(f_p[i][j]);
                if (error>big)
                    big = error;
            }
        }
        residual_flux=big*dt;*/


		// Solving pressure correction using Gauss Seidel iterations with red-black SOR:
		solverPGaussSeidelSOR(m, n, (float*)a_p, (float*)b_p, (float*)c_p, (float*)d_p, (float*)e_p, (float*)f_p, (float*)p_c, (float*)pA_old, (float*)p_c_red, (float*)p_c_black ,(float*)p_c_red_old, (float*)p_c_black_old, (float*)max_error3);


        //pressure correction
		#pragma omp parallel for default(shared) schedule(dynamic)
        for(int i=1; i<n-1; i++)
        {
            for (int j=1; j<m-1; j++)
            {
                pressure[i][j] = pressure[i][j]+ 0.0001*(p_c[i][j] - (dt/(2*Re))*(((((p_c[i][j+1]-p_c[i][j])/(x1[i][j+1]-x1[i][j]))-((p_c[i][j]-p_c[i][j-1])/(x1[i][j]-x1[i][j-1])))/(x[i][j]-x[i][j-1]))+((((p_c[i+1][j]-p_c[i][j])/(y1[i+1][j]-y1[i][j]))-((p_c[i][j]-p_c[i-1][j])/(y1[i][j]-y1[i-1][j])))/(y[i][j]-y[i-1][j]))));
            }
        }

        //u velocity correction
		#pragma omp parallel for default(shared) schedule(dynamic)
        for(int i=1; i<n_u-1; i++)
        {
            for (int j=1; j<m_u-1; j++)
            {
                u[i][j] = u[i][j]-dt*((p_c[i][j+1]-p_c[i][j])/(x1[i][j+1]-x1[i][j]));
            }
        }

        //v velocity correction
		#pragma omp parallel for default(shared) schedule(dynamic)
        for(int i=1; i<n_v-1; i++)
        {
            for (int j=1; j<m_v-1; j++)
            {
                v[i][j] = v[i][j]-dt*((p_c[i+1][j]-p_c[i][j])/(y1[i+1][j]-y1[i][j]));
            }
        }


		// #####   u- velocity correction at the outlet boundary to satisfy mass conservation for the entire coputational domain     #####
		for (int p=0; p<4; p++)
        {
            for(int i=0; i<n_u; i++)
            {
                for (int j=1; j<m_u; j++)
                {
                    if(i==0) // Bottom
                    {
                        u[i][j] = u[i+1][j];
                    }
                    else if(i==n_u-1) // Top
                    {
                        u[i][j] = u[i-1][j];
                    }
                    else if(j==m_u-1) // Right
                    {
                        u[i][j] = u[i][j-1];
                    }
                }
            }

			// Check the mass flow at the outlet and perform correction
            for(int ctr=1;ctr<=40;ctr++)
            {
                mass_out = 0.0;
				#pragma omp parallel for default(shared) reduction(+:mass_out) schedule(dynamic)
                for(int i=1;i<n_u-1;i++)
                {
                    mass_out += (u[i][m_u-1]*(y[i][m_u-1]-y[i-1][m_u-1]));
                }

				#pragma omp parallel for default(shared) schedule(dynamic)
                for(int i=1; i<n_u-1; i++)
                {
                    u[i][m_u-1] = u[i][m_u-1]+((mass_in-mass_out)/domain_width);
                }
            }
        }
		for(int i=1; i<n_v-1; i++)
			v[i][m_v-1] = v[i][m_v-2];


		// ###############    Actual fluid force calculation     ###############

		// Drag calculation

        temp=0.0;
        #pragma omp parallel for default(shared) private(uf,uf_old,uf_n,uf_s,un,us) reduction(+:temp) schedule(dynamic)
        for(int i=c3;i<c4;i++)
        {
            for (int j=c1;j<c2;j++)
            {
                if(flag1[i][j]==0)
                {
                    uf = ((x1[i][j]-x[i][j-1])*u[i][j]+(x[i][j]-x1[i][j])*u[i][j-1])/(x[i][j]-x[i][j-1]);
                    uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
                    uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
                    uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
                    un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
                    us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
                    temp += (((uf-uf_old)/dt)+1.5*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
                }
            }
        }
        temp1 = 0.0;
		#pragma omp parallel for default(shared) reduction(+:temp1) schedule(dynamic)
        for(int i=c3;i<c4;i++)
        {
            for(int j=c1;j<c2;j++)
            {
                if (flag_u[i][j] == 3)
                    temp1 += Bf_u[i][j]*(x1[i][j+1]-x1[i][j])*(y[i][j]-y[i-1][j]);
            }
        }
        drag3 = 2.0*(-temp1+temp+temp_D); //gives drag coeff Cd = D/(0.5*rho*u^2*(d*1))

		// Lift calculation

        temp=0.0;
        #pragma omp parallel for default(shared) private(vf,vf_old,vf_e,vf_w,ve,vw) reduction(+:temp) schedule(dynamic)
        for(int i=c3;i<c4;i++)
        {
            for (int j=c1;j<c2;j++)
            {
                if(flag1[i][j]==0)
                {
                    vf = ((y1[i][j]-y[i-1][j])*v[i][j]+(y[i][j]-y1[i][j])*v[i-1][j])/(y[i][j]-y[i-1][j]);
                    vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
                    vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
                    vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
                    ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
                    vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
                    temp += (((vf-vf_old)/dt)+1.5*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
                }
            }
        }
        temp1 = 0.0;
		#pragma omp parallel for default(shared) reduction(+:temp1) schedule(dynamic)
        for(int i=c3;i<c4;i++)
        {
            for(int j=c1;j<c2;j++)
            {
                if (flag_v[i][j] == 3)
                    temp1 += Bf_v[i][j]*(x[i][j]-x[i][j-1])*(y1[i+1][j]-y1[i][j]);
            }
        }
        lift3 = 2.0*(-temp1+temp+temp_L); //gives lift coeff Cl = L/(0.5*rho*u^2*(d*1))


		// Printing/ saving the lift and drag data at every time step
		if (file5c.is_open())
        {
            file5c << fixed << setprecision(6) << time1 << " " << pitch << " " << pitch_dot << " " << drag3 << " " << lift3 << endl;
            //file5c.close();
        }

        /*if (file6.is_open())
        {
            file6 << residual_flux << endl;
            //file6.close();
        }*/
        /*if (file6a.is_open())
        {
            file6a << time1 << " " << ((x1[probe_i][probe_j]-x[probe_i][probe_j-1])/(x[probe_i][probe_j]-x[probe_i][probe_j-1]))*u[probe_i][probe_j]+((x[probe_i][probe_j]-x1[probe_i][probe_j])/(x[probe_i][probe_j]-x[probe_i][probe_j-1]))*u[probe_i][probe_j-1] << " " << ((y1[probe_i][probe_j]-y[probe_i-1][probe_j])/(y[probe_i][probe_j]-y[probe_i-1][probe_j]))*v[probe_i][probe_j]+((y[probe_i][probe_j]-y1[probe_i][probe_j])/(y[probe_i][probe_j]-y[probe_i-1][probe_j]))*v[probe_i-1][probe_j] << endl;
            //file6a.close();
        }*/

        cout << k << endl;


		//   ##########  File handling. Printing out all the necessary fields    ##########

        if ( (k%writeInterval) == 0)
        {
            /*ofstream filew1("Cavity_u.txt");
            if (filew1.is_open())
            {
                for(int i=0; i<n_u; i++)
                {
                    for(int j=0; j<m_u; j++)
                    {
                        filew1 << u[i][j] << " ";
                    }
                    filew1 << endl;
                }
                filew1.close();
            }

            ofstream filew2("Cavity_v.txt");
            if (filew2.is_open())
            {
                for(int i=0; i<n_v; i++)
                {
                    for(int j=0; j<m_v; j++)
                    {
                        filew2 << v[i][j] << " ";
                    }
                    filew2 << endl;
                }
                filew2.close();
            }
            ofstream filew3("pressure.txt");
            if (filew3.is_open())
            {
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<m; j++)
                    {
                        filew3 << pressure[i][j] << " ";
                    }
                    filew3 << endl;
                }
                filew3.close();
            }
            */
            stringstream ss;
            ss << k;
            fileName = ss.str();

            fileName1 = "Cavity_u"+fileName+".txt";
            ofstream filew1(fileName1.c_str());

            fileName4 = "Cavity_v"+fileName+".txt";
            ofstream filew2(fileName4.c_str());

            fileName7 = "Cavity_u_old"+fileName+".txt";
            ofstream filew1a(fileName7.c_str());

            fileName8 = "Cavity_v_old"+fileName+".txt";
            ofstream filew2a(fileName8.c_str());

            fileName9 = "pressure"+fileName+".txt";
            ofstream filew3(fileName9.c_str());

            fileName10 = "p_c"+fileName+".txt";
            ofstream filew3a(fileName10.c_str());

            #pragma omp parallel default(shared)
            {
                #pragma omp single
                if (filew1.is_open())
                {
                    for(int i=0; i<n_u; i++)
                    {
                        for(int j=0; j<m_u; j++)
                        {
                            filew1 << u[i][j] << endl;
                        }
                    }
                    filew1.close();
                }
                #pragma omp single
                if (filew2.is_open())
                {
                    for(int i=0; i<n_v; i++)
                    {
                        for(int j=0; j<m_v; j++)
                        {
                            filew2 << v[i][j] << endl;
                        }
                    }
                    filew2.close();
                }
                #pragma omp single
                if (filew1a.is_open())
                {
                    for(int i=0; i<n_u; i++)
                    {
                        for(int j=0; j<m_u; j++)
                        {
                            filew1a << u_old[i][j] << endl;
                        }
                    }
                    filew1.close();
                }
                #pragma omp single
                if (filew2a.is_open())
                {
                    for(int i=0; i<n_v; i++)
                    {
                        for(int j=0; j<m_v; j++)
                        {
                            filew2a << v_old[i][j] << endl;
                        }
                    }
                    filew2a.close();
                }
                #pragma omp single
                if (filew3.is_open())
                {
                    for(int i=0; i<n; i++)
                    {
                        for(int j=0; j<m; j++)
                        {
                            filew3 << pressure[i][j] << endl;
                        }
                    }
                    filew3.close();
                }
                #pragma omp single
                if (filew3a.is_open())
                {
                    for(int i=0; i<n; i++)
                    {
                        for(int j=0; j<m; j++)
                        {
                            filew3a << p_c[i][j] << endl;
                        }
                    }
                    filew3a.close();
                }
            }

			// All data for every snap in one file
            for (int i=0;i<n;i++)
            {
                for (int j=0;j<m;j++)
                {
                    if (i==0 && j==0)//bottom row
                    {
                        uf1[i][j] = u[i][j];
                        vf1[i][j] = v[i][j];
                        pressure[i][j] = 2*pressure[i+1][j]-pressure[i+2][j];
                    }
                    else if (i==0 && j>0 && j <m-1)
                    {
                        uf1[i][j] = ((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u[i][j]+((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u[i][j-1];
                        vf1[i][j] = v[i][j];
                        pressure[i][j] = 2*pressure[i+1][j]-pressure[i+2][j];
                    }
                    else if (i==0 && j==m-1)
                    {
                        uf1[i][j] = u[i][j-1];
                        vf1[i][j] = v[i][j];
                        pressure[i][j] = 2*pressure[i+1][j]-pressure[i+2][j];
                    }
                    else if (i>0 && i<n-1 && j==0)//left column
                    {
                        uf1[i][j] = u[i][j];
                        vf1[i][j] = ((y1[i][0]-y[i-1][0])/(y[i][0]-y[i-1][0]))*v[i][j]+((y[i][0]-y1[i][0])/(y[i][0]-y[i-1][0]))*v[i-1][j];
                        pressure[i][j] = 2*pressure[i][j+1]-pressure[i][j+2];
                    }
                    else if (i>0 && i<n-1 && j==m-1)//right column
                    {
                        uf1[i][j] = u[i][j-1];
                        vf1[i][j] = ((y1[i][0]-y[i-1][0])/(y[i][0]-y[i-1][0]))*v[i][j]+((y[i][0]-y1[i][0])/(y[i][0]-y[i-1][0]))*v[i-1][j];
                        pressure[i][j] = 2*pressure[i][j-1]-pressure[i][j-2];
                    }
                    else if (i==n-1 && j==0)
                    {
                        uf1[i][j] = u[i][j];
                        vf1[i][j] = v[i-1][j];
                        pressure[i][j] = 2*pressure[i-1][j]-pressure[i-2][j];
                    }
                    else if (i==n-1 && j>0 && j <m-1)//top row
                    {
                        uf1[i][j] = ((x1[0][j]-x[0][j-1])/(x[0][j]-x[0][j-1]))*u[i][j]+((x[0][j]-x1[0][j])/(x[0][j]-x[0][j-1]))*u[i][j-1];
                        vf1[i][j] = v[i-1][j];
                        pressure[i][j] = 2*pressure[i-1][j]-pressure[i-2][j];
                    }
                    else if (i==n-1 && j==m-1)
                    {
                        uf1[i][j] = u[i][j-1];
                        vf1[i][j] = v[i-1][j];
                        pressure[i][j] = 2*pressure[i-1][j]-pressure[i-2][j];
                    }
                    else
                    {
                        uf1[i][j] = ((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u[i][j]+((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u[i][j-1];
                        vf1[i][j] = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v[i-1][j];
                        pressure[i][j] = pressure[i][j];
                    }
                }
            }
            fileName = "AllinOne"+fileName+".dat";
            ofstream file7(fileName.c_str());
            if (file7.is_open())
            {
                file7 << "ZONE T=DATA I=" << m << " J=" << n << endl;
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<m; j++)
                    {
                        file7 << x1[i][j] << " " << y1[i][j] << " " << uf1[i][j] << " " << vf1[i][j] << " " << pressure[i][j] << endl;
                    }
                }
                file7.close();
            }
        }
    }	// Time marching ends here

		/*ofstream file8("error_u.txt");
        if (file8.is_open())
        {
            for(int i=0; i<GSite; i++)
            {
                file8 << max_error1[i] << endl;
            }
            file8.close();
        }
        ofstream file9("error_v.txt");
        if (file9.is_open())
        {
            for(int i=0; i<GSite; i++)
            {
                file9 << max_error2[i] << endl;
            }
            file9.close();
        }
        ofstream file10("error_p_c.txt");
        if (file10.is_open())
        {
            for(int i=0; i<GSite2; i++)
            {
                file10 << max_error3[i] << endl;
            }
            file10.close();
        }*/
        /*for (int ctr=0; ctr<N_marker; ctr++)
        {
            min_dis=10;
            for (int i=1;i<n-1;i++)
            {
                for (int j=1;j<m-1;j++)
                {
                    if(j>c1 && j<c2 && i>c3 && i<c4)
                    {
                        if(flag1[i][j]==1)
                        {
                            distance = sqrt(pow((x_s[ctr]-x1[i][j]),2)+pow((y_s[ctr]-y1[i][j]),2));
                            if (distance<min_dis)
                            {
                                min_dis = distance;
                                i_min = i;
                                j_min = j;
                            }
                        }
                    }
                }
            }
            surface_pr[ctr] = pressure[i_min][j_min];
        }
        ofstream file11("surface_pressure.txt");
        if (file11.is_open())
        {
            for(int i=0; i<N_marker; i++)
            {
                file11 << 2*surface_pr[i] << endl;
            }
            file11 << 2*surface_pr[0] << endl;
            file11.close();
        }*/

	cout << "\nCongratulation. The simulation has completed successfully." << endl;

	return 0;
}
