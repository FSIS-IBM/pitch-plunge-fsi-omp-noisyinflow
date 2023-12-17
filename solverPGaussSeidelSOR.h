extern int GSite2;
extern float alpha_SOR, tol;
void solverPGaussSeidelSOR(int m, int n, float* a_p, float* b_p, float* c_p, float* d_p, float* e_p, float* f_p, float* p_c, float* pA_old, float* p_c_red, float* p_c_black, float* p_c_red_old, float* p_c_black_old, float* max_error3);