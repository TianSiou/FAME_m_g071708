function vec_y = FAME_Matrix_Vector_Production_Pr_Simple_fem(vec_x, Nx, Ny, Nz, N, Pi_Pr, Pi_Prs, D_k, D_ks, Lambda_fem, mode)
    if strcmp(mode,'normal') == 1
        vec_y = Pi_Pr*vec_x;
        vec_y = Lambda_fem.*vec_y;
        vec_y = FAME_Matrix_Vector_Production_IFFT_Triple_Simple(vec_y,Nx,Ny,Nz,N,D_k);
    elseif strcmp(mode,'hermitian') == 1
        vec_y = FAME_Matrix_Vector_Production_FFT_Triple_Simple(vec_x,Nx,Ny,Nz,N,D_ks);
        vec_y = conj(Lambda_fem).*vec_y;
        vec_y = Pi_Prs*vec_y;
    end
end