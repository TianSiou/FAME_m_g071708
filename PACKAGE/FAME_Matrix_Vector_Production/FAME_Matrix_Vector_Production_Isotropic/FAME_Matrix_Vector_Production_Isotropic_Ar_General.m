function vec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_General(vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz)
    vec_y = FAME_Matrix_Vector_Production_Qr_General(vec_x, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz, 'normal');
  
%     vec_y = invB_eps.*vec_y;
    vec_y = FAME_Matrix_Vector_Production_invB_Isotropic(vec_y, B);

    vec_y = FAME_Matrix_Vector_Production_Qr_General(vec_y, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz, 'hermitian');
end