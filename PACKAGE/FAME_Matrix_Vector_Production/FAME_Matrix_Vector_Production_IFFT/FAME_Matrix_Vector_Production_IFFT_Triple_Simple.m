function vec_y = FAME_Matrix_Vector_Production_IFFT_Triple_Simple(vec_x,Nx,Ny,Nz,N,D_k)

        vec_x_1 = vec_x(       1 :   N, 1 );
        vec_x_2 = vec_x(   N + 1 : 2*N, 1 );
        vec_x_3 = vec_x( 2*N + 1 : 3*N, 1 );
        % Reshape the input vector to do the inverse fast Fourier transform
        mat_x_1 = reshape(vec_x_1,Nx,Ny,Nz);
        mat_x_2 = reshape(vec_x_2,Nx,Ny,Nz);
        mat_x_3 = reshape(vec_x_3,Nx,Ny,Nz);
        % Inverse Fast Fourier Transform 
        vec_y_1 = reshape(ifftn(mat_x_1),N,1);
        vec_y_2 = reshape(ifftn(mat_x_2),N,1);
        vec_y_3 = reshape(ifftn(mat_x_3),N,1);
        %
        vec_y_1 = D_k.*vec_y_1;
        vec_y_2 = D_k.*vec_y_2;
        vec_y_3 = D_k.*vec_y_3;
        % Assemble output vector
        vec_y = sqrt(N)*[vec_y_1;vec_y_2;vec_y_3];            
end