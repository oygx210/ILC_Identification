% generate a toeplitz matrix from vector
    function H = convert_from_vector_to_matrix( h_vec , time)
        
        first_row_h = [h_vec(1), zeros(1,time.nt-1)];
        first_colum_h = h_vec;
        
        H = toeplitz(first_colum_h, first_row_h);
        
    end