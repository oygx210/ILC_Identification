    function F = convert_state_space_to_matrix(AD_in, BD_in, CD_in, time)
        
        f_vec = zeros(time.nt, 1);
        
        for k3 = 1:time.nt
            
            f_vec(k3) = CD_in*AD_in^(k3-1)*BD_in;
            
        end
        
        F = convert_from_vector_to_matrix( f_vec, time );
    end