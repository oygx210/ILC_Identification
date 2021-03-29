function g_m_update = update_fmincon_gm(sys, data, time)
                
                                options = optimoptions('fmincon','Display',...
                                    'iter-detailed','Algorithm','sqp', ...
                                    'MaxFunctionEvaluations', 1e5,...
                                    'StepTolerance', 1e-5,...
                                    'OptimalityTolerance', 1e-5);
                
                
                                g_m_update = fmincon( @costIdentGS, data.g_m_init , [],[],[],[],[],[], [] ,  options);
                
                
                                
                                    function cost = costIdentGS( g_m_estimate )
        
        % regularisation matrices
        Q_ident     = diag( linspace( 1, 1, time.nt ) );
        S_delta     = 1e1*diag( linspace( 1, 10, time.nt-1 ) );
        G_abs       = 1e0*diag( linspace( 0, 10, time.nt ) );
        
        % transform the input to a matrix
        G_M_estimate = convert_from_vector_to_matrix( g_m_estimate, time );
        
        e1 = (data.y_LTI - G_M_estimate * data.u);
        e2 = g_m_estimate(2:end) - g_m_estimate(1:end-1);
        e3 = g_m_estimate;
        % e4 = (err_curr - (eye(nt) - gamma*G_M_estimate*G_M^-1)*err_prev);
        e4 = 0*(sys.G_M *(data.e_current - data.e_prev) + data.gamma*G_M_estimate*data.e_prev);
        
        cost = e1'*Q_ident*e1 + e2'*S_delta*e2 + e3'*G_abs*e3 + e4'*Q_ident*e4;
        
    end
                
                
end
                