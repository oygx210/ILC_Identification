function g_m_update = update_fmincon_gm(sys, data, time)

%% settings
% regularisation matrices
Q_ident     = diag( linspace( 1, 1, time.nt ) );
S_delta     = 1e1*diag( linspace( 1, 10, time.nt-1 ) );
G_abs       = 1e0*diag( linspace( 0, 10, time.nt ) );


options = optimoptions('fmincon','Display',...
    'iter-detailed','Algorithm','sqp', ...
    'MaxFunctionEvaluations', 1e5,...
    'StepTolerance', 1e-5,...
    'OptimalityTolerance', 1e-5);



%% algorithm

% copy data
G_M         = sys.G_M;
e_current   = data.e_current;
e_prev      = data.e_prev;
gamma       = data.gamma;

g_m_update  = fmincon( @costIdentGS, data.g_m_init , [],[],[],[],[],[], [] ,  options);

    % fmincon cost function for conventional model identification
    function cost = costIdentGS( g_m_estimate )
        
        % transform the input to a matrix
        G_M_estimate = convert_from_vector_to_matrix( g_m_estimate, time );
        
        % calculate the error
        error_est = (G_M *(e_current - e_prev) + gamma*G_M_estimate*e_prev);
        
        % calculate impulse change rate for regularization
        g_delta_reg = g_m_estimate(2:end) - g_m_estimate(1:end-1);
        
        % calculate absolute impulse for regularization
        g_reg       = g_m_estimate;
        
        % calculate the cost
        cost = error_est'*Q_ident*error_est + g_delta_reg'*S_delta*g_delta_reg + g_reg'*G_abs*g_reg;
        
    end

end
