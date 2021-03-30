function g_m_update = update_fmincon_gm_io(data, time)

%% settings
% regularisation matrices
Q_ident     = diag( linspace( 1, 1, time.nt ) );
S_delta     = 1e-1*eye(time.nt -1);
G_abs       = 1e-2*diag( linspace( 1, 1, time.nt ) );
        
        
options = optimoptions('fmincon','Display',...
    'iter-detailed','Algorithm','sqp', ...
    'MaxFunctionEvaluations', 1e5,...
    'StepTolerance', 1e-5,...
    'OptimalityTolerance', 1e-5);


%% algorithm

g_m_update = fmincon( @costIdentConv, data.g_m_init , [],[],[],[],[],[], [] ,  options);

    % fmincon cost function for conventional model identification
    function cost = costIdentConv( g_m_estimate )
        
        % transform the input to a matrix
        G_M_estimate = convert_from_vector_to_matrix( g_m_estimate, time );
        
        % calculate the error
        error_est = (data.y_LTI - G_M_estimate * data.u);
        
        % calculate impulse change rate for regularization
        g_delta_reg = g_m_estimate(2:end) - g_m_estimate(1:end-1);
        
        % calculate absolute impulse for regularization
        g_reg       = g_m_estimate;
        
        % calculate the cost
        cost = error_est'*Q_ident*error_est + g_delta_reg'*S_delta*g_delta_reg + g_reg'*G_abs*g_reg;
        
    end

end

