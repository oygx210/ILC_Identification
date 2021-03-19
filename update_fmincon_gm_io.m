% if model update via nonlinear optimisation is used
function g_m_update = update_fmincon_gm_io(data, time)

options = optimoptions('fmincon','Display',...
    'iter-detailed','Algorithm','sqp', ...
    'MaxFunctionEvaluations', 1e5,...
    'StepTolerance', 1e-5,...
    'OptimalityTolerance', 1e-5);


g_m_update = fmincon( @costIdentConv, data.g_m_init , [],[],[],[],[],[], [] ,  options);


% fmincon cost function for conventional model identification
    function cost = costIdentConv( g_m_estimate )
        
        % regularisation matrices
        Q_ident     = diag( linspace( 1, 1, time.nt ) );
        S_delta     = 1e-1*eye(time.nt -1);
        G_abs       = 1e-2*diag( linspace( 1, 1, time.nt ) );
        
        % transform the input to a matrix
        G_M_estimate = convert_from_vector_to_matrix( g_m_estimate, time );
        
        e1 = (data.y_LTI - G_M_estimate * data.u);
        e2 = g_m_estimate(2:end) - g_m_estimate(1:end-1);
        e3 = g_m_estimate;
        
        cost = e1'*Q_ident*e1 + e2'*S_delta*e2 + e3'*G_abs*e3;
        
    end

end

