% if model update via nonlinear optimisation is used
function g_me_update = update_fmincon_gme(sys, time, data)

options = optimoptions('fmincon','Display',...
    'iter-detailed','Algorithm','sqp', ...
    'MaxFunctionEvaluations', 1e5,...
    'StepTolerance', 1e-5,...
    'OptimalityTolerance', 1e-5);


g_me_update = fmincon( @costIdent, data.g_me_init , [],[],[],[],[],[], [] ,  options);


% fmincon cost function for model identification
    function cost = costIdent( g_me_estimate )
        
        % regularisation matrices
        Q_ident     = diag( linspace( 1, 1, time.nt ) );
        S_delta     = 1e-1*diag( linspace( 1, 10, time.nt-1 ) );
        G_abs       = 1e-2*diag( linspace( 0, 10, time.nt ) );
        
        % transform the input to a matrix
        G_ME_estimate = convert_from_vector_to_matrix( g_me_estimate , time);
        
        right = G_ME_estimate * data.e_prev;
        left = sys.G_M * data.e_current;
        
        e1 = (left - right);
        e2 = g_me_estimate(2:end) - g_me_estimate(1:end-1);
        e3 = g_me_estimate;
        
        cost = e1'*Q_ident*e1 + e2'*S_delta*e2 + e3'*G_abs*e3;
        
    end

end