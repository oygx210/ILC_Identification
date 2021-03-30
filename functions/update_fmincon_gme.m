% if model update via nonlinear optimisation is used
function g_me_update = update_fmincon_gme(sys, time, data)
%% settings
options = optimoptions('fmincon','Display',...
    'iter-detailed','Algorithm','sqp', ...
    'MaxFunctionEvaluations', 1e5,...
    'StepTolerance', 1e-5,...
    'OptimalityTolerance', 1e-5);

% define regularisation matrices
Q_ident     = diag( linspace( 1, 1, time.nt ) );
S_delta     = 1e-1*diag( linspace( 1, 10, time.nt-1 ) );
G_abs       = 1e-2*diag( linspace( 0, 10, time.nt ) );


%% algorithm

% perform the optimization
g_me_update = fmincon( @costIdent, data.g_me_init , [],[],[],[],[],[], [] ,  options);


    % fmincon cost function for model identification
    function cost = costIdent( g_me_estimate )
        
        % transform the input to a matrix
        G_ME_estimate = convert_from_vector_to_matrix( g_me_estimate , time);
        
        % calculate right and lefthand side of equation
        right       = G_ME_estimate * data.e_prev;
        left        = sys.G_M * data.e_current;
        error_est   = (left - right);
        
        % calculate impulse change rate for regularization
        g_delta_reg = g_me_estimate(2:end) - g_me_estimate(1:end-1);
        
        % calculate absolute impulse for regularization
        g_reg       = g_me_estimate;
        
        % calculate the cost
        cost = error_est'*Q_ident*error_est + g_delta_reg'*S_delta*g_delta_reg + g_reg'*G_abs*g_reg;
        
    end

end