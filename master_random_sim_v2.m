function out = master_random_sim_v2( varargin )

addpath('inp')


if isempty(varargin)

    settings.enable_rand_model = false;
    
    

%% time parameters
% number of discrete time steps per iteration
nt = 200;
% start of each iteration in seconds
tstart = 0;
% end of each iteration in seconds
tend = 2;

% define a linear time vector
time = linspace(tstart, tend, nt);
% calcutate the time step length in seconds
dt = time(2) - time(1);



%% ILC
% number of iterations
nz = 10;

% maximum of model updates
gUpdateMax = 2;

% first cycle of the model update
gUpdateCycle = 2;

% learning rate
gamma = 0.5;

% number of cases to compare
nc = 3;





%% reference trajectory parameters
% define a reference trajectory, use a step function
% r = zeros(nt,1);
% r(80:end) = 10;
s = tf('s');


Range = [-1,1];
Band = [0 0.05];
prbs = idinput(2540,'prbs',Band,Range);
rand_traj = 0.1;
short_prbs = prbs( (1:nt) + floor( (length(prbs) - nt)*rand_traj ));
Tf = 0.05;
filter = tf(1,[Tf, 1])^7;

% 
% ILC.r         = idinput(nt,'rbs',[0 0.05],[-1,1]);
% 
% ILC.r_filter  = tf(1,[0.05, 1])^7;
% ILC.r         = lsim(ILC.r_filter, ILC.r, time);
% 





% [A,B,C,D] = tf2ss( filter.num{1}, filter.den{1} );
% filter_ss = ss(A,B,C,D);
% x0 = zeros(7,1);
% r = lsim(filter_ss, input, contTimeVecSim, x0)';
r = lsim(filter, short_prbs, time);


% accel = zeros(1,nt);
% accel(floor(0.1*nt):floor(0.25*nt)) = 1;
% accel(ceil(0.65*nt):end-floor(0.2*nt)) = -1;
% 
% r = lsim(1/s^2, accel, time);




%% algorithm parameters

% Enable model update algorithm.
enable_error_update = true;
% Enable model update each iteration, for "gUpdateMax" updates max.
% If disabled, model update is performed only once.
enable_error_update_every_iter = false;


% enable model update using the "real" system.
% this is only implemented for validation and cannot be applied to a real
% process.
enable_update_real_system = false;
% or enable model update using the MATLAB system identification toolbox
enable_update_tfest_gme = false;
% or enable model update using a constrained optimisation (recommended)
enable_update_fmincon = false;

enable_update_conv_fmincon = false;

enable_update_GS = true;

% update the input trajectory using quadprog
enable_control_quadprog = true;

% or update the input trajectory using fmincon
enable_control_fmincon = false;

enable_control_inverse = false;



enable_noise_y = true;
enable_noise_u = true;




end
    

if settings.enable_rand_model

% model
ILC.Sys.e_rel = 0.8;

% real System

ILC.Sys.s.n = randi([1 5]);
ILC.Sys.s.m = randi([0 ILC.Sys.s.n-1]);
ILC.Sys.s.n_ips = randi([0 floor(ILC.Sys.s.n/2)]);

ILC.Sys.s.ps = -1-9*rand(1,ILC.Sys.s.n-2*ILC.Sys.s.n_ips);
tmp.Re= -1-9*rand(1,ILC.Sys.s.n_ips);
tmp.Im= 10*rand(1,ILC.Sys.s.n_ips);
ILC.Sys.s.ps =[ILC.Sys.s.ps, tmp.Re+tmp.Im*1i, tmp.Re-tmp.Im*1i];
ILC.Sys.s.ns =-1-9*rand(1,ILC.Sys.s.m);


% time 

t.end   = 10/min(abs(real(ILC.Sys.s.ps)));              % end of each iteration in seconds
t.n     = 50*(floor(t.end)+1);                          % number of discrete time steps per iteration

t.vec   = linspace(0, t.end, t.n);  % define a linear time vector
t.d     = t.vec(2) - t.vec(1);          % calcutate the time step length in second



ILC.Sys.s.tf_s=tf(poly(ILC.Sys.s.ns), poly(ILC.Sys.s.ps));
ILC.Sys.s.tf_z=c2d(ILC.Sys.s.tf_s,t.d);
[ILC.Sys.s.dss.A, ILC.Sys.s.dss.b, ILC.Sys.s.dss.c, ~ ] = tf2ss(ILC.Sys.s.tf_z.Numerator{1, 1}, ILC.Sys.s.tf_z.Denominator{1, 1});

ILC.Sys.s.g = zeros(t.n, 1);
for k1 = 1:t.n           
    ILC.Sys.s.g(k1) = ILC.Sys.s.dss.c*ILC.Sys.s.dss.A^(k1-1)*ILC.Sys.s.dss.b;           
end
ILC.G_s=toeplitz(ILC.Sys.s.g, [ILC.Sys.s.g(1) zeros(1,t.n-1)]);


% Modell System

while true
    
    tmp.rand.n=1+ILC.Sys.e_rel*2*(rand([1,ILC.Sys.s.n+1])-0.5);
    tmp.rand.m=1+ILC.Sys.e_rel*2*(rand([1,ILC.Sys.s.m+1])-0.5);

    ILC.Sys.m.tf_s=tf([zeros(1,ILC.Sys.s.n-ILC.Sys.s.m) tmp.rand.m].*ILC.Sys.s.tf_s.Numerator{1, 1}, tmp.rand.n.*ILC.Sys.s.tf_s.Denominator{1, 1});

    if isstable(ILC.Sys.m.tf_s)
        break
    end
    
end

end






ILC.Sys.m.tf_z=c2d(ILC.Sys.m.tf_s,t.d);
[ILC.Sys.m.dss.A, ILC.Sys.m.dss.b, ILC.Sys.m.dss.c, ~ ] = tf2ss(ILC.Sys.m.tf_z.Numerator{1, 1}, ILC.Sys.m.tf_z.Denominator{1, 1});

ILC.Sys.m.g = zeros(t.n, 1);
for k1 = 1:t.n           
    ILC.Sys.m.g(k1) = ILC.Sys.m.dss.c*ILC.Sys.m.dss.A^(k1-1)*ILC.Sys.m.dss.b;           
end
ILC.G_m=toeplitz(ILC.Sys.m.g, [ILC.Sys.m.g(1) zeros(1,t.n-1)]);








%% Parameter 

% ILC
ILC.algebraic   = true;
ILC.n           = 10;                    % number of iterations

% opti

ILC.W=eye(t.n);%diag(flip(1./(1:t.n)));
para.fmincon_options = optimoptions('fmincon','Display','off','Algorithm','sqp', 'MaxFunctionEvaluations', 1e5);


%% reference trajectory 

ILC.r     = idinput(t.n,'rbs',[0 0.05],[-1,1]);

ILC.r_filter  = tf(1,[0.05, 1])^7;
ILC.r     = lsim(ILC.r_filter, ILC.r, t.vec);


function out = master_NOILC_fmincon_v6

% Jonas Weigand <jonas.weigand@mv.uni-kl.de>
% Gajanan Kanagalingam <gajanan.kanagalingam@gmail.com>
% 18-09-2020



% close all
clc



%% time parameters
% number of discrete time steps per iteration
nt = 200;
% start of each iteration in seconds
tstart = 0;
% end of each iteration in seconds
tend = 2;

% define a linear time vector
time = linspace(tstart, tend, nt);
% calcutate the time step length in seconds
dt = time(2) - time(1);



%% ILC
% number of iterations
nz = 10;

% maximum of model updates
gUpdateMax = 2;

% first cycle of the model update
gUpdateCycle = 2;

% learning rate
gamma = 0.5;

% number of cases to compare
nc = 3;





%% reference trajectory parameters
% define a reference trajectory, use a step function
% r = zeros(nt,1);
% r(80:end) = 10;
s = tf('s');


Range = [-1,1];
Band = [0 0.05];
prbs = idinput(2540,'prbs',Band,Range);
rand_traj = 0.1;
short_prbs = prbs( (1:nt) + floor( (length(prbs) - nt)*rand_traj ));
Tf = 0.05;
filter = tf(1,[Tf, 1])^7;

% 
% ILC.r         = idinput(nt,'rbs',[0 0.05],[-1,1]);
% 
% ILC.r_filter  = tf(1,[0.05, 1])^7;
% ILC.r         = lsim(ILC.r_filter, ILC.r, time);
% 





% [A,B,C,D] = tf2ss( filter.num{1}, filter.den{1} );
% filter_ss = ss(A,B,C,D);
% x0 = zeros(7,1);
% r = lsim(filter_ss, input, contTimeVecSim, x0)';
r = lsim(filter, short_prbs, time);


% accel = zeros(1,nt);
% accel(floor(0.1*nt):floor(0.25*nt)) = 1;
% accel(ceil(0.65*nt):end-floor(0.2*nt)) = -1;
% 
% r = lsim(1/s^2, accel, time);




%% algorithm parameters

% Enable model update algorithm.
enable_error_update = true;
% Enable model update each iteration, for "gUpdateMax" updates max.
% If disabled, model update is performed only once.
enable_error_update_every_iter = false;


% enable model update using the "real" system.
% this is only implemented for validation and cannot be applied to a real
% process.
enable_update_real_system = false;
% or enable model update using the MATLAB system identification toolbox
enable_update_tfest_gme = false;
% or enable model update using a constrained optimisation (recommended)
enable_update_fmincon = false;

enable_update_conv_fmincon = false;

enable_update_GS = true;

% update the input trajectory using quadprog
enable_control_quadprog = true;

% or update the input trajectory using fmincon
enable_control_fmincon = false;

enable_control_inverse = false;



enable_noise_y = true;
enable_noise_u = true;



%% model parameters
% We apply two models. The first one is defined as "real" or system model,
% which is used for estimating the next cycle and which simulates the real
% model. The second one is defined as "guessed" model, subscript m, which
% is used for the update of the input trajectory. It can be viewed as an
% initial guess of the real model. The goal of this algorithm is to match
% the "guessed" parameters to the "real" parameters.


% define a PT2 system with the "real" parameters
D = 0.5;
w0 = 3;
K = 1;

% % define a second PT2 system with model errors, "guessed" parameters
Km = 1.5 * K;
wm0 = 2 * w0;
Dm = 1.5 * D;
Tm = 0.5;





%% model definitions

% "real" transfer function definition
num = K*w0^2;
den = [1, 2*D*w0, w0^2];

% convert "real" model to cont-time and discrete-time state-space
[A,B,C,~] = tf2ss(num, den);
AD = expm(A*dt);
BD = (AD - eye(size(AD)))*A^-1*B;
CD = C;

% convert "real" model to ss variable
sys_real = ss( A,B,C,0 );

% convert "real" model to teoplitz matrix
G_S = convert_state_space_to_matrix(AD, BD, CD);


% "guessed" transfer function definition
% num_m = Km*wm0^2;
% den_m = [1, 2*Dm*wm0, wm0^2];
num_m = Km;
den_m = [1/Tm, 1];

% convert "guessed" model to cont-time and discrete-time state-space
[Am,Bm,Cm,~] = tf2ss(num_m, den_m);
ADm = expm(Am*dt);
BDm = (ADm - eye(size(ADm)))*Am^-1*Bm;
CDm = Cm;

% convert "guessed" model to ss variable
sys_model = ss( Am,Bm,Cm,0 );

% convert "guessed" model to teoplitz matrix
G_M_init = convert_state_space_to_matrix(ADm, BDm, CDm);



%% calculate and plot the model error without changes of the model

g_me_update_matrix = eye(nt) - G_S*G_M_init^-1;
previous_model_error = norm( g_me_update_matrix ) %#ok



%% initialise variables

cycle_error_each_iter   = zeros(nz, nc + 1); % all cases plus perfect measure
new_model_error         = zeros(nz, 1);
yIterAll                = zeros(nz, nt, nc); % all cycles, all timesteps, with and without model update

g_m_update              = zeros(nt,1);
g_me_update             = zeros(nt,1);
gUpdateCounter          = 0;


g_m_all = zeros(nt, 3); % initial, real update
g_m_all(:, 1) = G_M_init(:,1);
g_m_all(:, 2) = G_S(:,1);

%% Core Algorithm

% with and without model update
for k1 = 1:nc
    
    G_M = G_M_init;
    
    % the initial input equals the reference trajectory
    u = r;
    
    g_m_update = G_M_init(:,1);
    
    % for all cycles
    for k2 = 1:nz
        
        if enable_noise_u
            u = u + 0.01*randn(size(u));
            
            noise_norm = norm( 0.01*randn(size(u)) );
        end
        
        % calculate the "real" system
        y_LTI = G_S * u;
        
        if enable_noise_y
            
            y_LTI = y_LTI + 0.1 + 0.01*randn(size(y_LTI));
            
            noise_norm = norm (0.01*randn(size(y_LTI)));
        end
        
        % save statistics
        cycle_error_each_iter(k2,k1) = norm(r - y_LTI);
        yIterAll(k2, :, k1) = y_LTI;
        
        
        if k2 > 1
            err_prev = err_curr;
        else
            err_prev = 0*r;
        end
        err_curr = r - y_LTI;
        
        
        
        if k2 <= gUpdateCycle
            cycle_error_each_iter(k2,3) = norm(r - y_LTI);
            
        else
            cycle_error_each_iter(k2,3) = (1 - gamma)*cycle_error_each_iter(k2-1,3);
        end
        
        
        % model update
        if k1 == 1 && enable_error_update && (k2 == gUpdateCycle || enable_error_update_every_iter) && (gUpdateCounter < gUpdateMax) && k2 >= gUpdateCycle
            
            gUpdateCounter = gUpdateCounter+1;
            
            % get the current and previous error
            e_current   = r - yIterAll(k2,:,k1)';
            e_prev      = r - yIterAll(k2-1,:,k1)';
            
            
            
            % if information from the real system is used
            % (this cannot be applied in an online algorithm - use only for validation)
            if enable_update_real_system
                
                % calculate the inverse
                G_M_inv = G_M^-1;
                
                G_E = (eye(nt) - G_S*G_M_inv);
                g_me_update_matrix = G_M*G_E;
                g_me_update = g_me_update_matrix(:,1);
                
                figure
                hold on
                bodeplot(sys_real, 'b')
                bodeplot(sys_model, 'r')
                bodeplot(sys_error_real, 'm')
                
                sys_error_real = 1 - sys_real*sys_model^-1;
                
                figure
                hold on
                pzplot(sys_error_real, 'm')
                
                
            end
            
            % if model update via transfer function estimation is used
            if enable_update_tfest_gme
                
                
                em_current = G_M* e_current;
                
                data = iddata( em_current, e_prev, dt );
                
                % if the pure error model G_E is estimated, there should be a
                % feed-thourgh and therefore np = nz. Since G_ME is estimated,
                % a feed-through is not expected.
                np = 2;
                nz = 0;
                
                opt = tfestOptions;
                opt.Regularization.Lambda = 1e-3;
                
                sys_error_transfer = tfest( data, np, nz, opt );
                
                g_me_update = dt*impulse(sys_error_transfer, time);
                
                
                figure
                hold on
                bodeplot(sys_real, 'b')
                bodeplot(sys_model, 'r')
                bodeplot(sys_error_transfer, 'g')
                % bodeplot(sys_error_real, 'm')
                
                figure
                hold on
                pzplot(sys_error_transfer, 'g')
                % pzplot(sys_error_real, 'm')
                
                
            end

            % if model update via nonlinear optimisation is used
            if enable_update_fmincon
                
                options = optimoptions('fmincon','Display',...
                    'iter-detailed','Algorithm','sqp', ...
                    'MaxFunctionEvaluations', 1e5,...
                    'StepTolerance', 1e-5,...
                    'OptimalityTolerance', 1e-5);
                
                
                g_me_update = fmincon( @costIdent, g_me_update , [],[],[],[],[],[], [] ,  options);
                
            end
            
            % if model update via nonlinear optimisation is used
            if enable_update_conv_fmincon
                
                options = optimoptions('fmincon','Display',...
                    'iter-detailed','Algorithm','sqp', ...
                    'MaxFunctionEvaluations', 1e5,...
                    'StepTolerance', 1e-5,...
                    'OptimalityTolerance', 1e-5);
                
                
                g_m_update = fmincon( @costIdentConv, g_m_update , [],[],[],[],[],[], [] ,  options);
                
                % transfer update to a matrix
                G_M_update = convert_from_vector_to_matrix( g_m_update );
                
                % weight update
                delta = 0.1;
                
                % perform the update
                G_M =  delta*G_M + (1 - delta)*G_M_update;
                
            end
            
            % if model update via nonlinear optimisation is used
            if enable_update_GS
                
%                 options = optimoptions('fmincon','Display',...
%                     'iter-detailed','Algorithm','sqp', ...
%                     'MaxFunctionEvaluations', 1e5,...
%                     'StepTolerance', 1e-5,...
%                     'OptimalityTolerance', 1e-5);
%                 
%                 
%                 g_m_update = fmincon( @costIdentGS, g_m_update , [],[],[],[],[],[], [] ,  options);
%                 
%                 


% G_M *(err_curr - err_prev) + gamma*G_M_estimate*err_prev
% 1/gamma*G_M *(err_prev - err_curr) = G_M_estimate*err_prev

                ygs = 1/gamma*G_M *(err_prev - err_curr);
                ugs = err_prev;

                % ygs = (1/gamma*G_M *(err_prev - err_curr) + y_LTI)/2;
                % ugs = (err_prev + u)/2;
                    
                  ygs = y_LTI;
                  ugs = u;
                
                data = iddata( ygs, ugs, dt );

                % data = iddata( y_LTI, u, dt );
                
                % if the pure error model G_E is estimated, there should be a
                % feed-thourgh and therefore np = nz. Since G_ME is estimated,
                % a feed-through is not expected.
                np = 4;
                nze = 2;
                
                opt = tfestOptions;
                opt.Regularization.Lambda = 1e-8;
                
                sys_m = tfest( data, np, nze, opt );
                
                g_m_update = dt*impulse(sys_m, time);
                
                
                % transfer update to a matrix
                G_M_update = convert_from_vector_to_matrix( g_m_update );
                
                % weight update
                delta = 0;
                
                % perform the update
                G_M =  delta*G_M + (1 - delta)*G_M_update;
                
            else
                
                % transfer update to a matrix
                G_ME = convert_from_vector_to_matrix( g_me_update );
                
                % save the norm of the model update
                new_model_error(k2) = norm(G_ME, 2);
                
                % perform the update
                G_M =  1/gamma*(G_M - G_ME);
                
            end
            
            g_m_all(:, 3) = G_M(:,1);
            
        end
        
        % for each but the last cycle apply the update of the input trajectory
        if k2 < nz
            
            % calulate the update using a quadratic form
            if enable_control_quadprog
                
                options = optimset('Display','off','Algorithm','interior-point-convex');
                
                % weighting matrix for quadratic optimization
                Q = eye(nt,nt);
                
                f = -2*(r - y_LTI)'*Q*G_M;
                H = 2*G_M'*Q*G_M;
                u_update = quadprog(H,f, [],[],[],[],[],[],[], options);
                
                % or calculate the update using a nonlinear optimisation
            elseif enable_control_fmincon
                
                options = optimoptions('fmincon','Display',...
                    'off','Algorithm','sqp', ...
                    'MaxFunctionEvaluations', 1e5);
                
                u_update = fmincon( @costControl, u , [],[],[],[],[],[],[],  options);
                
                
            elseif enable_control_inverse
                
                u_update = G_M^-1 * (r - y_LTI);
                
            else
                
                u_update = 0;
            end
            
            % apply the update to the trajectory
            u = u + gamma * u_update;
            
        end
    end
end


%% Visualisation

figure
subplot(3,1,1)
hold on
plot(time, yIterAll(1,:,1),'b');
plot(time, yIterAll(2,:,1),'m');
plot(time, yIterAll(end,:,1),'k');
% plot(time, yIterAll(1,:,2),'b');
% plot(time, yIterAll(2,:,2),'m');
% plot(time, yIterAll(end,:,2),'k');
plot(time, r,'r-.');
legend('1st Cycle','2nd Cycle', 'Final Cycle', 'Reference', 'Location','eastoutside')
xlabel('Time')
ylabel('Output')
subplot(3,1,2)
hold on
plot(time, r - yIterAll(1,:,1)','b');
plot(time, r - yIterAll(2,:,1)','m');
plot(time, r - yIterAll(end,:,1)','k');
% plot(time, r - yIterAll(1,:,2)','b');
% plot(time, r - yIterAll(2,:,2)','m');
% plot(time, r - yIterAll(end,:,2)','k');
legend('1st Error','2nd Error', 'Final Error', 'Location','eastoutside')
xlabel('Time')
ylabel('Error')
subplot(3,1,3)
hold on
plot(time, g_m_all(:,2),'m');
plot(time, g_m_all(:,1),'b');
plot(time, g_m_all(:,3),'k-.');
legend('True', 'Initial', 'Updated', 'Location','eastoutside')
xlabel('Time')
ylabel('Impuse Response')

figure
grid on
hold on
plot(time, g_m_all(:,2),'m');
plot(time, g_m_all(:,1),'b');
plot(time, g_m_all(:,3),'k-.');
legend('True', 'Initial', 'Updated', 'Location','eastoutside')
xlabel('Time')
ylabel('Impuse Response')

out.g_m_all = g_m_all;
out.time = time;

figure
hold on
grid on
plot(cycle_error_each_iter(:,1), 'r-*')
plot(cycle_error_each_iter(:,2), 'k-*')
plot(cycle_error_each_iter(:,3), 'mo')
ax = gca;
ax.YScale = 'log';
legend('Model Update', 'Model Init', 'Perfect', 'Location','southwest')
xlabel('Cycle')
ylabel('Error in L2-Norm')










%% nested functions

% fmincon cost function for conventional model identification
    function cost = costIdentGS( g_m_estimate )
        
        % regularisation matrices
        Q_ident     = diag( linspace( 1, 1, nt ) );
        S_delta     = 1e1*diag( linspace( 1, 10, nt-1 ) );
        G_abs       = 1e0*diag( linspace( 0, 10, nt ) );
        
        % transform the input to a matrix
        G_M_estimate = convert_from_vector_to_matrix( g_m_estimate );
        
        e1 = (y_LTI - G_M_estimate * u);
        e2 = g_m_estimate(2:end) - g_m_estimate(1:end-1);
        e3 = g_m_estimate;
        % e4 = (err_curr - (eye(nt) - gamma*G_M_estimate*G_M^-1)*err_prev);
        e4 = 0*(G_M *(err_curr - err_prev) + gamma*G_M_estimate*err_prev);
        
        cost = e1'*Q_ident*e1 + e2'*S_delta*e2 + e3'*G_abs*e3 + e4'*Q_ident*e4;
        
    end

% fmincon cost function for conventional model identification
    function cost = costIdentConv( g_m_estimate )
        
        % regularisation matrices
        Q_ident     = diag( linspace( 1, 1, nt ) );
        S_delta     = 1e-1*eye(nt -1);
        G_abs       = 1e-2*diag( linspace( 1, 1, nt ) );
        
        % transform the input to a matrix
        G_M_estimate = convert_from_vector_to_matrix( g_m_estimate );
        
        e1 = (y_LTI - G_M_estimate * u);
        e2 = g_m_estimate(2:end) - g_m_estimate(1:end-1);
        e3 = g_m_estimate;
        
        cost = e1'*Q_ident*e1 + e2'*S_delta*e2 + e3'*G_abs*e3;
        
    end

% fmincon cost function for model identification
    function cost = costIdent( g_me_estimate )
        
        % regularisation matrices
        Q_ident     = diag( linspace( 1, 1, nt ) );
        S_delta     = 1e-1*diag( linspace( 1, 10, nt-1 ) );
        G_abs       = 1e-2*diag( linspace( 0, 10, nt ) );
        
        % transform the input to a matrix
        G_ME_estimate = convert_from_vector_to_matrix( g_me_estimate );
        
        right = G_ME_estimate * e_prev;
        left = G_M * e_current;
        
        e1 = (left - right);
        e2 = g_me_estimate(2:end) - g_me_estimate(1:end-1);
        e3 = g_me_estimate;
        
        cost = e1'*Q_ident*e1 + e2'*S_delta*e2 + e3'*G_abs*e3;
        
    end

% fmincon cost function for input trajectory optimisation
    function cost = costControl( u_use )
        
        Q_abs       = eye(nt,nt);
        Q_delta_u   = 1e-2*eye(nt-1,nt-1);
        
        y_sim       = G_M * u_use;
        
        e1          = r - y_LTI - y_sim;
        e2          = u_use(2:end) - u_use(1:end-1);
        
        cost = e1'*Q_abs*e1 + e2'*Q_delta_u*e2;
    end

% generate a toeplitz matrix model from a discrete-time state space
% model
    function F = convert_state_space_to_matrix(AD_in, BD_in, CD_in)
        
        f_vec = zeros(nt, 1);
        
        for k3 = 1:nt
            
            f_vec(k3) = CD_in*AD_in^(k3-1)*BD_in;
            
        end
        
        F = convert_from_vector_to_matrix( f_vec );
    end

% generate a toeplitz matrix from vector
    function H = convert_from_vector_to_matrix( h_vec )
        
        first_row_h = [h_vec(1), zeros(1,nt-1)];
        first_colum_h = h_vec;
        
        H = toeplitz(first_colum_h, first_row_h);
        
    end


end

end