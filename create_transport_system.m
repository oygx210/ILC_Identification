function [sys, time] = create_transport_system




%% time parameters
% number of discrete time steps per iteration
time.nt = 1000;
% start of each iteration in seconds
tstart = 0;
% end of each iteration in seconds
time.tend = 15;

% define a linear time vector
time.vec = linspace(tstart, time.tend, time.nt);
% calcutate the time step length in seconds
time.dt = time.vec(2) - time.vec(1);




%% model parameters
% We apply two models. The first one is defined as "real" or system model,
% which is used for estimating the next cycle and which simulates the real
% model. The second one is defined as "guessed" model, subscript m, which
% is used for the update of the input trajectory. It can be viewed as an
% initial guess of the real model. The goal of this algorithm is to match
% the "guessed" parameters to the "real" parameters.


P = 100;
D = 20;
mass = 12;
fv = 10;

mass_m = 0.5*mass;
fv_m = 2*fv;


%% model definitions

% "real" transfer function definition
num = [D, P];
den = [mass, (D + fv), P];

% convert "real" model to cont-time and discrete-time state-space
[A,B,C,~] = tf2ss(num, den);
AD = expm(A*time.dt);
BD = (AD - eye(size(AD)))*A^-1*B;
CD = C;


% convert "real" model to teoplitz matrix
sys.G_S = convert_state_space_to_matrix(AD, BD, CD, time);


% "guessed" transfer function definition
num_m = [D, P];
den_m = [mass_m, (D + fv_m), P];


% convert "guessed" model to cont-time and discrete-time state-space
[Am,Bm,Cm,~] = tf2ss(num_m, den_m);
ADm = expm(Am*time.dt);
BDm = (ADm - eye(size(ADm)))*Am^-1*Bm;
CDm = Cm;

% convert "guessed" model to ss variable
% ss_model = ss( Am,Bm,Cm,0 );

% convert "guessed" model to teoplitz matrix
sys.G_M = convert_state_space_to_matrix(ADm, BDm, CDm, time);



end