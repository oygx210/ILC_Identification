function [sys, time] = create_single_system




%% time parameters
% number of discrete time steps per iteration
time.nt = 100;
% start of each iteration in seconds
tstart = 0;
% end of each iteration in seconds
time.tend = 2;

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


% define a PT2 system with the "real" parameters
D = 0.5;
w0 = 3;
K = 0.5;

% % define a second PT1 system with model errors, "guessed" parameters
Km = 1.8;
Tm = 0.5;





%% model definitions

% "real" transfer function definition
num = K*w0^2;
den = [1, 2*D*w0, w0^2];

% convert "real" model to cont-time and discrete-time state-space
[A,B,C,~] = tf2ss(num, den);
AD = expm(A*time.dt);
BD = (AD - eye(size(AD)))*A^-1*B;
CD = C;

% convert "real" model to teoplitz matrix
sys.G_S = convert_state_space_to_matrix(AD, BD, CD, time);


% "guessed" transfer function definition
% num_m = Km*wm0^2;
% den_m = [1, 2*Dm*wm0, wm0^2];
num_m = Km;
den_m = [1/Tm, 1];

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