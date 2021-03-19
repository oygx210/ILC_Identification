function [sys, time] = create_rand_system()


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

time.tend   = 10/min(abs(real(ILC.Sys.s.ps)));              % end of each iteration in seconds
time.nt     = 50*(floor(time.tend)+1);                          % number of discrete time steps per iteration

time.vec   = linspace(0, time.tend, time.nt);  % define a linear time vector
time.dt     = time.vec(2) - time.vec(1);          % calcutate the time step length in second



ILC.Sys.s.tf_s=tf(poly(ILC.Sys.s.ns), poly(ILC.Sys.s.ps));
ILC.Sys.s.tf_z=c2d(ILC.Sys.s.tf_s,time.dt);
[ILC.Sys.s.dss.A, ILC.Sys.s.dss.b, ILC.Sys.s.dss.c, ~ ] = tf2ss(ILC.Sys.s.tf_z.Numerator{1, 1}, ILC.Sys.s.tf_z.Denominator{1, 1});

ILC.Sys.s.g = zeros(time.nt, 1);
for k1 = 1:time.nt           
    ILC.Sys.s.g(k1) = ILC.Sys.s.dss.c*ILC.Sys.s.dss.A^(k1-1)*ILC.Sys.s.dss.b;           
end


sys.G_S = toeplitz(ILC.Sys.s.g, [ILC.Sys.s.g(1) zeros(1,time.nt-1)]);


% Modell System

while true
    
    tmp.rand.n=1+ILC.Sys.e_rel*2*(rand([1,ILC.Sys.s.n+1])-0.5);
    tmp.rand.m=1+ILC.Sys.e_rel*2*(rand([1,ILC.Sys.s.m+1])-0.5);

    ILC.Sys.m.tf_s=tf([zeros(1,ILC.Sys.s.n-ILC.Sys.s.m) tmp.rand.m].*ILC.Sys.s.tf_s.Numerator{1, 1}, tmp.rand.n.*ILC.Sys.s.tf_s.Denominator{1, 1});

    if isstable(ILC.Sys.m.tf_s)
        break
    end
    
end



ILC.Sys.m.tf_z=c2d(ILC.Sys.m.tf_s,time.dt);
[ILC.Sys.m.dss.A, ILC.Sys.m.dss.b, ILC.Sys.m.dss.c, ~ ] = tf2ss(ILC.Sys.m.tf_z.Numerator{1, 1}, ILC.Sys.m.tf_z.Denominator{1, 1});

ILC.Sys.m.g = zeros(time.nt, 1);
for k1 = 1:time.nt           
    ILC.Sys.m.g(k1) = ILC.Sys.m.dss.c*ILC.Sys.m.dss.A^(k1-1)*ILC.Sys.m.dss.b;           
end



sys.G_M = toeplitz(ILC.Sys.m.g, [ILC.Sys.m.g(1) zeros(1,time.nt-1)]);




end