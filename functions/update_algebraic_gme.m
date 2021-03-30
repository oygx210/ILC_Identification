function g_me_update = update_algebraic_gme(sys, data, time)
%% settings

% weight matrix
W = 1.0 * eye(time.nt);

% minimal condition regarding inverse
cond_min   = 10^-8;

% regularization candidates
reg_cand   = 10.^(-20:20);


%% algorithm

% copy data
e_prev     = data.e_prev;
e_current  = data.e_current;
G_M        = sys.G_M;

% transform error to toeplitz matrix
E       = toeplitz(e_prev, [e_prev(1) zeros(1,time.nt-1)]);

% create an identity matrix
I = eye(time.nt);    

% find a sufficient regularization
for k1=1:length(reg_cand)
    
    M = (E'*W*E + reg_cand(k1)*I);
    
    % break if condition is sufficient
    if rcond(M) > cond_min
        break
    end
end

% calculate impulse response
g_me_update = M^-1*W*E'*G_M*e_current;


end