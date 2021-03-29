function [tik_reg, tik_reg_error, tik_ok] = tikhonov_reg(sys, time, r, u, y, tik_local_search, tik_reg)

G_M        = sys.G_M;
I          = eye(time.nt);



if tik_local_search
    
    tik_cand = [tik_reg, 0.9 * tik_reg];
    
else
    
    tik_cand   = 10.^(-20:20);
end

tik_cand   = 10.^(-20:20);

gm1 = 1;
e_k = r - y;

A = (G_M/gm1*( r - G_M*u ))' * e_k;

tik_reg_alg = 2*(e_k' * e_k)/A - gm1^2


tik_err = inf(1,length(tik_cand));
tik_ok  = false;

for k1=1:length(tik_cand)
    
    % M = (E'*W*E+tik_cand(k1)*I);
    M = (G_M' * G_M + tik_cand(k1)*I);
    
    if rcond(M) > 10^-15
        
        tik_ok      = true;
        
              
        tik_err(k1) = norm(r - G_M*u - M^-1 * G_M' * (r - y));
        
        % tik_err(k1) = norm(G_M*(r-y) - E*(M^-1*W*E'*G_M*e_current));
    end
end

hold on
semilogx(tik_err)

if tik_ok
    
    % take smallest candidate
    [tik_reg_error, ind_tik]=min(tik_err);
    
    tik_reg = tik_cand(ind_tik);
    
else
    
    warning('No tikhonov candidate found.')
    tik_reg = 1;
    tik_reg_error = inf;
end

end