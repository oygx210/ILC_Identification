function g_me_update = update_algebraic_gme(sys, data, time)


% tik_reg    = data.tik_reg;
e_prev     = data.e_prev;
e_current  = data.e_current;
G_M        = sys.G_M;


E       = toeplitz(e_prev, [e_prev(1) zeros(1,time.nt-1)]);


W = eye(time.nt);
I = eye(time.nt);                  
  

tik_cand   = 10.^(-20:20);

for k1=1:length(tik_cand)
    
    M = (E'*W*E + tik_cand(k1)*I);
    
    if rcond(M) > 10^-8
       
        % tik = tik_cand(k1);
        
        break
    end
end


g_me_update = M^-1*W*E'*G_M*e_current;


end