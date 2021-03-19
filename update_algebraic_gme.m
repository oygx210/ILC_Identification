function g_me_update = update_algebraic_gme(sys, data, time)

r = data.r;
y = data.y_LTI;
e_prev= data.e_prev;
e_current= data.e_current;
G_M = sys.G_M;

W = eye(time.nt);
I = eye(time.nt);            
tik_cand = 10.^(-20:20);           
            
E=toeplitz(e_prev, [e_prev(1) zeros(1,time.nt-1)]);
e_tik=inf(1,length(tik_cand));

for k1=1:length(tik_cand)
    
    M = (E'*W*E+tik_cand(k1)*I);
    if rcond(M) > 10^-15
        
        e_tik(k1)=norm(G_M*(r-y)-E*(M^-1*W*E'*G_M*e_current));
    end
end

[~,ind_tik]=min(e_tik);

g_me_update = (E'*W*E + tik_cand(ind_tik)*I)^-1*W*E'*G_M*e_current;

end