% fmincon cost function for input trajectory optimisation
function cost = cost_control( u_use , sys, time, r, y_LTI)

nt = time.nt;

Q_abs       = eye(nt,nt);
Q_delta_u   = 1e-2*eye(nt-1,nt-1);

y_sim       = sys.G_M * u_use;

e1          = r - y_LTI - y_sim;
e2          = u_use(2:end) - u_use(1:end-1);

cost = e1'*Q_abs*e1 + e2'*Q_delta_u*e2;

end