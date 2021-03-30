function g_m_update = update_tfest_gm(sys, data, time)
% find the impusle response using the system identification toolbox.

%% settings

np = 2;
nz = 0;

opt = tfestOptions;
opt.Regularization.Lambda = 1e-3;



%% algorithm

% copy data
G_M         = sys.G_M;
e_current   = data.e_current;
e_prev      = data.e_prev;
gamma       = data.gamma;

                ygs = 1/gamma *G_M *(e_prev - e_current);
                ugs = e_prev;
                
id_data = iddata( ygs, ugs, time.dt );

sys_m = tfest( id_data, np, nz, opt );

g_m_update = time.dt*impulse(sys_m, time.vec);

end