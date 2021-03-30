function  g_me_update = update_tfest_gme(sys, time, data)
% find the impusle response using the system identification toolbox.

%% settings

np = 2;
nz = 0;

opt = tfestOptions;
opt.Regularization.Lambda = 1e-3;



%% algorithm

em_current = sys.G_M* data.e_current;

id_data = iddata( em_current, data.e_prev, time.dt );

sys_error_transfer = tfest( id_data, np, nz, opt );

g_me_update = time.dt*impulse(sys_error_transfer, time.vec);

end