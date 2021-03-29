function  g_me_update = update_tfest_gme(sys, time, data)


em_current = sys.G_M* data.e_current;

id_data = iddata( em_current, data.e_prev, time.dt );

% if the pure error model G_E is estimated, there should be a
% feed-thourgh and therefore np = nz. Since G_ME is estimated,
% a feed-through is not expected.
np = 2;
nz = 0;

opt = tfestOptions;
opt.Regularization.Lambda = 1e-3;

sys_error_transfer = tfest( id_data, np, nz, opt );

g_me_update = time.dt*impulse(sys_error_transfer, time.vec);

end