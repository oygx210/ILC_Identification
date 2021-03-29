function set = noise_sim_50( k1 )



% load default settings and only override changes
set = default();


% type of model update
if k1 <= 25
    set.model_update_type = 'tfest_gm_io';
else
    set.model_update_type = 'tfest_gm';
end


% noise settings
set.enable_noise = true;
set.mu = 0.5  *2*(rand-0.5);
set.sigma = 0.01  *2*(rand-0.5);



end