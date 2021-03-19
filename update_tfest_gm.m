function g_m_update = update_tfest_gm(sys, data, time)
                                
                ygs = 1/data.gamma*sys.G_M *(data.e_prev - data.e_current);
                ugs = data.e_prev;
                                
                id_data = iddata( ygs, ugs, time.dt );
                
                np = 2;
                nz = 0;
                
                opt = tfestOptions;
                opt.Regularization.Lambda = 1e-8;
                
                sys_m = tfest( id_data, np, nz, opt );
                
                g_m_update = time.dt*impulse(sys_m, time.vec);
                
end