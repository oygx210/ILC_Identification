function g_m_update = update_tfest_gm_io(data, time)
                                
                ygs = data.y_LTI;
                ugs = data.u;
                
                id_data = iddata( ygs, ugs, time.dt );
                
                np = 2;
                nz = 0;
                
                opt = tfestOptions;
                opt.Regularization.Lambda = 1e-8;
                
                sys_m = tfest( id_data, np, nz, opt );
                
                g_m_update = time.dt*impulse(sys_m, time.vec);
                
end