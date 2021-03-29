function set = no_update_2( k_run )


% load default settings and only override changes
set = default();



% enable model update
if k_run == 1
    set.enable_error_update = true;
else
    set.enable_error_update = false;
end

end