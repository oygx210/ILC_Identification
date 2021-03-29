function MASTER_multirun


addpath('functions', 'settings')

close all
clear
clc

n_runs = 50;
all_out = cell(n_runs, 1);






for k1 = 1:n_runs
    
    disp(['Start run ',num2str(k1),' of ',num2str(n_runs),'.'])
    
    % set = default();
    % set = no_update_2(k1);
    % set = noise_sim_50(k1);
    set = rand_sim_50(k1);
    
    
    % perform simulation
    all_out{k1} = MASTER_simulation(set);
end











%% PLOT

figure
for k1 = 1:n_runs
    hold on
    plot(all_out{k1}.time.vec, all_out{k1}.yIterAll(end,:),'k');
    plot(all_out{k1}.time.vec, all_out{k1}.r,'r-.');
    legend('Final Cycle', 'Reference', 'Location','eastoutside')
    xlabel('Time')
    ylabel('Output')
end


figure
for k1 = 1:n_runs
    subplot(3,1,1)
    hold on
    plot(all_out{k1}.time.vec, all_out{k1}.yIterAll(1,:),'b');
    plot(all_out{k1}.time.vec, all_out{k1}.yIterAll(2,:),'m');
    plot(all_out{k1}.time.vec, all_out{k1}.yIterAll(end,:),'k');
    plot(all_out{k1}.time.vec, all_out{k1}.r,'r-.');
    legend('1st Cycle','2nd Cycle','Final Cycle', 'Reference', 'Location','eastoutside')
    % xlabel('Time')
    ylabel('Output')
    subplot(3,1,2)
    hold on
    plot(all_out{k1}.time.vec, all_out{k1}.r - all_out{k1}.yIterAll(1,:)','b');
    plot(all_out{k1}.time.vec, all_out{k1}.r - all_out{k1}.yIterAll(2,:)','m');
    plot(all_out{k1}.time.vec, all_out{k1}.r - all_out{k1}.yIterAll(end,:)','k');
    legend('1st Error','2nd Error', 'Final Error', 'Location','eastoutside')
    % xlabel('Time')
    ylabel('Error')
    subplot(3,1,3)
    hold on
    plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,2),'m');
    plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,1),'b');
    plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,3),'k-.');
    legend('True', 'Initial', 'Updated', 'Location','eastoutside')
    xlabel('Time')
    ylabel('Impuse Response')
end

figure    
grid on
    hold on
for k1 = 1:n_runs

    if k1 <= n_runs/2
        plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,3),'r');
    else
        plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,3),'b');
    end
end
    plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,3),'k-.');
    legend('IO', 'Con', 'True', 'Location','eastoutside')
    xlabel('Time')
    ylabel('Impuse Response')

% nt = all_out{1}.time.nt;
% g_m = zeros(n_runs, nt);
% for k1 = 1:n_runs
% 
%     g_m(k1, :) = all_out{k1}.g_m_all(:,3);
% end
%     
% n       = floor(n_runs/2);
% g_m_io  = g_m(1:n, :);
% g_m_con = g_m(n+1:end, :);
% 
% % g_m_io_p = max(g_m_io);
% % g_m_io_m = min(g_m_io);
% % g_m_con_p = max(g_m_con);
% % g_m_con_m = min(g_m_con);
% g_m_io_p = mean(g_m_io) + var(g_m_io);
% g_m_io_m = mean(g_m_io) - var(g_m_io);
% g_m_con_p = mean(g_m_con) + var(g_m_con);
% g_m_con_m = mean(g_m_con) - var(g_m_con);
% x = all_out{1}.time.vec;
% 
% figure
%     grid on
%     hold on
%     fill( [x fliplr(x)],  [g_m_io_p fliplr(g_m_io_m)], 'ro-');
%     fill( [x fliplr(x)],  [g_m_con_p fliplr(g_m_con_m)], 'bo-');
%     plot(x, all_out{k1}.g_m_all(:,2),'k');
%     legend('IO', 'Con','True', 'Location','eastoutside')
%     xlabel('Time')
%     ylabel('Impuse Response')





figure
hold on
grid on
for k1 = 1:n_runs
    if k1 <= n_runs/2
        plot(all_out{k1}.cycle_error_each_iter, 'ro-')
    else
        plot(all_out{k1}.cycle_error_each_iter, 'bo-')
    end
    % plot(all_out{k1}.cycle_error_perfect, 'm')
end
ax = gca;
ax.YScale = 'log';
legend('Model Update', 'Model Init', 'Location','eastoutside')
xlabel('Cycle')
ylabel('Error in L2-Norm')


figure
hold on
grid on
for k1 = 1:n_runs
        plot(all_out{k1}.cycle_error_each_iter, 'ro-')
end
ax = gca;
ax.YScale = 'log';
xlabel('Cycle')
ylabel('Error in L2-Norm')





end