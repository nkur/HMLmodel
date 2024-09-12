% This script generates the results for Satisficing motor learning behavior for subject "01"
clc;    clearvars;  

mc_reps = 128*10;

load data_save_01.mat

syn         = 4;
sessions    = 8;
trials      = 60;
joint_dim   = 19;

% Targets
targets = [0.5, 2.5, 2.5, 4.5;
           4.5, 0.5, 2.5, 4.5];
tsIndex = [1, 2, 3, 4];

% Defining C, synergies and weights
C       = data_save.C;
Phi     = data_save.Phi;
W_orig  = C * pinv(Phi);

% Parameter initialization
params = [0.06638312 3.17414901 2.45812435 1.30978218 0.87637697 0.13695541];
a      = 10;

ex_array = [0.3, 0.5, 0.7, 0.9];
eW_threshold_array = 0.29 : 0.01 : 0.90; 
sigma_off = 0;

dt       = 0.01;
time_arr = 1 : 0.2 : 3;

% prob_ex = zeros(length(eW_threshold_array), length(time_arr), length(ex_array));

maxIterArr     = time_arr/dt;

for eW_iter = 1:length(eW_threshold_array)
    eW_threshold = eW_threshold_array(eW_iter);
    [prob_ex(eW_iter, :, :, :)] = sufficingUpFun(dt, maxIterArr, sessions, trials, joint_dim, C, targets, tsIndex, params, a, ex_array, Phi, W_orig, syn, eW_threshold, sigma_off, mc_reps);
end


%% New plot code
close all;
figTextSz = 30;

query_times = [1.2, 2.0, 3.0, 1.2]';
% [~, qt_idx] = find(time_arr == query_times);

% sbp1 = subplot(1,4,[1,2,3]);
ax = tiledlayout(1, length(query_times), "TileSpacing","tight");
for time_idx_i = 1:length(query_times)

    [~, time_idx] = find(time_arr == query_times(time_idx_i));

    prob_plot(:, :, :) = prob_ex(:, time_idx, :, :);
    time = time_arr(time_idx);


    nexttile;
    plot_AvgConfBds(eW_threshold_array, prob_plot);

    for i = 1:length(ex_array)
        lgdStr{2*i-1} = '';
        lgdStr{2*i} = sprintf('\\rho_x=%.1f', ex_array(i));
    end
%     legend(lgdStr, fontsize=10);

    titleStr = sprintf('Trial time: %.1fs', time);
    title(titleStr)
    ylim([0, 1])
    xlim([0.29, eW_threshold_array(end)])

    if time_idx_i == 1
        ylabel('Probability of Success', 'FontSize', figTextSz);
        legend(lgdStr, Location="northeast")

        % Select the zoom region
        rectangle('Position',[.29 .01 .31 .24],'EdgeColor','r', 'LineWidth',1)
        % daspect([1,1,1])
    else
        if time_idx_i ~= length(query_times)
            set(gca,'YTickLabel',[]);
        end
    end

    if time_idx_i == length(query_times)
        rectangle('Position',[.29 .01 .31 .24],'EdgeColor','r', 'LineWidth',2)
        % set(gca, 'YTickLabel', [0.1, 0.3, 0.5])
        xlim([0.29 0.6])
        ylim([0.01 0.25])
        titleStr = sprintf('Trial time: %.1fs (zoomed)', time);
        title(titleStr)
        % daspect([1,0.21/0.5,1])
    end

end
xlabel(ax, 'FME Threshold (Learning Accuracy)', 'FontSize', figTextSz, 'FontName','Times New Roman', 'interpreter','none')


% lg = legend(lgdStr, 'Orientation', 'Vertical');
% lg.Layout.Tile = 'east';
set(findall(gcf,'-property','FontSize'),'FontSize', figTextSz)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
set(gcf, "Position", [100         249        1757         606])