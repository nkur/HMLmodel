% This script generates the results and associated figures for Flexibility vs Performance
% motor learning behavior for subject "01". Use the data_save_01_syn.mat
% file as it contains all the synergies required to generate the results.

clc; close all; clearvars;

rng('shuffle');
load data_save_01_syn.mat

rand_length = 1;
rand_iter = 1;

syn_range   = 2:2:19;
time        = 50;
dt          = 0.01;
t           = dt : dt : time;
maxIter     = length(t)-1;
joint_dim   = 19;
mcReps      = 10;

sessions    = 4;
trials      = 60;

% Loading C and synergies
C_orig = data_save.C;
Synergies = data_save.Synergies;

% Perturbing C
perC    = 0.01;
C_perb  = -perC * ones(2, joint_dim) + 2*perC * rand(2, joint_dim) + C_orig;


% Parameters
params = [0.06638312 3.17414901 2.45812435 1.30978218 0.87637697 0.13695541];
gamma   = params(1);              
eta     = params(2);
mu      = params(3);
k_p     = params(4);
sigma_u = params(5);
sigma_q = params(6);

sigma_u_array = 0.1:0.1:2.0;
a        = 10;
epsilonBall = 0.15;
sigma_u_array_len = length(sigma_u_array);

parfor sigu_i = 1:sigma_u_array_len

    sigma_u = sigma_u_array(sigu_i);
    [FME_multiRun_sess(rand_iter, sigu_i, :)] = flexibilityHML(syn_range, dt, maxIter, sessions, trials, Synergies, C_perb, gamma, eta, mu, k_p, sigma_u, sigma_q, a, epsilonBall, mcReps);    
    disp(rand_iter)

end

FME_sess(:, :) = mean(FME_multiRun_sess, 1);

%% Heatmap
close all;
figure;
FME_sess_plot = FME_sess;
[val, ~] = min(FME_sess_plot, [], 2);
FME_sess_plot(FME_sess_plot == val) = NaN;
h = heatmap(syn_range, sigma_u_array, FME_sess_plot, 'Colormap', summer, 'FontSize', 24, 'FontName', 'Times New Roman');
h.MissingDataLabel = 'min(FME)';
h.MissingDataColor = [0.1 0.3 0.6];
h.CellLabelColor = 'none';
titleStr = ['FME Variation'];
h.Title = titleStr;
h.XLabel = 'Synergies';
h.YLabel = '\sigma_u';

 set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
% set(gcf, 'Position', [1402         861         978         850])
h.NodeChildren(3).YDir='normal';