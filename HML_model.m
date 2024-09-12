% This script runs the HML model for all subjects and compares the model generated cursor trajectories, RE, and SOT with the
% subject data. It also generates the FME plots for all subjects.
clc; clearvars;

rng('shuffle');
addpath(genpath(pwd));

subjects = ["01", "02", "03", "04", "05", "06"];
paramsList_all  =  [0.06638312 3.17414901 2.45812435 1.30978218 0.87637697 0.13695541;
                    3.02844779e-03 3.14480612e+00 3.30565657e+00 1.59646841e+00 1.01647669e+00 5.45089277e-01;
                    0.04560625 1.53826924 3.30723253 3.27134663 1.00814878 0.05076906;
                    0.13975294 1.9856537  3.57353247 1.89765196 0.9555619  0.01692007;
                    1.27167635e-03 2.49159204e+00 3.53825520e+00 1.55688299e+00 1.97490096e+00 7.11824456e-01;
                    0.12519652 0.71310447 3.97435408 2.25149752 0.92981486 0.00643546];

for subID_idx = 1:length(subjects)

    clearvars -except subjects paramsList_all FME_sub_trials subID_idx data RE RE_true SOT SOT_true;
    
    subID = subjects(subID_idx);
    loadSubStr = "data_save_" + subID + ".mat";
    load(loadSubStr)
    data_all = data_save.data_all;
    
    dt          = 0.01;
    syn         = 4;
    sessions    = 8;
    trials      = 60;
    joint_dim   = 19;
    
    C       = data_save.C;
    Phi     = data_save.Phi;
    W_orig  = C / Phi;
    
    W_hat(:, :, 1) = 0.0001* [zeros(2,syn-2), eye(2)];
    
    % Parameter initialization
    params  = paramsList_all(subID_idx, :);
    gamma   = params(1);
    eta     = params(2);
    mu      = params(3);
    k_p     = params(4);
    sigma_u = params(5);
    sigma_q = params(6);
    a       = 10;
    
    % -------------------------------------------------------------------------
    for session = 1:sessions    
        for trial = 1:trials
    
            trialTime = data_all(session, trial).time(end);
            t_trial = dt: dt: trialTime;
            maxIter = length(t_trial) - 1;
                    
            startPt = data_all(session, trial).startPt;
            endPt = data_all(session, trial).endPt;
    
            if ~(session==1 && trial==1)
                clearvars x u W_hat C_hat q subsamp_x_hml subsamp_q_hml subsamp_Q_hml deltaQ plot_time
                W_hat(:, :, 1) = W_hat_old(:, :, end);
                C_hat(:, :, 1) = C_hat_old(:, :, end);
            end
    
            if trial == 1
                  x(:, 1) = startPt;
                  q(:, 1) = pinv(C) * x(:,1);
            else
                x(:, 1) = x_old(:, end);
                q(:, 1) = q_old(:, end);
            end
    
            u(:, 1)         = zeros(joint_dim,1);
            deltaQ(:, 1)    = zeros(joint_dim,1);
    
            k = 1;
    
            %%%%% Trial run starts  ---------------------------------------
            while k<maxIter
    
                x(:, k+1) = x(:, k) + dt * C * u(:, k);
    
                u(:, k+1) = u(:, k) - dt * eta * ( (Phi'*W_hat(:, :, k)' * W_hat(:, :, k)*Phi + mu * eye(joint_dim)) * u(:, k) - k_p * Phi' * W_hat(:, :, k)' * (endPt - x(:, k)) ) + sqrt(dt) * sigma_u * randn(joint_dim, 1); 
    
                % Adding PE noise (\xi) to deltaQ
                deltaQ(:, k+1) = deltaQ(:, k) + dt * (-a*deltaQ(:, k) + u(:, k)) + sqrt(dt) * sigma_q * randn(joint_dim, 1);
    
                W_hat(:, :, k+1) = W_hat(:, :, k) - dt * gamma * (W_hat(:, :, k) - W_orig) * Phi * deltaQ(:, k) * deltaQ(:, k)' * Phi';
                C_hat(:, :, k+1) = W_hat(:, :, k) * Phi;
    
                q(:, k+1) = q(:, k) + dt * u(:, k);
                k = k+1;          
    
            end
            %%%%%% --------------------------------------------------------
    
            % Prepping for next trial
            W_hat_old = W_hat;      C_hat_old = C_hat;      x_old = x;      q_old = q;
    
            plot_time = dt:dt:k*dt;   
    
            %%%%% Cost function and FME computation
            FME = norm(C - C_hat(:, :, end))/(norm(C));
            
            %%%%% Storing data
            data(subID_idx,session,trial).x         = x;
            data(subID_idx,session,trial).q_dot     = u;
            data(subID_idx,session,trial).W_hat     = W_hat;
            data(subID_idx,session,trial).time      = plot_time;
            data(subID_idx,session, trial).startPt  = startPt;
            data(subID_idx,session, trial).endPt    = endPt;

            trial_idx = (session-1)*trials + trial;
            FME_sub_trials(subID_idx, trial_idx) = FME;
            
            %%% Subsampling for Comp %%%%%%%%%%%%%%%%%%%%%
            p_hml = x;
            subsamp_x_hml(1,:) = interp1(plot_time, p_hml(1,:), data_all(session,trial).time);
            subsamp_x_hml(2,:) = interp1(plot_time, p_hml(2,:), data_all(session,trial).time);
            subsamp_x_hml(isnan(subsamp_x_hml)) = 0;

            % Subsampling the joint angles
            for joint = 1:joint_dim
                subsamp_q_hml(joint,:) = interp1(plot_time, q(joint, :), data_all(session,trial).time);
            end
    
            %%%%% RE at the end of ballistic phase (2s)
            tEval = 2;
            idx_true = find(data_all(session,trial).time >= tEval, 1);            
            if isempty(idx_true)
                idx_true = maxIter;
            else
                idx_true = min(maxIter, idx_true);
            end
            idx = idx_true;
            
            RE(subID_idx, trial_idx) = norm(subsamp_x_hml(:, idx) - endPt);
            RE_true(subID_idx, trial_idx) = norm(data_all(session,trial).x(:, idx_true) - endPt);
    
            %%%%% Computing SOT for the current trial
            SOT_true(subID_idx, trial_idx) = trajStraightness(data_all(session,trial).x, startPt, endPt);
            SOT(subID_idx, trial_idx)      = trajStraightness(subsamp_x_hml, startPt, endPt); 
    
    
            dispp = ['Subject: ', num2str(subID_idx), '- Session: ', num2str(session)];
            disp(dispp)
    
        end    
    end
end


%% Plots for visualization
close all;
% Choose subject to display the figures for
subID_idx = 1;

subID = subjects(subID_idx);
loadSubStr = "data_save_" + subID + ".mat";
load(loadSubStr)
data_all = data_save.data_all;

figTextSize = 30;

% Plotting trajectory tiles with RE -------------------------------------
firstPlot_flag = 1;

figure;
t1 = tiledlayout(1,4);

for ii = [1,8]
    nexttile;
    for jj =1:trials
        plot(data(subID_idx, ii, jj).x(1, :), data(subID_idx, ii, jj).x(2, :), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
        plot(data_all(ii, jj).startPt(1), data_all(ii, jj).startPt(2), 'r.', data_all(ii, jj).endPt(1), data_all(ii, jj).endPt(2), 'r.', 'MarkerSize',50)
        hold on
    end 
    axis([-2 8 -1 6.5])
    titleStr = ['Session ', num2str(ii), ' Trials (Model)'];
    title(titleStr)

    set(findall(gcf,'-property','FontSize'),'FontSize',figTextSize)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

    if firstPlot_flag
        firstPlot_flag = 0;
    else
        set(gca,'YTickLabel',[]);
    end
    
    nexttile;
    for jj = 1:trials
        plot(data_all(ii, jj).x(1, :), data_all(ii, jj).x(2, :), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
        plot(data_all(ii, jj).startPt(1), data_all(ii, jj).startPt(2), 'r.', data_all(ii, jj).endPt(1), data_all(ii, jj).endPt(2), 'r.', 'MarkerSize',50)
        hold on
    end 
    axis([-2 8 -1 6.5])
    titleStr = ['Session ', num2str(ii), ' Trials (Subject)'];
    title(titleStr)
    set(findall(gcf,'-property','FontSize'),'FontSize',figTextSize)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

    set(gca,'YTickLabel',[]);
    
end
titleStr = ['Subject ', subID];
title(t1, titleStr, fontsize=40, fontname='Times New Roman');
t1.Padding = 'none';
t1.TileSpacing = 'tight';

figure;
plot(smootherFun(RE(subID_idx,:), 10))
hold on; plot(smootherFun(RE_true(subID_idx,:), 10), '-r');
legend('Model', 'Subject')
grid on;
xlabel('Trials')
ylabel('Reaching Error')
title('Reaching Error')
ylim([0, 3])
xlim([0, 480-10])
set(findobj(gca,'type','line'),'linew',2)
set(findall(gcf,'-property','FontSize'),'FontSize',figTextSize)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

%% FME Plot all subjects
figure;
t2 = tiledlayout(1,6);

for subID_idx = 1:length(subjects)
    nexttile;
    plot(FME_sub_trials(subID_idx, :)*100, 'linewidth', 2)
    if subID_idx == 1
        ylabel('Forward Modeling Error(%)')
    else
        set(gca,'YTickLabel',[]);
    end
    xlabel('Trials')
    titleStr = ['Subject ', num2str(subID_idx)];
    title(titleStr, fontsize=figTextSize, fontname='Times New Roman');
    grid on;
    ylim([0, 100])
    xlim([0, 480])
    set(findall(gcf,'-property','FontSize'),'FontSize',figTextSize)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

end
t2.Padding = 'tight';
t2.TileSpacing = 'compact';

