% This script generates the RE and SOT plots presented in the paper. Most
% of the code here is also contained in HML_model.m

clc; clearvars;

rng('shuffle');

subjects = ["01", "05", "08", "12", "14", "16"];
params_all = [0.06638312 3.17414901 2.45812435 1.30978218 0.87637697 0.13695541;
              3.02844779e-03 3.14480612e+00 3.30565657e+00 1.59646841e+00 1.01647669e+00 5.45089277e-01;
              0.04560625 1.53826924 3.30723253 3.27134663 1.00814878 0.05076906;
              0.13975294 1.9856537  3.57353247 1.89765196 0.9555619  0.01692007;
              1.27167635e-03 2.49159204e+00 3.53825520e+00 1.55688299e+00 1.97490096e+00 7.11824456e-01;
              0.12519652 0.71310447 3.97435408 2.25149752 0.92981486 0.00643546];
    
RE = zeros(length(subjects), 8*60);    RE_true = RE;
SOT = RE;   SOT_true = RE;

for subID = 1:length(subjects)

    clearvars -except subjects params_all RE RE_true SOT SOT_true subID
    subject = subjects(1,subID);

    loadFileStr = "data_save_" + subject + ".mat";
    load(loadFileStr)
    data_all = data_save.data_all;

    dt = 0.01;
    syn     = 4;
    sessions  = 8;
    trials = 60;
    joint_dim   = 19;
    
    C = data_save.C;
    Phi = data_save.Phi;
    W_orig = C / Phi;
    
    W_hat(:, :, 1) = 0.0001* [zeros(2,syn-2), eye(2)];
    
    % Parameter initialization
    params  = params_all(subID, :);
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
                clearvars x u W_hat C_hat subsamp_p_hml subsamp_Q_hml deltaQ plot_time
                W_hat(:, :, 1) = W_hat_old(:, :, end);
                C_hat(:, :, 1) = C_hat_old(:, :, end);
            end
    
            if trial == 1
                  x(:, 1) = startPt;
            else, x(:, 1) = x_old(:, end);
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
                k = k+1;          
    
            end
            %%%%%% --------------------------------------------------------
    
            % Prepping for next trial
            W_hat_old = W_hat;      C_hat_old = C_hat;      x_old = x;
    
            plot_time = dt:dt:k*dt;   
    
            %%%%% Cost function and FME computation
            FME = norm(W_orig - W_hat(:, :, end))/(norm(W_orig));
            
            FME_trials(session,trial)     = FME;
            
            %%% Subsampling for Comp %%%%%%%%%%%%%%%%%%%%%
            p_hml = x; % - x(:, 1);        % Not converting to delta_x
            subsamp_p_hml(1,:) = interp1(plot_time, p_hml(1,:), data_all(session,trial).time);
            subsamp_p_hml(2,:) = interp1(plot_time, p_hml(2,:), data_all(session,trial).time);
            subsamp_p_hml(isnan(subsamp_p_hml)) = 0;
    
            %%%%% RE at the end of ballistic phase
            tEval = 2;
            idx_true = find(data_all(session,trial).time >= tEval, 1);
            idx = idx_true; %find(plot_time >= tEval, 1);

            trialNum = (session-1)*60 + trial;
            
            RE(subID, trialNum)      = norm(subsamp_p_hml(:, idx) - endPt);
            RE_true(subID, trialNum) = norm(data_all(session,trial).x(:, idx_true) - endPt);
    
            %%%%% Computing SOT for the current trial
            SOT_true(subID, trialNum) = maxPerpendicularDist(data_all(session,trial).x, startPt, endPt)/norm(endPt - startPt);
            SOT(subID, trialNum)      = maxPerpendicularDist(subsamp_p_hml, startPt, endPt)/norm(endPt - startPt); 
    
            dispp = [num2str(session), '-', num2str(trial)];
            disp(dispp)
    
        end    
    end
end

%% Plots
close all;
figure;

figTextSize = 25;
t = tiledlayout(1,6);
for subID = 1:6

    nexttile;
    plot(smootherFun(RE(subID,:), 10))
    hold on; plot(smootherFun(RE_true(subID,:), 10), '-r');
    legend('Model', 'Subject')
    grid on;
    titleStr = ['Subject ', num2str(subID)];
    title(titleStr)
    ylim([0, 3])
    yticks([0:0.5:3])
    xlim([0, 480-10])
    set(findobj(gca,'type','line'),'linew',2)
    set(findall(gcf,'-property','FontSize'),'FontSize',figTextSize)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

    if subID == 1
        ylabel('Reaching Error');
    else
        yticklabels([]);
        yticks([0:0.5:3])
    end
    xlabel('Trials')
end
set(gcf, 'Position', [1711        1812        2852         400])
t.Padding = 'none';
t.TileSpacing = 'tight';


figure;
figTextSize = 25;
t = tiledlayout(1,6);
for subID = 1:6

    nexttile;
    plot(smootherFun(SOT(subID,:), 10))
    hold on; plot(smootherFun(SOT_true(subID,:), 10), '-r');
    legend('Model', 'Subject')
    grid on;
    titleStr = ['Subject ', num2str(subID)];
    title(titleStr)
    yticks([0:0.3:1.5])
    xlim([0, 480-10])
    set(findobj(gca,'type','line'),'linew',2)
    set(findall(gcf,'-property','FontSize'),'FontSize',figTextSize)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

    if subID == 1
        ylabel('Straightness of Trajectory');
    
    else
        yticklabels([]);
        yticks([0:0.3:1.5])
    end
    xlabel('Trials')
end
set(gcf, 'Position', [1711        812        2852         400])
t.Padding = 'none';
t.TileSpacing = 'tight';







