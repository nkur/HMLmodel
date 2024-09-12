function [RE, RE_true, SOT, SOT_true, varargout] = HMLmodel_RESOT(params, data_save)

    rng('shuffle');

    data_all = data_save.data_all;
    
    C = data_save.C;
    Phi = data_save.Phi;
    W_orig = C / Phi;
    
    dt          = 0.01;
    syn         = size(Phi, 1);
    sessions    = 8;
    trials      = 60;
    joint_dim   = size(Phi, 2);
    
    % Parameter initialization
    gamma   = params(1);
    eta     = params(2);
    mu      = params(3);
    k_p     = params(4);
    sigma_u = params(5);
    sigma_q = params(6);
    a       = 10;
    
    W_hat(:, :, 1) = 0.0001* [zeros(2,syn-2), eye(2)];
    
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
            
            RE(trialNum)      = norm(subsamp_p_hml(:, idx) - endPt);
            RE_true(trialNum) = norm(data_all(session,trial).x(:, idx_true) - endPt);
    
            %%%%% Computing SOT for the current trial
            SOT_true(trialNum) = maxPerpendicularDist(data_all(session,trial).x, startPt, endPt)/norm(endPt - startPt);
            SOT(trialNum)      = maxPerpendicularDist(subsamp_p_hml, startPt, endPt)/norm(endPt - startPt); 
    
            FME(trialNum) = norm(C - C_hat(:,:,end))/norm(C);
    
    %         dispp = [num2str(session), '-', num2str(trial)];
    %         disp(dispp)
    
        end    
    end

    varargout{1} = FME;


end