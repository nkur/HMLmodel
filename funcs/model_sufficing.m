function [count_ex, nl_trial]=model_sufficing(dt, maxIterArr, sessions, trials, joint_dim, C, targets, tsIndex, params, a, epsilonBall_arr, Phi, W_orig, syn, eW_threshold, sigma_off)

    rng('shuffle', 'twister');

    gamma   = params(1);
    eta     = params(2);
    mu      = params(3);
    k_p     = params(4);
    sigma_u = params(5);
    sigma_q = params(6);

    W_hat(:, :, 1) = 0.001* [zeros(2,syn-2), eye(2)];

    count_ex = zeros(length(maxIterArr), length(epsilonBall_arr));  % Counter for cursor reaching each ex-ball
    nl_trial = 0;%zeros(length(maxIterArr), 1);  %nl_trial_flag = 0;              % Counters to calculate trials for which no learning occurs


    for session = 1:sessions
        for trial = 1:trials
            
            % Trial targets
            if trial == 1
                spIndex = 3;
            else, spIndex = epIndex;
            end

            epRand = randi(3);     nextIndex = tsIndex;   nextIndex(nextIndex == spIndex) = [];
            epIndex = nextIndex(epRand);
            startPt = targets(:, spIndex);
            endPt = targets(:, epIndex);
            
            % Initializing variables -----------
            if ~(session==1 && trial==1)
                clearvars x u W_hat C_hat
                W_hat(:, :, 1) = W_hat_old(:, :, end);
                C_hat(:, :, 1) = C_hat_old(:, :, end);
            end
    
            if trial == 1,     x(:, 1) = startPt;
            else,           x(:, 1) = x_old(:, end);   end
    
            u(:, 1)      = zeros(joint_dim,1);
            deltaQ(:, 1) = zeros(joint_dim,1);
            
            k = 1;
            
            %%%%% Trial starts --------------------------------------------
            while norm(x(:, k) - endPt) > 0.15   
                x(:, k+1) = x(:, k) + dt * C * u(:, k);

                u(:, k+1) = u(:, k) - dt * eta * ( (Phi'*W_hat(:, :, k)' * W_hat(:, :, k)*Phi + mu * eye(joint_dim)) * u(:, k) - k_p * Phi' * W_hat(:, :, k)' * (endPt - x(:, k)) ) + sqrt(dt) * sigma_u * randn(joint_dim, 1); 

                % Adding PE noise (\xi) to deltaQ
                deltaQ(:, k+1) = deltaQ(:, k) + dt * (-a*deltaQ(:, k) + u(:, k)) + sqrt(dt) * sigma_q * randn(joint_dim, 1);

                W_hat(:, :, k+1) = W_hat(:, :, k) - dt * gamma * (W_hat(:, :, k) - W_orig) * Phi * (deltaQ(:, k)) * (deltaQ(:, k))' * Phi';
                C_hat(:, :, k+1) = W_hat(:, :, k) * Phi;

                % Checking which epsilon balls are entered for which
                % trial-time indices after learning stops
                if  norm(W_orig - W_hat(:, :, k))/norm(W_orig) <= eW_threshold        % Sufficing condition
%                     sigma_q = sigma_off;
%                     sigma_u = sigma_off;
                    gamma = sigma_off;
%                     mu = sigma_off;

                    % Checking if any trial-time indices are hit
                    if any(k==maxIterArr)   % if any trial-time index is hit                    %%%%%%%%%%%% changed from <= to ==
                        time_idx = find(k==maxIterArr);  % Which trial-time index is hit
                        temp_cnt_arr = norm(x(:, end) - endPt) <= epsilonBall_arr; %which epsilon_arr is entered
                        for idx = time_idx
                            count_ex(idx, :) = count_ex(idx, :) + temp_cnt_arr;
%                             nl_trial = nl_trial + 1;
                        end
                    end
                end
                k = k+1;     
                if k >= 50/dt, break; end
            end
            %%%%% Trial ends ----------------------------------------------

            


            % Updating nl_trial if learning has stopped
            if gamma == sigma_off
                nl_trial = nl_trial + 1;

                % Add to the trial times which haven't been hit yet
                % Ex., if trajectory entered within 2s, it would have definitely entered in >2s
                time_idx = find(k<maxIterArr);
                temp_cnt_arr = norm(x(:, end) - endPt) <= epsilonBall_arr; %which epsilon_arr is entered
                for idx = time_idx
                    count_ex(idx, :) = count_ex(idx, :) + temp_cnt_arr;
                end
            end
    
            % Prepping for next trial
            W_hat_old = W_hat;      C_hat_old = C_hat;      x_old = x;            
        end        
    end

end