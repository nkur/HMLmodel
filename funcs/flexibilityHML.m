function [FME] = flexibilityHML(syn_range, dt, maxIter, sessions, trials, Synergies_orig, C_perb, gamma, eta, mu, k_p, sigma_u, sigma_q, a, epsilonBall, mcReps)

    % Targets
    targets = [0.5, 2.5, 2.5, 4.5;
               4.5, 0.5, 2.5, 4.5];
    tsIndex = [1, 2, 3, 4];
    joint_dim = 19;
        
    for syn_iter = 1:length(syn_range)
        FME_cum = 0;
        for mcRep = 1:mcReps
            
            syn = syn_range(syn_iter);

            clearvars x u delta_u W_hat C_hat Synergies Phi W_orig

            % Defining synergies and weights
            phi_orig            = Synergies_orig(1:syn, :);
            W_opt               = C_perb * pinv(phi_orig);
            W_hat(:, :, 1)      = 0.001* [zeros(2,syn-2), eye(2)];

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
                        clearvars x u delta_u W_hat C_hat
                        W_hat(:, :, 1) = W_hat_old(:, :, end);
                        C_hat(:, :, 1) = C_hat_old(:, :, end);            
                    end

                    if trial == 1
                          x(:, 1) = startPt;
                    else, x(:, 1) = x_old(:, end);
                    end

                    u(:, 1)    = zeros(joint_dim,1);
                    deltaQ(:, 1)    = zeros(joint_dim,1);

                    k = 1;

                    while norm(x(:, k) - endPt) > epsilonBall

                        x(:, k+1) = x(:, k) + dt * (C_perb) * u(:, k);

                        u(:, k+1) = u(:, k) - dt * eta * ( (phi_orig'*W_hat(:, :, k)' * W_hat(:, :, k)*phi_orig + mu * eye(joint_dim)) * u(:, k) - k_p * phi_orig' * W_hat(:, :, k)' * (endPt - x(:, k)) ) + sqrt(dt) * sigma_u * randn(joint_dim, 1); 

                        deltaQ(:, k+1) = deltaQ(:, k) + dt * (-a*deltaQ(:, k) + u(:, k)) + sqrt(dt) * sigma_q * randn(joint_dim, 1);

                        W_hat(:, :, k+1) = W_hat(:, :, k) - dt * gamma * (W_hat(:, :, k) - W_opt) * phi_orig * deltaQ(:, k) * deltaQ(:, k)' * phi_orig';
                        C_hat(:, :, k+1) = W_hat(:, :, k) * phi_orig;
%                         if k == maxIter
%                             k = 1;
%                             x(:, 2:end) = [];
%                             u(:, 2:end) = [];
%                             deltaQ(:, 2:end) = [];
%                             W_hat(:, :, 2:end) = [];
%                             C_hat(:, :, 2:end) = [];
%                             disp('Trial restart...');
% %                             restartTrial = restartTrial + 1;
%                         else
                        k = k+1;
%                         end

                    end

                    % Prepping for next trial
                    W_hat_old = W_hat;      C_hat_old = C_hat;      x_old = x;   

                end           
            end

%             FME_cum = FME_cum + norm(C_perb - C_hat(:, :, end))/norm(C_perb);
            FME_cum = FME_cum + norm(W_opt - W_hat(:, :, end))/norm(W_opt);
        end        
        FME(1, syn_iter) = FME_cum/mcReps;
        
%         if FME(1, syn_iter) > 1, FME(1, syn_iter) = NaN;    end
        
        
        dispp = [num2str(sigma_u), '-', num2str(syn)];
        disp(dispp)
    end
    
end

