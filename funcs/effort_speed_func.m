function [drvEffPlot, expEffPlot, speedPlot, accPlot] = effort_speed_func(params)

    load data_save_01.mat
    
    % Defining synergies and weights
    C = data_save.C;
    Phi = data_save.Phi;
    W_orig = C * pinv(Phi);
    
    gamma   = params(1);
    eta     = params(2);
    mu      = params(3);
    k_p     = params(4);
    sigma_u = params(5);
    sigma_q = params(6);
    a       = 10;

    epsilonBall = 0.15;
    
    syn         = 4;
    time        = 50;
    dt          = 0.01;
    t           = dt : dt : time;
    maxIter     = length(t)-1;
    sessions    = 8;
    trials      = 60;
    joint_dim   = 19;
    restartTrial = 1;
    
    % Targets
    targets = [0.5, 2.5, 2.5, 4.5;
               4.5, 0.5, 2.5, 4.5];
    tsIndex = [1, 2, 3, 4];
    
    W_hat(:, :, 1) = 0.001* [zeros(2,syn-2), eye(2)];
    

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
            deltaQ    = zeros(joint_dim,1);
            
            k = 1;
            
            while norm(x(:, k) - endPt) > epsilonBall
    
                x(:, k+1) = x(:, k) + dt * C * u(:, k);
    
                u(:, k+1) = u(:, k) - dt * eta * ( (Phi'*W_hat(:, :, k)' * W_hat(:, :, k)*Phi + mu * eye(joint_dim)) * u(:, k) - k_p * Phi' * W_hat(:, :, k)' * (endPt - x(:, k)) ) + sqrt(dt) * sigma_u * randn(joint_dim, 1); 

                % Adding PE noise (\xi) to deltaQ
                deltaQ(:, k+1) = deltaQ(:, k) + dt * (-a*deltaQ(:, k) + u(:, k)) + sqrt(dt) * sigma_q * randn(joint_dim, 1);
    
                W_hat(:, :, k+1) = W_hat(:, :, k) - dt * gamma * (W_hat(:, :, k) - W_orig) * Phi * deltaQ(:, k) * deltaQ(:, k)' * Phi';
                C_hat(:, :, k+1) = W_hat(:, :, k) * Phi;

                if k == maxIter
                    k = 1;
                    x(:, 2:end) = [];
                    u(:, 2:end) = [];
                    deltaQ(:, 2:end) = [];
                    W_hat(:, :, 2:end) = [];
                    C_hat(:, :, 2:end) = [];
                    disp('Trial restart...');
                    restartTrial = restartTrial + 1;
                else
                    k = k+1;
                end
    
            end
    
            % Prepping for next trial
            W_hat_old = W_hat;      C_hat_old = C_hat;      x_old = x;
            
            plot_time = dt:dt:k*dt;
    
            % Driving Effort vs Exporatory Effort
            [~, maxInd] = max(vecnorm(u));
            [drivingEff(session, trial), explorEff(session, trial)] = effortComp(C_hat(:, :, maxInd), u(:, maxInd));
    
            % Accuracy vs Speed
            [accInv(session, trial), speedInv(session, trial)] = accComp(x, startPt, endPt, plot_time);
        
        end                   
    end
    
    % Performance Measures - Driving vs Exploratory Effort
    tmpDrvEff = drivingEff';    tmpExpEff = explorEff';
    drvEffPlot = tmpDrvEff(:);  expEffPlot = tmpExpEff(:);
    
    
    % Accuracy vs Speed
    tmpAccInv = accInv';          tmpSpeedInv = speedInv';
    accPlot = 1./tmpAccInv(:);    speedPlot = 2000./tmpSpeedInv(:);

end