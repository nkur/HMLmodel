% This script fits the model proposed by "Pierella etal. (2019) The
% dynamics of motor learning through the formation of internal models. PLoS
% Comput Biol" to our experiment data, and then compares the error in RE
% from model and data.

% Run first section to generate fit files. Then run second section to
% compare RE error between two models
clc;
clearvars;
close all;

rng('shuffle');

subID = '01';
loadFileStr = "data_save_" + subID + ".mat";
load(loadFileStr)

data_all = data_save.data_all;

sessions = 8;
trials = 60;

C = data_save.C;
Phi = data_save.Phi;
W_orig = C / Phi;

RE_sub = [];    IME_sub = [];
Q = [];         U = [];         G = [];

% Getting RE and IME from exp data
for session = 1:sessions
    for trial = 1:trials
        
        endPt = data_all(session, trial).endPt;
        tEval = 2;
        idx = find(data_all(session, trial).time >= tEval, 1);
        RE_sub(end+1) = norm(data_all(session, trial).x(:, idx) - endPt);

        % Computing inv(C_hat) = G by performing LS on windowed data
        if session == 1 && trial <= 12
            Q = [Q, data_all(session, trial).q(:,idx)];
            U = [U, endPt];
        else
            G = cat(3, G, Q * U' / (U*U'));
            Q = [Q(:, 2:end), data_all(session, trial).q(:,idx)];
            U = [U(:, 2:end), endPt];
            
            % Computing IME
            IME_sub(end+1) = norm(C*G(:,:,end)-eye(2));
        end

        
    end
end
RE_sub_smoo = smootherFun(RE_sub, 10);
% RE_sub_smoo = RE_sub;
trialVec = 1:length(RE_sub_smoo);

% Creating LSfit model and setting options
fo = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-5, 0, 0], 'Upper', [10, 10, 10]);
ft =  fittype('a*exp(-b*x) + c', 'options', fo);

[curve, gof] = fit(trialVec', RE_sub_smoo', ft)

figure;
plot(curve, '-r', trialVec, RE_sub_smoo, '-k')

figure;
plot(IME_sub)

% Fitting sigma and epsilon by using sum_FME as cost and fminunc as optimizer

reps = 128;
foptions = optimset("MaxIter", 1000, "MaxFunEvals", 10000);

fun = @pierella_model;
parfor rep = 1:reps
    x0 = 0.1 * rand(1,2);
    [params(rep, :), fval(rep)] = fminsearch(fun, x0, foptions);
end

[~, min_idx] = min(fval);
params_opt = params(min_idx, :);

% Propagating the Pierella Model based on fitted param values
epsilon = params_opt(1);
sigma = params_opt(2);
eta = curve.b;

H  = 0.01 * rand(2,19);        H_ = H;
G = 0.01 * pinv(H);            G_ = G;

RE_true = [];       RE = [];
FME = [];
for session = 1:sessions
    for trial = 1:trials
        endPt = data_all(session, trial).endPt;

        tEval = 2;
        idx = find(data_all(session, trial).time >= tEval, 1);
        RE_true(end+1) = norm(data_all(session, trial).x(:, idx) - endPt);

        body_q = G(:, :, end) * endPt + sigma * randn(19, 1);

        cursor_p = C * body_q;

        H_ = H_ + epsilon * (cursor_p - H_ * body_q) * body_q';
        H = cat(3, H, H_);

        G_ = G_ - eta * H_' * (cursor_p - endPt) * endPt';
        G = cat(3, G, G_);

        RE(end+1) = norm(cursor_p - endPt);

        FME(end+1) = norm(C - H_)/norm(C);

    end
end

figure;
plot(smootherFun(RE, 10), '-b');
hold on; plot(smootherFun(RE_true, 10), '-r');
grid on;
legend('Model', 'Subject');
titleStr = "RE Trend" + subID;
title(titleStr);

figure;
plot(FME);
titleStr = "FME Trend" + subID;
title(titleStr);
grid on;

norm(smootherFun(RE, 10)- smootherFun(RE_true, 10), 2)

%
saveFileStr = "pierella_model_sub_" + subID;
save(saveFileStr)

%% Computing RE fitting error
clearvars;
close all;

subjects = ["01" "02" "03" "04" "05" "06"];

for sub = 1:length(subjects)

    subID = subjects(sub);
    loadFileStr = "pierella_model_sub_" + subID + ".mat";
    load(loadFileStr)
    
    epsilon = params_opt(1);
    sigma = params_opt(2);
    eta = curve.b;
    
    %%% Multiple runs for evaluating the RE fitting error
    mcreps = 500;
    parfor mcrep = 1:mcreps
        H = 0.001 * rand(2,19);        H_ = H;
        G = 0.001 * pinv(H);           G_ = G;
        
        RE_true = [];       RE = [];
        FME = [];
        for session = 1:sessions
            for trial = 1:trials
                endPt = data_all(session, trial).endPt;
        
                tEval = 2;
                idx = find(data_all(session, trial).time >= tEval, 1);
                RE_true(end+1) = norm(data_all(session, trial).x(:, idx) - endPt);
        
                body_q = G(:, :, end) * endPt + sigma * randn(19, 1);
        
                cursor_p = C * body_q;
        
                H_ = H_ + epsilon * (cursor_p - H_ * body_q) * body_q';
                H = cat(3, H, H_);
        
                G_ = G_ - eta * H_' * (cursor_p - endPt) * endPt';
                G = cat(3, G, G_);
        
                RE(end+1) = norm(cursor_p - endPt);
        
                FME(end+1) = norm(C - H_)/norm(C);
        
            end
        end
    
        RE_error(sub, mcrep) = norm(smootherFun(RE, 10)- smootherFun(RE_true, 10), 2)/norm(smootherFun(RE_true, 10), 2);
        [corcoff, ~, rl, ru] = corrcoef(smootherFun(RE, 10), smootherFun(RE_true, 10));
        Rsq(sub, mcrep) = corcoff(1,2);
        Rsq_l(sub, mcrep) = rl(1,2);
        Rsq_u(sub, mcrep) = ru(1,2);
    end
end
mean(RE_error, 2)
mean(Rsq, 2)

figure;
errorbar(mean(Rsq, 2), mean(Rsq_u, 2) - mean(Rsq_l, 2), 'b', 'LineWidth',2)
xlim([0, 7])