% This script analyzes the exploration-exploitation, and speed-accuracy
% tradeoffs for subject "01". To change the subject under study, go to the
% effort_speed_func.m
% Various statistical tests are also run at the end of script to support
% the hypothesis of maxima and plateauing around the respective fit values.

clc; clearvars;  

rng('shuffle');
addpath(genpath(pwd));

mc_reps = 128;

param_test = 0.5:0.5:30;
len_param_test = length(param_test);

% Plot flags
effort = 1;
spedAcc = 1;

study_parameter = '\eta';

parfor mc_rep = 1:mc_reps
    for test_var = 1:len_param_test
        params = [0.06638312 3.17414901 2.45812435 1.30978218 0.87637697 0.13695541];
    
        % Change the params index to the corresponding parameter under study -------------
        params(2) = param_test(test_var);
    
        [drvEffPlot_mc(mc_rep, :, test_var), expEffPlot_mc(mc_rep, :, test_var), speedPlot_mc(mc_rep, :, test_var), accPlot_mc(mc_rep, :, test_var)] = effort_speed_func(params);    
    
    end
end

%% Mean-trial metrics plots with CI across mcreps on 
close all;
drvEffPlot = squeeze(mean(drvEffPlot_mc, 1));        expEffPlot = squeeze(mean(expEffPlot_mc, 1));
speedPlot = squeeze(mean(speedPlot_mc, 1));          accPlot = squeeze(mean(accPlot_mc, 1));

% Effort
dEff_mean   = mean(drvEffPlot, 1);
dEff_ci     = 1.96 * std(drvEffPlot, 0)/sqrt(size(drvEffPlot, 1));
dEff_lb     = dEff_mean - dEff_ci;
dEff_ub     = dEff_mean + dEff_ci;

eEff_mean   = mean(expEffPlot, 1);
eEff_ci     = 1.96 * std(expEffPlot, 0)/sqrt(size(expEffPlot, 1));
eEff_lb     = eEff_mean - eEff_ci;
eEff_ub     = eEff_mean + eEff_ci;

figure;
ha = tight_subplot(2,1, 0.05, [.15 .05], [.15 .1]);
axes(ha(1))
plot(param_test, dEff_mean, 'LineWidth',2);
hold on;
p_drv = fill([param_test, param_test(end:-1:1)], [dEff_lb, dEff_ub(end:-1:1)], 'blue');
p_drv.FaceColor = [0.8 0.8 1];    p_drv.EdgeColor = 'none';     p_drv.FaceAlpha = 0.8;
hold off;
ylabel('Driving Effort')
xticks([param_test(1):1:param_test(end)])
xticklabels('')
xlim([param_test(1) param_test(end)]);
% ylim([0.6 1.5])
grid on;

axes(ha(2))
plot(param_test, eEff_mean, 'LineWidth',2)
hold on;
p_exp = fill([param_test, param_test(end:-1:1)], [eEff_lb, eEff_ub(end:-1:1)], 'blue');
p_exp.FaceColor = [0.8 0.8 1];    p_exp.EdgeColor = 'none';     p_exp.FaceAlpha = 0.8;
hold off;
ylabel('Exploratory Effort')
xlabel(study_parameter)
xticks([param_test(end):1:param_test(end)])
xlim([param_test(1) param_test(end)]);
% ylim([0.5 3.5])
grid on;

set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

% Plotting differences distribution
clearvars pVal_d pVal_e
qp_dt = 0.1;
qp_1 = 1:0.1:param_test(end)-qp_dt;
[~, idx1] = find(abs(param_test - qp_1') < 1e-5);
[~, idx2] = find(abs(param_test - [qp_1 + qp_dt]') < 1e-5);
tmp = [idx1'; idx2'];
idx_dtDiff = tmp(:);

dEff_diff = abs(drvEffPlot(:, idx2) - drvEffPlot(:, idx1));
eEff_diff = abs(expEffPlot(:, idx2) - expEffPlot(:, idx1));

for i = 1:length(idx1)-1
%     [~, pVal_d(i)] = ttest2(dEff_diff(:,i), dEff_diff(:, i+1), 'Tail', 'right');    
    i_odd = 2*i-1;      i_even = 2*i;
    [~, pVal_d(i)] = ttest2(drvEffPlot(:, idx_dtDiff(i_odd)), drvEffPlot(:, idx_dtDiff(i_even)), 'Tail', 'both');
end
figure;
% boxplot(dEff_diff)
pt_plot = param_test(idx1);
plot(pt_plot(1:end-1), pVal_d, '.-k', 'MarkerSize',20);
xlblStr = "\eta - (\eta+" + num2str(qp_dt) + ") value pairs";
xlabel(xlblStr)
% ylabel("p-values (MD_{\eta}(Eff_{drv}) > MD_{\eta^+}(Eff_{drv}))")
ylabel("p-values (H_0: Eff_{drv}(\eta) < Eff_{drv}(\eta^+))")
xlim([pt_plot(1) pt_plot(end-1)])

for i = 1:length(idx1)-1
%     [~, pVal_e(i)] = ttest2(eEff_diff(:,i), eEff_diff(:, i+1), 'Tail', 'right');
    i_odd = 2*i-1;      i_even = 2*i;
    [~, pVal_e(i)] = ttest2(expEffPlot(:, idx_dtDiff(i_odd)), expEffPlot(:, idx_dtDiff(i_even)), 'Tail', 'both');

end
figure;
% boxplot(eEff_diff)
plot(pt_plot(1:end-1), pVal_e, '.-k', 'MarkerSize',20)
xlblStr = "\eta - (\eta+" + num2str(qp_dt) + ") value pairs";
xlabel(xlblStr)
% ylabel("p-values (MD_{\eta}(Eff_{exp}) > MD_{\eta^+}(Eff_{exp}))")
ylabel("p-values (H_0: Eff_{exp}(\eta) > Eff_{exp}(\eta^+))")
xlim([pt_plot(1) pt_plot(end-1)])

%% ANOVA on Driving Effort
qp_arr = [3 : 1 : 10]';
[~, idx] = find(param_test == qp_arr);
drvEff_qp = drvEffPlot(:, idx);
aov = anova1(drvEff_qp)

%% Speed-Accuracy

accu_mean   = mean(accPlot, 1);
accu_ci     = 1.96 * std(accPlot, 0)/sqrt(size(accPlot, 1));
accu_lb     = accu_mean - accu_ci;
accu_ub     = accu_mean + accu_ci;

sped_mean   = mean(speedPlot, 1);
sped_ci     = 1.96 * std(speedPlot, 0)/sqrt(size(speedPlot, 1));
sped_lb     = sped_mean - sped_ci;
sped_ub     = sped_mean + sped_ci;

figure;
ha = tight_subplot(2,1, 0.05, [.15 .05], [.15 .1]);
axes(ha(1))
plot(param_test, accu_mean, 'LineWidth',2);
hold on;
p_acu = fill([param_test, param_test(end:-1:1)], [accu_lb, accu_ub(end:-1:1)], 'blue');
p_acu.FaceColor = [0.8 0.8 1];    p_acu.EdgeColor = 'none';     p_acu.FaceAlpha = 0.8;
hold off;
ylabel('Accuracy')
xticks([0:1:param_test(end)])
xticklabels('')
xlim([param_test(1) param_test(end)]);  ylim([3 7])
grid on;

axes(ha(2))
plot(param_test, sped_mean, 'LineWidth',2)
hold on;
p_spd = fill([param_test, param_test(end:-1:1)], [sped_lb, sped_ub(end:-1:1)], 'blue');
p_spd.FaceColor = [0.8 0.8 1];    p_spd.EdgeColor = 'none';     p_spd.FaceAlpha = 0.8;
hold off;
ylabel('Speed')
xlabel(study_parameter)
xticks([0:1:param_test(end)])
xlim([param_test(1) param_test(end)]);  ylim([0 13])
grid on;

set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

% Plotting differences distribution
clearvars pVal_s pVal_a
query_param = [1.3, 0.3, 2.3, 3.3]';
for qp_i = 1:length(query_param)
    qp = query_param(qp_i);
    [~, idxs(qp_i)] = find(abs(param_test - qp) < 1e-5);
end
% acc_diff = abs(accPlot(:, indices(2:end)) - accPlot(:, indices(1:end-1)));
% sped_diff = abs(speedPlot(:, indices(2:end)) - speedPlot(:, indices(1:end-1)));

for i = 1 : length(idxs)-1
    [~, pVal_a(i)] = ttest2(accPlot(:,idxs(1)), accPlot(:, idxs(i+1)), 'Tail', 'right');
end
figure;
plot(pVal_a, '.-k', 'MarkerSize',20)


