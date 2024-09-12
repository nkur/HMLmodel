function [prob_ex_mc] = sufficingUpFun(dt, maxIterArr, sessions, trials, joint_dim, C, targets, tsIndex, params, a, epsilonBall_arr, Phi, W_orig, syn, eW_threshold, sigma_off, mc_reps)

    prob_ex_mc = zeros(length(maxIterArr), length(epsilonBall_arr), mc_reps);

    parfor mc_rep = 1:mc_reps

        rng('shuffle', 'twister')

        [count_ex, nl_trial] = model_sufficing(dt, maxIterArr, sessions, trials, joint_dim, C, targets, tsIndex, params, a, epsilonBall_arr, Phi, W_orig, syn, eW_threshold, sigma_off);        

        % Calculating per mc rep probabilities and summing up
%         nl_trial_arr = repmat(nl_trial, 1, size(count_ex, 2));
        prob_ex_mc(:, :, mc_rep) = count_ex./nl_trial;

    end

    % Calculating probabilities across all mc reps
%     prob_ex = prob_ex_mc/mc_reps;
    % pause;

    dispp = [num2str(eW_threshold)];
    disp(dispp)
end