function [cost] = pierella_model(params)

    epsilon = params(1);
    sigma   = params(2);

    load data_save_01.mat
    data_all = data_save.data_all;
    C = data_save.C;
    
    sessions = 8;
    trials = 60;

    H  = 0.01 * randn(2,19);
    H_ = H;
    cost = 0;

    for session = 1:sessions
        for trial = 1:trials
            endPt = data_all(session, trial).endPt;

            tEval = 2;
            idx = find(data_all(session, trial).time >= tEval, 1);

            cursor_p = data_all(session, trial).x(:, idx);
            
            cursor_err = cursor_p - endPt;

            body_q = data_all(session, trial).q(:, idx) + sigma * randn(19, 1);

            H_ = H_ + epsilon * (cursor_p - H_ * body_q) * body_q';
            H = cat(3, H, H_);

            FME = norm(C - H_)/norm(C);

            cost = cost + FME;

        end
    end
end