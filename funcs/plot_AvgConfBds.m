function [] = plot_AvgConfBds(x_var, y_var)

    cmap = parula(6);
%     cmap = [1, 0, 0];

    hold on;
    for varNum = 1:size(y_var, 2)  % change the dimension
        var = squeeze(y_var(:,varNum, :))';
        var(isnan(var)) = 0;

        var_mean = mean(var, 1);
        var_ci   = 1.96 * std(var, 1)/sqrt(size(var, 1));
        var_lb   = var_mean - var_ci;
        var_ub   = var_mean + var_ci;

        p = fill([x_var, x_var(end:-1:1)], [var_lb, var_ub(end:-1:1)], cmap(varNum, :) - 0.1);
        p.EdgeColor = 'none'; 
        p.FaceColor = cmap(varNum, :) - 0.1;
        p.FaceAlpha = 0.7;

        plot(x_var, var_mean, 'color', cmap(varNum, :), LineWidth=2)
    end
    hold off;
    grid on;
    set(gca, 'XDir', 'reverse')

end