function pf_plot_stats(time, weight_avgs, pos_estimates, pos_groundtruths, pos_errs, particle_stddevs)
%     MSE = (1/n) * Σ(y_i - ŷ_i)^2
%     mse_x = mean(sum(pos_groundtruths(1,:) - pos_estimates(1,:))^2);
%     mse_y = mean(sum(pos_groundtruths(2,:) - pos_estimates(2,:))^2);
% % 
% % %     variance = (1/(n-1)) * Σ(y_i - mean(y))^2
% % %     normalized MSE = MSE / variance
% % %     var_x = var(pos_groundtruths(1,:));
% % %     var_y = var(pos_groundtruths(2,:));
% % 
% %     var_x = var(pos_errs(1,:));
% %     var_y = var(pos_errs(2,:));
% % 
% %     n_mse_x = mse_x / var_x;
% %     n_mse_y = mse_y / var_y;
%     n_mse_x = 0;
%     n_mse_y = 0;
% 
%     fprintf('The normalized mean-squared error for x is %0.2f, normalized is %0.2f\n', mse_x, n_mse_x);
%     fprintf('The normalized mean-squared error for y is %0.2f, normalized is %0.2f\n', mse_y, n_mse_y)

    err = pos_groundtruths - pos_estimates;
    maex = mean(abs(err(1,:)));
    maey = mean(abs(err(2,:)));
    fprintf('mean absolute err x %0.4f\n', maex);
    fprintf('mean absolute err y %0.4f\n', maey);


%     mse = immse(pos_estimates(1,:),pos_groundtruths(1:x));
%     var = var(pos_errs)
%     fprintf('The mean-squared error is %0.4f\n', mse);
% 
%     i = find(pos_errs < 30 ,1,'first');
%     mse = immse(pos_estimates(i:end),pos_groundtruths(i:end));
%     fprintf('The mean-squared error after convergence is %0.4f\n', mse);

    figure('Name','Weight Average');
    plot(time, weight_avgs);
    title('Weight average');
    xlabel('Time (s)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','Postion');
    subplot(2,1,1)
    plot(time, pos_estimates(1,:));
    hold on
    plot(time, pos_groundtruths(1,:));
    hold off
    
    title('X position');
    legend('estimate x','groundtruth x')
    xlabel('Time (s)');

    subplot(2,1,2)
    plot(time, pos_estimates(2,:));
    hold on
    plot(time, pos_groundtruths(2,:));
    hold off

    title('Y position');
    legend('estimate y','groundtruth y')
    xlabel('Time (s)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','Postion Error');
    subplot(2,1,1)
    plot(time, pos_groundtruths(1,:) - pos_estimates(1,:));
    
    title('x error');
    xlabel('Time (s)');

    subplot(2,1,2)
    plot(time, pos_groundtruths(2,:) - pos_estimates(2,:));

    title('y error');
    xlabel('Time (s)');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','Postion error distance');
    plot(time, pos_errs);

    title('Position error distance');
    xlabel('Time (s)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','Particle Standard Deviation');
    subplot(2,1,1)
    plot(time, particle_stddevs(1,:));
    
    title('X Standard Deviation');
    xlabel('Time (s)');

    subplot(2,1,2)
    plot(time, particle_stddevs(2,:))

    title('Y Standard Deviation');
    xlabel('Time (s)');



end