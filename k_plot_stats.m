function k_plot_stats(time, pos_estimates, pos_groundtruths, pos_errs, meas, meas_predict, Px, Py)
%     mse = immse(pos_estimates,pos_groundtruths);
%     fprintf('The mean-squared error is %0.4f\n', mse);
% 
%     i = find(pos_errs < 30 ,1,'first');
%     mse = immse(pos_estimates(i:end),pos_groundtruths(i:end));
%     fprintf('The mean-squared error after convergence is %0.4f\n', mse);

    err = pos_groundtruths - pos_estimates;

    maex = mean(abs(err(1,:)));
    maey = mean(abs(err(2,:)));
    fprintf('mean absolute err x %0.4f\n', maex);
    fprintf('mean absolute err y %0.4f\n', maey);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','Position');
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','Postion error dist');
    plot(time, pos_errs);

    title('Position error dist');
    xlabel('Time (s)');

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','Innovation');
    innovation = (meas - meas_predict);
    plot(time, innovation);

    title('Innovation');
    xlabel('Time (s)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','Covariance');
    subplot(2,1,1)
    plot(time, Px);
    
    title('Variance x');
    xlabel('Time (s)');

    subplot(2,1,2)
    plot(time, Py);

    title('Variance y');
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

end