function pf_plot_stats(time, weight_avgs, pos_estimates, pos_groundtruths, pos_errs, particle_stddevs)

    err = pos_groundtruths - pos_estimates;
    maex = mean(abs(err(1,:)));
    maey = mean(abs(err(2,:)));
    fprintf('mean absolute err x %0.4f\n', maex);
    fprintf('mean absolute err y %0.4f\n', maey);



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