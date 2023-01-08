function plot_particles(S)
    hold on
    plot(S.X(1,:),S.X(2,:),'.','Color','green') %plot the particles 
    drawnow
    hold off

    hold on
    plot(mean(S.X(1,:)),mean(S.X(2,:)),'x','Color','red','MarkerSize', 10, 'LineWidth', 5) %plot the particles avg 
    drawnow
    hold off
end