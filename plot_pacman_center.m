function pacman_center = plot_pacman_center(frame, params)
    pacman_center = get_pacman_center(frame, params);
    hold on
    plot(pacman_center(1), pacman_center(2),'x','Color','magenta','MarkerSize', 10, 'LineWidth', 5) %plot the pacman center 
    drawnow
    hold off
end