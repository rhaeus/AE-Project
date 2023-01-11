% return the center position of pacman based on color
% mean of all pixel positions considered belonging to pacman

function pos = get_pacman_center(frame, params)
    positions = find_pacman_positions(frame, params);
    pos = mean(positions,2);
end