% return pos (2xn) with n = number of pixel with pacman color

function pos = get_pacman_center(frame, params)
    positions = find_pacman_positions(frame, params);
    pos = mean(positions,2);
end