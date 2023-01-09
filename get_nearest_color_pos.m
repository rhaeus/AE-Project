
function pos = get_nearest_color_pos(frame, params, x_hat)
    % calc distance of each pixel in RGB space to target color
    % mark each pixel with distance greater than cutoff as 0, others as 1
    % nearest neighbor: use position of nearest pixel with color marked
    % with 1 as measurement of pacman position

    dist_image = color_dist(frame, params.pcm_colour); %980x1920

    bin_image = dist_image < params.cutoff_dist;

    [D,idx] = bwdist(bin_image); % idx 980x1920

    nearest_nbr_idx = idx(floor(x_hat(2)), floor(x_hat(1)));

    [row,col] = ind2sub([params.state_space_bound(2), params.state_space_bound(1)],nearest_nbr_idx);

    pos = [col; row];
    
end