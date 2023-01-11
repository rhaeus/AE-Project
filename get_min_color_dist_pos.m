% return a pixel position with minimum color distance to target color

function pos = get_min_color_dist_pos(frame, params)
    frame = frame(params.bounds(1,1):params.bounds(1,2), params.bounds(2,1):params.bounds(2,2),:);
    dist_mat = color_dist(frame, params.pcm_colour);
    [~, i] = min(dist_mat, [], "all");
    
    [c,r] = ind2sub(size(dist_mat),i);
    pos = [r;c];
end