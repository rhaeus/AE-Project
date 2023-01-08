function positions = find_pacman_positions(frame, params)
    frame = frame(params.bounds(1,1):params.bounds(1,2), params.bounds(2,1):params.bounds(2,2),:);
    dist_image = color_dist(frame, params.pcm_colour); %980x1920
    bin_image = dist_image < params.cutoff_dist;

    ps = find(bin_image);

    n = size(ps,1);
    %     pos = ind2sub([params.bounds(1,2), params.bounds(2,2)], pos);
    %     pos = ind2sub(size(bin_image), pos);
    
    positions = zeros([2 n]);
    for i = 1 : n
        [c,r] = ind2sub(size(bin_image), ps(i));
        positions(1,i) = r;
        positions(2,i) = c;
    end
end