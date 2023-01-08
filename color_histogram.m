function histogram = color_histogram(frame, colour)
    [H, W, D] = size(frame); %height, width, dimension of video frame matrix
    RGB_matrix = double(reshape(frame,[H*W, D])); % create matrix with R G B values listed in separate columns
    histogram = pdist2(RGB_matrix, colour, "euclidean");
    histogram = 1./histogram ; 
    % histogram = histogram / sum(sum(histogram)) ;
    
    histogram = reshape(histogram, H, W);
    % histogram_size = size(histogram);
    % frame_proc = reshape(histogram, H, W);
    % disp("histogram filter successful")

end