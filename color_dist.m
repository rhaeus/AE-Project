% calculate euclidean distance between color of each pixel and target color
% treating color in RGB space

function dist_mat = color_dist(frame, colour)
    [H, W, D] = size(frame); %height, width, dimension of video frame matrix
    RGB_matrix = double(reshape(frame,[H*W, D])); % create matrix with R G B values listed in separate columns
    dist_mat = pdist2(RGB_matrix, colour, "euclidean");
    dist_mat = reshape(dist_mat, H, W);
end