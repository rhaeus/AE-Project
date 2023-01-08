% histogram 1080x1920
% S.W 1xM
% S.X 2x1000
function [S_bar, w_avg] = pf_weight(S, params, histogram)
    S_bar = S;
    x_coords = ceil(S_bar.X(1, :)) ; % 1x1000
    y_coords = ceil(S_bar.X(2, :)) ; % 1x1000


    for m=1:params.M
        x = x_coords(1,m);
        y = y_coords(1,m);

        if (x > params.bounds(2,2) || x <= params.bounds(2,1) || y > params.bounds(1,2) || y <= params.bounds(1,1))
            S_bar.W(m) = 0;
        else
            S_bar.W(m) = histogram(y, x);
        end
    end

    w_avg = mean(S_bar.W);

    % normalize weights
    S_bar.W = S_bar.W ./ sum(S_bar.W) ;
end