video = VideoReader("video/pacman.mp4") ;
video.CurrentTime = 5;

% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];
params.state_space_bound = [video.Width; video.Height-100]; %1920 1080 
params.bounds = [1, video.Height - 100; 1, video.Width]; % height bounds; width bounds
params.cutoff_dist = 25;


% xhat = [800;850]; % position
% use first frame to initialize start position 
vidFrame = readFrame(video); %read video frame of pacmans, class: uint8
frame = vidFrame(1:params.state_space_bound(2),1:params.state_space_bound(1),:);
dist_image = color_dist(frame, params.pcm_colour); %980x1920
bin_image = dist_image < params.cutoff_dist;
[r, c] = find(bin_image, 1, 'first');
xhat = [r;c];

% 'velocity' of pacman in pixel per frame
v = 30;


A = [1 0; 0 1];
B = [v 0; 0 v];
C = [1 0; 0 1];

P = eye(2)*1; % uncertainty of position
G = eye(2); % dim correction for process noise
D = eye(2); % dim correction for measurement noise
R = diag([10 10]); %process noise
Q = diag([4 4]); % measurement noise





%% Main LOOP %%
while hasFrame(video) 
    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8

    frame = vidFrame(1:params.state_space_bound(2),1:params.state_space_bound(1),:);

    % do kalman
    u = getControl();
    % Prediction
    xhat = A * xhat + B * u;
    xhat = min(max(xhat, [0;0]), params.state_space_bound);
    P = A * P * A' + G * R * G';

    % measurement
    y = doMeasurement(frame, params, xhat);

    % Measurement update
    K = P * C' * inv(C * P * C' + D * Q * D');

    xhat = xhat + K * (y - C * xhat);
    xhat = min(max(xhat, [0;0]), params.state_space_bound);
    P = P - K * C * P;



%     subplot(2,2,1);
    image(vidFrame,Parent=gca);

    hold on
    plot(xhat(1), xhat(2),'x','Color','red','MarkerSize', 10, 'LineWidth', 5) %plot the estimate 
    drawnow
    hold off

    plot_pacman_center(vidFrame, params);



%     dist = color_dist(frame, params.pcm_colour);
%     bin = dist < params.cutoff_dist;
%     subplot(2,2,2);
%     imshow(mat2gray(bin),Parent=gca);

%     subplot(2,2,3);
%     imshow(vidFrame,Parent=gca)


%     hold on
%     plot(y(1), y(2),'x','Color','red') %plot the measurement 
%     drawnow
%     hold off

    



    pause(0.5)
end

function u = getControl()
    us = [[1;0], [0;1], [-1;0], [0;-1]];
    i = randi([1,4],1);
    u = us(:,i);
end

function pos = doMeasurement(frame, params, x_hat)
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