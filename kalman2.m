startTime = 5;
stopTime = 29;

video = VideoReader("video/pacman.mp4") ;
video.CurrentTime = startTime;
% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];
params.state_space_bound = [video.Width; video.Height-100]; %1920 1080 
params.bounds = [1, video.Height - 100; 1, video.Width]; % height bounds; width bounds
params.cutoff_dist = 25;


% xhat = [800;850]; % position
% use first frame to initialize start position 
vidFrame = readFrame(video); %read video frame of pacmans, class: uint8
frame = vidFrame(1:params.state_space_bound(2),1:params.state_space_bound(1),:); %vidframe: height, width, dim
dist_image = color_dist(frame, params.pcm_colour); %980x1920
bin_image = dist_image < params.cutoff_dist; %mask
[r, c] = find(bin_image, 1, 'first'); %find the first index where mask is true

vrow = 1; 
vcol = 1; 
delta_t = 1 ; %change in time
xhat = [r;c ; vrow ; vcol]; %4x1  position and velocity

%Initialize matrices
A = [1 0 delta_t 0 ; 0 1 0 delta_t ; 0 0 1 0 ; 0 0 0 1] ; % A 4x4 matrix

%set covariances as variable in case you want to change
cov_pos = eye(2)*1 ;
cov_vel = eye(2)*1 ; 
P = [cov_pos, eye(2)*0; eye(2)*0, cov_vel]; % 4x4, process covariance matrix,assume initially, no covariance between velocity and position
Q_x = eye(2)*2; % process noise position
Q_v = eye(2)*2; % processs noise velocity
Q = [Q_x*delta_t, eye(2)*0 ; eye(2)*0, Q_v*delta_t]; % process noise moition model
C = [1 0 0 0; 0 1 0 0 ;0 0 0 0 ; 0 0 0 0] ; 
R = eye(4)*2; %sensor noise

pos_estimates = [];
pos_groundtruths = [];
pos_errs = [];
meas = [];
meas_predict = [];
time = [];

%% Main LOOP %%
while hasFrame(video) && video.CurrentTime < stopTime
    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8

    frame = vidFrame(1:params.state_space_bound(2),1:params.state_space_bound(1),:); %height, width, so y, x

    % do kalman
%     u = getControl();
    % Prediction
    xhat = A * xhat;
    xhat = min(max(xhat, [1;1;0;0]), [params.state_space_bound(1); params.state_space_bound(2) ;100;100]); %implement clipping a second time
    P = A * P * A' + Q; %add process noise

    % measurement
%     y = doMeasurement(frame, params, xhat(1:2, :));
%     y = get_pacman_center(frame, params);
    x = [xhat(1); xhat(2)];
    y = get_nearest_color_pos(frame, params, x);
%     y = get_min_color_dist_pos(frame,params);

    y = [y ; 0 ; 0] ; % correct dimension

    yhat = C*xhat;

    % Measurement update
    K = P * C' * inv(C * P * C' + R);

    xhat = xhat + K * (y - C * xhat);
    xhat = min(max(xhat, [1;1;0;0]), [params.state_space_bound(1); params.state_space_bound(2) ;100;100]);
    P = P - K * C * P;



%     subplot(2,2,1);
    image(vidFrame,Parent=gca);

    hold on
    plot(xhat(1), xhat(2),'x','Color','red','MarkerSize', 10, 'LineWidth', 5) %plot the estimate 
    drawnow
    hold off

    pos_groundtruth = plot_pacman_center(vidFrame, params);
    x = [xhat(1); xhat(2)]; 


    time = [time video.CurrentTime];
    pos_estimates = [pos_estimates x];
    pos_groundtruths = [pos_groundtruths pos_groundtruth];
    pos_err = norm(x - pos_groundtruth);
    pos_errs = [pos_errs pos_err];
    meas = [meas y];
    meas_predict = [meas_predict yhat];




%     dist = color_dist(frame, params.pcm_colour);
%     bin = dist < params.cutoff_dist;
%     subplot(2,2,2);
%     imshow(mat2gray(bin),Parent=gca);

%     subplot(2,2,3);
%     imshow(vidFrame,Parent=gca)


%     hold on
%     plot(y(2), y(1),'x','Color','magenta') %plot the particles 
%     drawnow
%     hold off



%     pause(0.5)
end

k_plot_stats(time, pos_estimates, pos_groundtruths, pos_errs, meas, meas_predict);

function u = getControl()
%randomly select control
    us = [[1;0], [0;1], [-1;0], [0;-1]]; %down, right, up, left
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

    [D,idx] = bwdist(bin_image); % idx 980x1920, euclidean distance of the binary image

    nearest_nbr_idx = idx(floor(x_hat(1)), floor(x_hat(2)));

    [row,col] = ind2sub([params.state_space_bound(2), params.state_space_bound(1)],nearest_nbr_idx); %width, height // col, row

    pos = [row; col]; %pos = [x, y] = [col, row]
    
end
