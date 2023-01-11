%Kalman filter variant 2

clc
clear all
close all

start = 50.5;
video = VideoReader("video/pacman_full.mp4") ;
params.Sigma_R = diag([20 20]);

% start = 5;
% video = VideoReader("video/pacman.mp4") ;
% params.Sigma_R = diag([200 200]);

startTime = start;
stopTime = start + 25;
video.CurrentTime = startTime;

bottom_cut = video.Height / 15;

% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];
params.state_space_bound = [video.Width; video.Height-bottom_cut]; %1920 1080 
params.bounds = [1, video.Height - bottom_cut; 1, video.Width]; % height bounds; width bounds
params.cutoff_dist = 25;

video_writer = VideoWriter("k2.avi"); %create the video object
video_writer.FrameRate = video.FrameRate;
open(video_writer); %open the file for writing


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
R = eye(4)*1; %sensor noise

pos_estimates = [];
pos_groundtruths = [];
pos_errs = [];
meas = [];
meas_predict = [];
time = [];
Px = [];
Py = [];

%% Main LOOP %%
while hasFrame(video) && video.CurrentTime < stopTime
    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8

    frame = vidFrame(1:params.state_space_bound(2),1:params.state_space_bound(1),:); %height, width, so y, x

    % do kalman
%     u = getControl();
    % Prediction
    xhat = A * xhat;
    xhat = min(max(xhat, [1;1;-100;-100]), [params.state_space_bound(1); params.state_space_bound(2) ;100;100]); %implement clipping a second time
    P = A * P * A' + Q; %add process noise

    % measurement
    y = get_pacman_center(frame, params);
%     x = [xhat(1); xhat(2)];
%     y = get_nearest_color_pos(frame, params, x);
%     y = get_min_color_dist_pos(frame,params);

    y = [y ; 0 ; 0] ; % correct dimension

    yhat = C*xhat;

    % Measurement update
    K = P * C' * inv(C * P * C' + R);

    xhat = xhat + K * (y - C * xhat);
    xhat = min(max(xhat, [1;1;-100;-100]), [params.state_space_bound(1); params.state_space_bound(2) ;100;100]);
    P = P - K * C * P;



%     subplot(2,2,1);
    imshow(vidFrame,Parent=gca);

    hold on
    plot(xhat(1), xhat(2),'x','Color','red','MarkerSize', 10, 'LineWidth', 3) %plot the estimate 
    drawnow
    hold off

    pos_groundtruth = get_pacman_center(vidFrame, params);
%     plot_pacman_center(vidFrame, params);
    x = [xhat(1); xhat(2)]; 

    F = getframe(gcf);
    writeVideo(video_writer,F); % Write the image to file.


    time = [time video.CurrentTime];
    pos_estimates = [pos_estimates x];
    pos_groundtruths = [pos_groundtruths pos_groundtruth];
    pos_err = norm(x - pos_groundtruth);
    pos_errs = [pos_errs pos_err];
    meas = [meas y];
    meas_predict = [meas_predict yhat];
    Px = [Px P(1,1)];
    Py = [Py P(2,2)];




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



%     pause(0.05)
end
close(video_writer);
k_plot_stats(time, pos_estimates, pos_groundtruths, pos_errs, meas, meas_predict, Px, Py);


