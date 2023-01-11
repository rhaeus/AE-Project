% Kalman filter variant 1

clc
clear all
close all

start = 50.5;
% start = 0;

startTime = start;
stopTime = start + 25;

video = VideoReader("video/pacman_full.mp4") ;
video.CurrentTime = startTime;

bottom_cut = video.Height / 15;

% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];
params.state_space_bound = [video.Width; video.Height-bottom_cut]; %1920 1080 
params.bounds = [1, video.Height-bottom_cut ; 1, video.Width]; % height bounds; width bounds
params.cutoff_dist = 25;

video_writer = VideoWriter("k.avi"); %create the video object
video_writer.FrameRate = video.FrameRate;
open(video_writer); %open the file for writing

% xhat = [800;850]; % position
% use first frame to initialize start position 
vidFrame = readFrame(video); %read video frame of pacmans, class: uint8
frame = vidFrame(1:params.state_space_bound(2),1:params.state_space_bound(1),:);
dist_image = color_dist(frame, params.pcm_colour); %980x1920
bin_image = dist_image < params.cutoff_dist;
% imshow(mat2gray(bin_image),Parent=gca);
[r, c] = find(bin_image, 1, 'first');
xhat = [c;r];


% 'velocity' of pacman in pixel per frame
% v = 30;
v = 5;


A = [1 0; 0 1];
B = [v 0; 0 v];
C = [1 0; 0 1];

P = eye(2)*1; % uncertainty of position
G = eye(2); % dim correction for process noise
D = eye(2); % dim correction for measurement noise
R = diag([2 2]); %process noise
Q = diag([1 1]); % measurement noise


pos_estimates = [];
pos_groundtruths = [];
pos_errs = [];
meas = [];
meas_predict = [];
time = [];
Px = [];
Py = [];



%% Main LOOP %%
while hasFrame(video)  && video.CurrentTime < stopTime
    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8

    frame = vidFrame(1:params.state_space_bound(2),1:params.state_space_bound(1),:);

    % do kalman
    u = getControl();
    % Prediction
    xhat = A * xhat + B * u;
    xhat = min(max(xhat, [0;0]), params.state_space_bound);
    P = A * P * A' + G * R * G';

    % measurement
    y = get_pacman_center(frame, params);
%     y = get_nearest_color_pos(frame, params, xhat);
%     y = get_min_color_dist_pos(frame,params);

    % Measurement update
    K = P * C' * inv(C * P * C' + D * Q * D');

    yhat = C*xhat;

    xhat = xhat + K * (y - C * xhat);
    xhat = min(max(xhat, [1;1]), params.state_space_bound);
    P = P - K * C * P;

%     dist = color_dist(frame, params.pcm_colour);
%     bin = dist < params.cutoff_dist;
% %     subplot(2,2,2);
%     imshow(mat2gray(bin),Parent=gca);

%     subplot(2,2,1);
    imshow(vidFrame,Parent=gca);

    hold on
    plot(xhat(1), xhat(2),'x','Color','red','MarkerSize', 15, 'LineWidth', 3) %plot the estimate 
    drawnow
    hold off

    pos_groundtruth = get_pacman_center(vidFrame, params);
%     plot_pacman_center(vidFrame, params);

    F = getframe(gcf);
    writeVideo(video_writer,F); % Write the image to file.


    time = [time video.CurrentTime];
    pos_estimates = [pos_estimates xhat];
    pos_groundtruths = [pos_groundtruths pos_groundtruth];
    pos_err = norm(xhat - pos_groundtruth);
    pos_errs = [pos_errs pos_err];
    meas = [meas y];
    meas_predict = [meas_predict yhat];
    Px = [Px P(1,1)];
    Py = [Py P(2,2)];
    

%     subplot(2,2,3);
%     imshow(vidFrame,Parent=gca)


%     hold on
%     plot(y(1), y(2),'x','Color','cyan','MarkerSize', 15, 'LineWidth', 5) %plot the measurement 
%     drawnow
%     hold off


%     pause(0.5)
end

close(video_writer);

k_plot_stats(time, pos_estimates, pos_groundtruths, pos_errs, meas, meas_predict, Px, Py);

function u = getControl()
    us = [[1;0], [0;1], [-1;0], [0;-1]];
    i = randi([1,4],1);
    u = us(:,i);
end