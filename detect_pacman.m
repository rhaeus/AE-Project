import clib.opencv.*;
import vision.opencv.*;

videoSample = VideoReader("video/pacman.mp4");
videoSample.CurrentTime = 2;
display(videoSample.FrameRate)

% frame = imread("video/frame00.png");
% image(frame)


history = 100;
threshold = 400;
shadow = false;

cvPtr = cv.createBackgroundSubtractorKNN(history,threshold,shadow);
kNNBase = util.getBasePtr(cvPtr);
% 
foregroundmask = zeros(videoSample.Height,videoSample.Width,videoSample.NumFrames);
while hasFrame(videoSample)
    frame = readFrame(videoSample);

    subplot(2,2,1);
    image(frame,Parent=gca);

    [inMat,imgInput] = util.createMat(frame);
    [outMat,outImg] = util.createMat();
    kNNBase.apply(imgInput,outImg);
    
    foregroundmask = util.getImage(outImg);
%     image(foregroundmask);
        
    foregroundmask = rescale(foregroundmask);
    foregroundmask = cast(foregroundmask,"like",frame);

%     image(foregroundmask);
    
    foreground(:,:,1) = frame(:,:,1).*foregroundmask;
    foreground(:,:,2) = frame(:,:,2).*foregroundmask;
    foreground(:,:,3) = frame(:,:,3).*foregroundmask;

    subplot(2,2,2);
    image(foreground,Parent=gca);

%     [inputMat,inputArray] = util.createMat(foreground);
%     [outMat,contours] = util.createMat();
%     d = edge(outImg, 'canny', .6); 

%     cont = cv.findContours(frame );

    pause(0.5);
end