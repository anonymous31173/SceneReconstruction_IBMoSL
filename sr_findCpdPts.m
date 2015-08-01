function [Pos1,Pos2] = sr_findCpdPts(img1,img2,scoreThreshold,scale,displayFlag)
% sr_findCpdPts() finds corresponding points between two images img1 and img2
% using vl_SIFT algorithm and return their coordinates in orginal images.

% Pos1 is a n*2 matrix with first and second column is x and y coordinates,
% respectively, of corresponding points in Img1; Same to Pos2 as
% coordinates of Img2.

% scoreThreshold is a threshold to filtrate matched points, the smaller
% value it is, the better match result and less corresponding points.

% scale is used to resize original image for compute efficiency. the
% greater it is, the faster to get result.

% set displayFlag 1 to display image pair with corresponding.
if ~exist('scoreThreshold')
    scoreThreshold = 0.5;
end
% ensure the image to be implemented has width between [800 1600] pixels
[height,width,channel] = size(img1);
if ~exist('scale')
    scale = 2.^(ceil(log(ceil(width/1600))/log(2)));
end
if ~exist('displayFlag')
    displayFlag = 0;
end
if prod(size(img1)==size(img2)) == 0
    error('The dimension of two input images must be identical.')
end
% the input of vl_sift must be gray-scale image in single precision
if channel > 1
    Img1 = rgb2gray(img1);
    Img2 = rgb2gray(img2);
end
if ~isa(Img1,'single')
    Img1 = single(Img1);
    Img2 = single(Img2);
end
Img1 = imresize(Img1,[height/scale width/scale]);
Img2 = imresize(Img2,[height/scale width/scale]);
% find feature points location and descriptors using vl_feat toolbox
[Frame1,Dscpt1] = vl_sift(Img1);
[Frame2,Dscpt2] = vl_sift(Img2);
[Matches,Scores] = vl_ubcmatch(Dscpt1,Dscpt2,3); % match features between two images
BadMatches = Scores>scoreThreshold*max(Scores); % discard matches with poor scores
Matches(:,BadMatches) = [];
if length(Matches)<100
    error('Too few matches points, increase scoreThreshold value or use another image pair.')
end
MatchedFrame1 = Frame1(:,Matches(1,:));
MatchedFrame2 = Frame2(:,Matches(2,:));
Pos1 = scale*MatchedFrame1(1:2,:)'; % get matched points location
Pos2 = scale*MatchedFrame2(1:2,:)';
% for visualization
if displayFlag
    ImgCombine = [img1,img2];
    Pos2_temp = Pos2(:,1)+width;
    figure(1);imshow(ImgCombine);
    X = [Pos1(:,1)';Pos2_temp'];
    Y = [Pos1(:,2)';Pos2(:,2)'];
    line(X,Y,'color',[.2 .8 .2]);
    figure(2);
    imshow(img1);hold on;scatter(Pos1(:,1),Pos1(:,2),[],[.1 .9 .1]);
    figure(3);
    imshow(img2);hold on;scatter(Pos2(:,1),Pos2(:,2),[],[.1 .9 .1]);
end

