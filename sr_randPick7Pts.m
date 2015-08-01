function [RndPts1,RndPts2] = sr_randPick7Pts(Pos1,Pos2)
% sr_randPick7Pts() random pick 7 corresponding points over whole image range
% from all putative correspoing points.
% 7 corresponding points should have a good spatial distribution over the
% image to improve robustness, so randomly select 20 points from whole
% image range and find 20 corresponding points closest to them.

% Pos1 and Pos2 are all putative points in image1 and image2.

% RndPts1 is a 7*2 matrix with columns containing coordinates of points in
% image1.

Pos = [Pos1,Pos2];
width = max(max(Pos(:,1)),max(Pos(:,3)));
height = max(max(Pos(:,2)),max(Pos(:,4)));
rand_x = randi([1 round(width)],20,1);
rand_y = randi([1 round(height)],20,1);
for i = 1:20
    [~,order(i)] = min(abs((Pos(:,1)-rand_x(i)) + (Pos(:,2)-rand_y(i))));
end
order = unique(order);
if length(order) < 7
    error('Can not pick 7 points, try one more time.');
end
order = order(1:7);
RndPts1 = Pos(order,1:2);
RndPts2 = Pos(order,3:4);
