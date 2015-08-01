function [Pos1_norm,Pos2_norm,NormMat1,NormMat2] = sr_dataNorm(Pos1,Pos2)
% sr_dataNorm() normalize points position so that they obey properties
% follow: 1) their centroid is at the origin; 2) the average distance from
% the origin is equal to sqrt(2).

% Reference: Multiple View Geometry in Computer Vision, second edition, page107.

Pos1 = [Pos1,ones(size(Pos1,1),1)]; % to homogeneous coordinate
Pos2 = [Pos2,ones(size(Pos2,1),1)];
x_offset1 = mean(Pos1(:,1)); % translation of Pos1
y_offset1 = mean(Pos1(:,2));
x_offset2 = mean(Pos2(:,1)); % translation of Pos2
y_offset2 = mean(Pos2(:,2));
TrsltMat1 = [1,0,-x_offset1;0,1,-y_offset1;0,0,1];
TrsltMat2 = [1,0,-x_offset2;0,1,-y_offset2;0,0,1];
% square of average distance of points in Pos1 from the origin
scale1 = mean((Pos1(:,1)-x_offset1).^2+(Pos1(:,2)-y_offset1).^2);
scale2 = mean((Pos2(:,1)-x_offset2).^2+(Pos2(:,2)-y_offset2).^2);
ScaleMat1 = [sqrt(2/scale1),0,0;0,sqrt(2/scale1),0;0,0,1];
ScaleMat2 = [sqrt(2/scale2),0,0;0,sqrt(2/scale2),0;0,0,1];
NormMat1 = ScaleMat1*TrsltMat1;
NormMat2 = ScaleMat2*TrsltMat2;
Pos1_norm = (NormMat1*Pos1')';
Pos2_norm = (NormMat2*Pos2')';