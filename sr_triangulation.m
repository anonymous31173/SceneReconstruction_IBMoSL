function [X1_opt, X2_opt] = sr_triangulation(X1,X2,FdmMat)
% sr_triangulation() find the corrected correspondences X1_opt and X1_opt,
% given measured points X1, X2 and fundamental matrix FdmMat,that minimize
% the geometric error d^2(X1,X1_opt) + d^2(X2,X2_opt) subject to the
% epipolar constrain X2_opt'*FdmMat*X1_opt = 0

% X1 and X2 are measured points in image1 and image2, assumed in 
% homogeneous coordinate, and X1 = [x1, y1, 1], X2 = [x2, y2, 1].

% Reference: Multiple View Geometry in Computer Vision, second edition, page318.

% apply a rigid transformation to each image in order to place both points
% X1 and X2 at the origin (0;0;1)
TrsfMat1 = [1 0 -X1(1); 0 1 -X1(2); 0 0 1];
TrsfMat2 = [1 0 -X2(1); 0 1 -X2(2); 0 0 1];
FdmMat = (TrsfMat2^-1)'*FdmMat*TrsfMat1^-1; % fundamental matrix after translation
Epipole1 = null(FdmMat); % epipole e in image1
Epipole2 = null(FdmMat'); % epipole e' in image2
% normalize e such that e1^2+e2^2=1. same to e'.
Epipole1 = Epipole1/sqrt(Epipole1(1)^2 + Epipole1(2)^2);
Epipole2 = Epipole2/sqrt(Epipole2(1)^2 + Epipole2(2)^2);
% rotate both image to make sure that e has coordinate (1 0 e3) and e' has
% coordinate (1 0 e3')
RotMat1 = [Epipole1(1) Epipole1(2) 0;-Epipole1(2) -Epipole1(1) 0;0 0 1];
RotMat2 = [Epipole2(1) Epipole2(2) 0;-Epipole2(2) -Epipole2(1) 0;0 0 1];
FdmMat = RotMat2*FdmMat*RotMat1'; % fundamental matrix after rotation
% build polynomial equation to find minimum of cost function
f1 = Epipole1(3);f2 = Epipole2(3);
a = FdmMat(2,2);b = FdmMat(2,3);c = FdmMat(3,2);d = FdmMat(3,3);
syms t_sym
eqn = t_sym*((a*t_sym+b)^2+f2^2*(c*t_sym+d)^2)^2 - ...
      (a*d-b*c)*(1+f1^2*t_sym^2)^2*(a*t_sym+b)*(c*t_sym+d) == 0;
t_root = double(vpasolve(eqn,t_sym));
t_realroot = real(t_root);
% also consider cost function has minimum when t tend to be infinity
CostFunction = [t_realroot.^2./(1+f1^2*t_realroot.^2) + ...
               (c*t_realroot+d).^2./((a*t_realroot+b).^2+f2^2*(c*t_realroot+d).^2);...
               1/f1^2 + c^2/(a^2+f2^2*c^2)];
% last term of t_realroot is a representation that t tend to be infinity
t_realroot = [t_realroot;99999];
[~,minOrder] = min(CostFunction);
t = t_realroot(minOrder); % optimal parameter t of 6 root plus infinity
EpipolarLine1 = [t*f1 1 -t]; % epipolar line L in image1 in homogeneous coordinate
EpipolarLine2 = [-f2*(c*t+d);a*t+b;c*t+d]; % epipolar line L' in image2
% closest point(X1_opt) on epipolar line L to the origin(X1)
X1_opt = [-EpipolarLine1(1)*EpipolarLine1(3) -EpipolarLine1(2)*EpipolarLine1(3)...
           EpipolarLine1(1)^2+EpipolarLine1(3)^2];
X2_opt = [-EpipolarLine2(1)*EpipolarLine2(3) -EpipolarLine2(2)*EpipolarLine2(3)...
           EpipolarLine2(1)^2+EpipolarLine2(3)^2];
% transform back to the original coordinates before translation and rotation       
X1_opt = (TrsfMat1^-1 * RotMat1' * X1_opt')';
X2_opt = (TrsfMat2^-1 * RotMat2' * X2_opt')';
% 
X1_opt = X1_opt/X1_opt(3);
X2_opt = X2_opt/X2_opt(3);