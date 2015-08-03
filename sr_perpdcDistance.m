function [X1_opt, X2_opt, perpdcDst,flag] = sr_perpdcDistance(X1,X2,FdmMat,method)
% sr_perpdcDistance() calculate perpendicular (square)distance between
% corrected correspondences X1_opt, X2_opt and measured points X1, X2,
% given fundamental matrix FdmMat.

% When method == 'optimal', X1_opt and X2_opt are optimal solution that
% minimize the geometric error d^2(X1,X1_opt) + d^2(X2,X2_opt) subject
% to the epipolar constrain X2_opt'*FdmMat*X1_opt = 0.

% For faster computation, consider using method = 'simple'.
% simple method use d^2(X2,F*X1) + d^2(X1,F'*X2) as distance, which would
% be greater than optimal distance d^2(X1,X1_opt) + d^2(X2,X2_opt), but
% need much less computation time.

% X1 and X2 are measured points in image1 and image2, assumed in 
% homogeneous coordinate, and X1 = [x1, y1, 1], X2 = [x2, y2, 1].

% Reference: Multiple View Geometry in Computer Vision, second edition, page318.
flag = 0;
% to homogeneous coordinate and convert them to column vector
X1 = [X1,1]'; X2 = [X2,1]';

if ~exist('method')
    method = 'optimal';
end
if strcmp(method,'simple')
    L2_estm = FdmMat*X1;
    L1_estm = FdmMat'*X2;
    % X1_opt is the closest point on L1_estm to X1
    X1_opt = [L1_estm(1)*L1_estm(2)*X1(2)-L1_estm(2)^2*X1(1)+L1_estm(1)*L1_estm(3);...
              L1_estm(2)*L1_estm(3)-L1_estm(1)^2*X1(2)+L1_estm(1)*L1_estm(2)*X1(1);...
             -L1_estm(1)^2-L1_estm(2)^2];
    X2_opt = [L2_estm(1)*L2_estm(2)*X2(2)-L2_estm(2)^2*X2(1)+L2_estm(1)*L2_estm(3);...
              L2_estm(2)*L2_estm(3)-L2_estm(1)^2*X2(2)+L2_estm(1)*L2_estm(2)*X2(1);...
             -L2_estm(1)^2-L2_estm(2)^2];
    % back to inhomogeneous coordinate
    X1_opt = X1_opt/(X1_opt(3));
    X2_opt = X2_opt/(X2_opt(3));
    % perpdcDst is sum of distance from points to lines rather points to
    % points.
    perpdcDst = (L1_estm'*X1)^2/(L1_estm(1)^2+L1_estm(2)^2) + ...
                (L2_estm'*X2)^2/(L2_estm(1)^2+L2_estm(2)^2);
    return
end
% apply a rigid transformation to each image in order to place both points
% X1 and X2 at the origin (0;0;1)
TrsfMat1 = [1 0 -X1(1); 0 1 -X1(2); 0 0 1];
TrsfMat2 = [1 0 -X2(1); 0 1 -X2(2); 0 0 1];
FdmMat = (TrsfMat2^-1)'*FdmMat*TrsfMat1^-1; % fundamental matrix after translation
% fundamental matrix should always be singular, but in case of noise this
% constrain may be invalid, so use svd to singularize FdmMat.
if rank(FdmMat) == 3
    FdmMat = singlrz(FdmMat);
end
Epipole1 = null(FdmMat); % epipole e in image1
Epipole2 = null(FdmMat'); % epipole e' in image2
if numel(Epipole1) == 0 || numel(Epipole2) == 0
    X1_opt = 0;X2_opt = 0;perpdcDst = 0;flag = 1;
    return
end
% normalize e such that e1^2+e2^2=1. same to e'.
Epipole1 = Epipole1/sqrt(Epipole1(1)^2 + Epipole1(2)^2);
Epipole2 = Epipole2/sqrt(Epipole2(1)^2 + Epipole2(2)^2);
% rotate both image to make sure that e has coordinate (1 0 e3) and e' has
% coordinate (1 0 e3')
RotMat1 = [Epipole1(1) Epipole1(2) 0;-Epipole1(2) Epipole1(1) 0;0 0 1];
RotMat2 = [Epipole2(1) Epipole2(2) 0;-Epipole2(2) Epipole2(1) 0;0 0 1];
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
X1_opt = [-EpipolarLine1(1)*EpipolarLine1(3); -EpipolarLine1(2)*EpipolarLine1(3); ...
           EpipolarLine1(1)^2+EpipolarLine1(2)^2];
X2_opt = [-EpipolarLine2(1)*EpipolarLine2(3); -EpipolarLine2(2)*EpipolarLine2(3); ...
           EpipolarLine2(1)^2+EpipolarLine2(2)^2];
% transform back to the original coordinates before translation and rotation       
X1_opt = TrsfMat1^-1 * RotMat1' * X1_opt;
X2_opt = TrsfMat2^-1 * RotMat2' * X2_opt;
% back to inhomogeneous coordinate for distance calculation
X1_opt = X1_opt./X1_opt(3);
X2_opt = X2_opt./X2_opt(3);

perpdcDst = (X1(1) - X1_opt(1))^2 + (X1(2) - X1_opt(2))^2 + ...
            (X2(1) - X2_opt(1))^2 + (X2(2) - X2_opt(2))^2;
end

function FdmMat_singlrz = singlrz(FdmMat)
% singlrz() enforce fundamental matrix FdmMat is singular and hence
% rank(FdmMat) == 2.
[U,D,V] = svd(FdmMat);
D(3,3) = 0;
FdmMat_singlrz = U*D*V';
end