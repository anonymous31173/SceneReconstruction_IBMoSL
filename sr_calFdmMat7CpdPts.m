function FdmMat = sr_calFdmMat7CpdPts(CpdPts1,CpdPts2)
% sr_calFdmMat7CpdPts() calculate fundamental matrix from 7 correspinding
% points.

% Because of degeneracy, there will be one or three real solution, so
% FdmMat is a 3*3*1 or 3*3*3 matrix.

% CpdPts1 is a 7*2 matrix with columns cotaining x and y coordinates of
% putative corresponding points in image1. CpdPts2 has the form in a
% similar way.

% Reference: Multiple View Geometry in Computer Vision, second edition, page281.
A = [CpdPts2(:,1).*CpdPts1(:,1) , CpdPts2(:,1).*CpdPts1(:,2) , CpdPts2(:,1),...
     CpdPts2(:,2).*CpdPts1(:,1) , CpdPts2(:,2).*CpdPts1(:,2) , CpdPts2(:,2),...
     CpdPts1(:,1) , CpdPts1(:,2) , ones(7,1)];
 if  rank(A) ~= 7
     error('Rank of A is not equal to 7.')
 end
 NullSpace = null(A);
 F1 = reshape(NullSpace(:,1),3,3)';
 F2 = reshape(NullSpace(:,2),3,3)';
 % fundamental matrix has property that det(F) = 0.
 % solve a cubic polynomial equation and find the real root alpha.
 syms alpha_sym;
 F = alpha_sym*F1 + (1-alpha_sym)*F2;
 eqn = det(F)==0;
 alpha_root = double(vpasolve(eqn,alpha_sym));
 if length(alpha_root) ~= 3
     error('Can not find 3 roots for cubic polynomial equation.')
 end
 for iRoot = 1:3
     if isreal(alpha_root(iRoot))
         alpha = alpha_root(iRoot);
         FdmMat(:,:,iRoot) = alpha*F1 + (1-alpha)*F2;
     end
 end
 