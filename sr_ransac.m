function [FdmMat_final,InlierOrder_final,sampleCount] = sr_ransac(Pos1,Pos2,method)
if ~exist('method')
    method = 'optimal';
end
% probability that at least one of the random samples of 7 points is inlier
p = 0.95;
% correspondece with perpendicular distance less than inlierThreshold
% will be regarded as inlier
inlierThreshold = 8;
totalCrpPtsNum = size(Pos1,1);

% [Pos1_norm,Pos2_norm,NormMat1,NormMat2] = sr_dataNorm(Pos1,Pos2);
Pos1_norm = Pos1;Pos2_norm = Pos2;

% the fundamental matrix hasing the max number of inlier is chosen
maxInlierNum = 0;
minSTD = 99999;
sampleNum = 99999;
sampleCount = 0;
while sampleCount < sampleNum
    % randomly pick 7 point-pairs from all putative points
    [RndPts1,RndPts2] = sr_randPick7CrpPts(Pos1_norm,Pos2_norm);
    % calculate fundamental matrix from these 7 pairs
    FdmMat = sr_calFdmMat7CpdPts(RndPts1,RndPts2);
    % calculate the perpendicular distance for each putative correspondence
    if size(FdmMat,3) == 1 % case that only one real solution for FdmMat
        for iCrp = 1:totalCrpPtsNum
            X1 = Pos1_norm(iCrp,:); X2 = Pos2_norm(iCrp,:);
            [X1_opt, X2_opt, perpdcDst(iCrp,1),flag] = sr_perpdcDistance(X1,X2,FdmMat,method);
        end
        if flag
            continue
        end
        % the order of all putative points that has distance smaller than threshold
        InlierOrder = perpdcDst < inlierThreshold;
        inlierNum = nnz(InlierOrder); % the number of inlier
        % choose the FdmMat with most inlier.
        if inlierNum > maxInlierNum
            FdmMat_final = FdmMat;
            InlierOrder_final = InlierOrder;
            maxInlierNum = inlierNum;
        % in case of same number of inlier, choose the FdmMat that has
        % the lowest standard deviation of inlier.
        elseif inlierNum == maxInlierNum && std(perpdcDst(InlierOrder)) < minSTD
            FdmMat_final = FdmMat;
            InlierOrder_final = InlierOrder;
            minSTD = std(perpdcDst(InlierOrder));
        end
        epsilon = 1 - (inlierNum/totalCrpPtsNum);
    elseif size(FdmMat,3) == 3 % case that three real solutions for FdmMat
        for iFdmMat = 1:3
            FdmMat_temp = FdmMat(:,:,iFdmMat);
            for iCrp = 1:totalCrpPtsNum
                X1 = Pos1_norm(iCrp,:); X2 = Pos2_norm(iCrp,:);
                [X1_opt, X2_opt, perpdcDst(iCrp,1),flag] = sr_perpdcDistance(X1,X2,FdmMat_temp,method);
            end
            if flag
                continue
            end
            InlierOrder(:,iFdmMat) = perpdcDst < inlierThreshold;
            inlierNum(iFdmMat) = nnz(InlierOrder(:,iFdmMat));
        end
        [inlierNum,maxFdmMatOrder] = max(inlierNum);
        if inlierNum > maxInlierNum
            FdmMat_final = FdmMat(:,:,maxFdmMatOrder);
            InlierOrder_final = InlierOrder(:,maxFdmMatOrder);
            maxInlierNum = inlierNum;
        elseif inlierNum == maxInlierNum && std(perpdcDst(InlierOrder(:,maxFdmMatOrder))) < minSTD
            FdmMat_final = FdmMat(:,:,maxFdmMatOrder);
            InlierOrder_final = InlierOrder(:,maxFdmMatOrder);
            minSTD = std(perpdcDst(InlierOrder(:,maxFdmMatOrder)));
        end
        epsilon = 1 - (inlierNum/totalCrpPtsNum);
    end
    %tbd   
    % s = 7, s is the minimal sample size to calculate fundamental matrix
    sampleNum_temp = log(1-p)/log(1-(1-epsilon)^7);
    if sampleNum_temp < sampleNum
        sampleNum = sampleNum_temp;
    end
    sampleCount = sampleCount + 1;
end