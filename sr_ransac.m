function sr_ransac
% probability that at least one of the random samples of 7 points is inlier
p = 0.99;
totalPtsNum = size(Pos1,1);
[Pos1_norm,Pos2_norm,NormMat1,NormMat2] = sr_dataNorm(Pos1,Pos2);
sampleNum = 9999;
sampleCount = 1;
while sampleCount < sampleNum
    [RndPts1,RndPts2] = sr_randPick7Pts(Pos1_norm,Pos2_norm);
    FdmMat = sr_calFdmMat7CpdPts(RndPts1,RndPts2);
    %tbd
    epsilon = 1 - (numInlier/totalPtsNum);
    % s = 7, s is the minimal sample size to calculate fundamental matrix
    sampleNum = log(1-p)/log(1-(1-epsilon)^7);
    sampleCount = sampleCount + 1;
end