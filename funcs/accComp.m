function [dist, time] = accComp(traj, startPt, endPt, timeVec)
    for i = 1: size(traj, 2)
        shoDist(i) = abs( (endPt(1)-startPt(1)) * (startPt(2)-traj(2, i)) - (startPt(1)-traj(1, i)) * (endPt(2)-startPt(2)) ) / sqrt((endPt(1)-startPt(1))^2 + (endPt(2)-startPt(2))^2);
    end

    dist = sqrt(sum(shoDist.^2)/length(shoDist));
    time = size(timeVec, 2);
end