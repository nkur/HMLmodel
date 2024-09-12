function [smoothedVec] = smootherFun(origVec, windowLen)
    for REiter = 1: length(origVec)-windowLen
       tmp1 = 0;
       for AViter = REiter : (REiter-1)+windowLen
           tmp1 = tmp1 +  origVec(AViter);
       end
       smoothedVec(REiter) = tmp1/windowLen;
    end
end