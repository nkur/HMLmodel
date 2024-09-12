function [maxDistance, maxDistanceIndex] = maxPerpendicularDist(x, startPt, endPt)

    distances = point_to_line(x, startPt, endPt);
    [maxDistance, maxDistanceIndex] = max(distances);

    function d = point_to_line(pt, v1, v2)
        % pt should be nx3
        % v1 and v2 are vertices on the line (each 1x3)
        % d is a nx1 vector with the orthogonal distances

        %prepare inputs
        pt = pt'; %2xn to nx2
        v1=v1(:)';%force 1x3 or 1x2
        v2=v2(:)';%force 1x3 or 1x2
        if length(v1)==2,v1(3)=0;  end %extend 1x2 to 1x3 if needed
        if length(v2)==2,v2(3)=0;  end %extend 1x2 to 1x3 if needed
        if size(pt,2)==2,pt(1,3)=0;end %extend nx2 to nx3 if needed
        v1_ = repmat(v1,size(pt,1),1);
        v2_ = repmat(v2,size(pt,1),1);

        % Caluculating perpendicular distances
        a = v1_ - v2_;
        b = pt - v2_;
        d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
    end

end