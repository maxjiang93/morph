function M = SphRotation(n,theta)
% Find matrix for unit sphere rotation given sample vector transform
% in cartesian coordinates
% Input: n: normal vector fo rotation(3d)
%        theta: rotation angle in radians
% Output: Linear Transform Matrix

% By Chiyu 'Max' Jiang, May 7, 2016

c = eye(3);
M = zeros(3,3);
M(:,1) = Rot(n,theta,c(:,1));
M(:,2) = Rot(n,theta,c(:,2));
M(:,3) = Rot(n,theta,c(:,3));
end



function d = Rot(n,theta,c)

if norm(c)~=1
    n = n/norm(n);
end

cp = (dot(c,n))*n;
ct = c-cp;
if norm(ct)~=0
    f = norm(ct)^2*cos(theta);
    e = norm(ct)^2*sin(theta)*n;
    rhs1 = [e(1)+e(2);e(3);f];
    LHS1 = [ ct(3),-ct(3), ct(2)-ct(1);...
        -ct(2), ct(1), 0        ;...
        ct(1), ct(2), ct(3)      ];
    rhs2 = [e(1)+e(3);e(2);f];
    LHS2 = [ -ct(2),ct(1)-ct(3), ct(2);...
        ct(3), 0         ,-ct(1);...
        ct(1), ct(2), ct(3)      ];
    rhs3 = [e(2)+e(3);e(1);f];
    LHS3 = [ ct(3)-ct(2),ct(1), -ct(1);...
        0          ,-ct(3),     2;...
        ct(1), ct(2), ct(3)      ];
    if det(LHS1)~=0
        LHS = LHS1;
        rhs = rhs1;
    elseif det(LHS2)~=0
        LHS = LHS2;
        rhs = rhs2;
    elseif det(LHS3)~=0
        LHS = LHS3;
        rhs = rhs3;
    else
        error('FIX BUG')
    end
    dt = LHS\rhs;
    d = cp+dt;
else
    d = c;
end
end