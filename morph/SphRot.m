function M = SphRot(a,b)
% Find matrix for unit sphere rotation given sample vector transform
% in cartesian coordinates
% Input: a: old position vector(3d)
%        b: new position vector(3d)
% Output: Linear Transform Matrix

% By Chiyu 'Max' Jiang, May 7, 2016

if norm(a-b)>=1e-6
    c = eye(3);
    M = zeros(3,3);
    M(:,1) = Rot(a,b,c(:,1));
    M(:,2) = Rot(a,b,c(:,2));
    M(:,3) = Rot(a,b,c(:,3));
else
    M = ones(3,3);
end
end



function d = Rot(a,b,c)

if norm(a)~=1
    a = a/norm(a);
end
if norm(b)~=1
    b = b/norm(b);
end
n = cross(a,b)/norm(cross(a,b));
theta = acos(dot(a,b)/(norm(a)*norm(b)));
cp = (dot(c,n))*n;
ct = c-cp;
if norm(ct)~=0
    f = norm(ct)^2*(dot(a,b)/(norm(a)*norm(b)));
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