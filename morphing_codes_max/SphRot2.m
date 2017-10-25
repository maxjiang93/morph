function M = SphRot2(a,b,c)
% Find matrix for unit sphere rotation given sample vector transform
% in cartesian coordinates
% Input: a: position of pivot point(3d)
%        b: position of 2nd feature point to match(3d)
%        c: 2nd feature point to rotate
% Output: Linear Transform Matrix

% By Chiyu 'Max' Jiang, May 7, 2016

n1 = cross(a,b-a); n1 = n1/norm(n1);
n2 = cross(a,c-a); n2 = n2/norm(n2);

theta = acos(dot(n1,n2));
M1 = SphRotation(a,theta);
M2 = SphRotation(a,-theta);
c1 = M1*c;
c2 = M2*c;
nc1 = cross(a,c1-a); nc1 = nc1/norm(nc1);
nc2 = cross(a,c2-a); nc2 = nc2/norm(nc2);

theta1 = acos(dot(n1,nc1));
theta2 = acos(dot(n1,nc2));

if abs(real(theta1))<1e-6
    M = M1;
elseif abs(real(theta2))<1e-6
    M = M2;
else
    M = eye(3,3);
    fprintf('Rotation 2 failed\n')
end