function vo = SphNorm(vi)
% Normalize coordinates to unit sphere
% vi and vo are coordinates in cartesian, n by 3

vs = zeros(size(vi)); vo = vs;
[vs(:,1),vs(:,2),vs(:,3)] = cart2sph(vi(:,1),vi(:,2),vi(:,3));
[vo(:,1),vo(:,2),vo(:,3)] = sph2cart(vs(:,1),vs(:,2),ones(size(vs(:,3))));
end