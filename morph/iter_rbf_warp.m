function [U,Ustep] = iter_rbf_warp(V,b,bc,n)

% Function to warp mesh on sphere using radial basis function
% Reference: "Mesh deformation based on radial basis function interpolation"
% Boer et.al 2007
% process is subdivided into n smaller steps
% Input:  V : Vertices # by 3, list of vertices coordinates
%         b : #b by 1, list of control point index
%         bc: #b by 3, list of terminal position for control points
%         n : number of iterations
% Output: U : V# by 3, list of final position for all vertices

% By Chiyu 'Max' Jiang, May 8, 2016

%% Normalize the inianial vertex and final control positions
V = SphNorm(V);
bc = SphNorm(bc);

%% Find trajectory of control points
Ubs = V(b',:); % Start positions of control points
Ube = bc;     % End positions of control points
Theta = real(acos(dot(Ubs,Ube,2)));
N = n+1;
alpha = zeros(length(b),N);
if mod(N,2)==0 % even
    for i = 1:N/2
        alpha(:,i) = (sin(Theta/2)-tan(Theta/2-Theta/n*(i-1)).*cos(Theta/2))./(2*sin(Theta/2));
        alpha(find(Theta==0),i) = zeros(size(find(Theta==0)));
        alpha(:,end:N/2+1) = 1-alpha(:,1:N/2);
    end
else
    for i = 1:floor(N/2)
        alpha(:,i) = (sin(Theta/2)-tan(Theta/2-Theta/n*(i-1)).*cos(Theta/2))./(2*sin(Theta/2));
        alpha(find(Theta==0),i) = zeros(size(find(Theta==0)));
        alpha(:,end+1-i) = 1-alpha(:,i);
    end
    alpha(:,ceil(N/2)) = 0.5;
end
traj_line = zeros(length(b),N,3); traj_curve = zeros(length(b),n,3);
for i = 1:N
    traj_line(:,i,:) = Ubs+alpha(:,i).*(Ube-Ubs);
end
traj_line = traj_line(:,2:end,:);
for i = 1:n
    traj_curve(:,i,:) = SphNorm(traj_line(:,i,:));
end

%% Step through trajectory
Ustep = zeros(size(V,1),size(V,2),n+1);
Ustep(:,:,1) = V;
for i = 2:n+1
    Ustep(:,:,i) = rbf_warp(Ustep(:,:,i-1),b,reshape(traj_curve(:,i-1,:),[length(b),3]));
end
U = Ustep(:,:,end);


end

%% Normalize to unit sphere subroutine
function vo = SphNorm(vi)
% Normalize coordinates to unit sphere
% vi and vo are coordinates in cartesian, n by 3

vs = zeros(size(vi)); vo = vs;
[vs(:,1),vs(:,2),vs(:,3)] = cart2sph(vi(:,1),vi(:,2),vi(:,3));
[vo(:,1),vo(:,2),vo(:,3)] = sph2cart(vs(:,1),vs(:,2),ones(size(vs(:,3))));
end

%% Warp by radial basis function subroutine
function U = rbf_warp(V,b,bc)

Vb = V(b',:);
nb = length(b);
Mbb = zeros(nb,nb);
for i = 1:nb
    for j = 1:nb
        Mbb(i,j) = rbf(norm(Vb(i,:)-Vb(j,:)));
    end
end
Pb = [ones(nb,1),Vb];
db = bc-Vb;
db = [db;zeros(4,3)];
M = [[Mbb;Pb'],[Pb;zeros(4,4)]];
ab = zeros(nb+4,3);
for i = 1:3
    ab(:,i) = M\db(:,i);
end
alpha = ab(1:nb,:);
beta = ab(nb+1:end,:);
d = zeros(size(V));
for i = 1:size(V,1)
    delta = V(i,:)-Vb;
    delta = sqrt(delta(:,1).^2+delta(:,2).^2+delta(:,3).^2);
    d(i,:) = sum(rbf(delta).*alpha,1)+beta(1,:)+sum(V(i,:)'.*beta(2:4,:),1);
end
U = SphNorm(V+d);
end

%% Radial Basis Function
function phi = rbf(xi)
r = 1;
xi = xi/r;
% Radial Basis Function: CTPS C_a^2
z = double(xi==0);
xi = xi+z;
phi = double(xi==0)+(xi>0).*(xi<=1).*(1-30*xi.^2-10*xi.^3+45*xi.^4-6*xi.^5-60*xi.^3.*log(xi));
phi = phi.*(z==0)+z;
end