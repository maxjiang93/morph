function [F,varargout] = morph_multi(varargin)
% Morph two shapes given triangle mesh
% Acceptable formats: .off, .ply, .wrl, .obj, .m, .gim.
% Format: [F,VS1,VS2,...] = morph_multi(n_feat,filename1,filename2,...)
% Inputs:
%                        n_feat : number of feature points to match
%                                 order of importance: 1>2>3,4,5,6...
%     filename1, filename2, ... : a string for the file name of 1st mesh file
% Outputs:
%
% By Chiyu 'Max' Jiang, May 7, 2016

addpath(genpath('~/Desktop/Codes/gptoolbox/'))
assert(nargin>=3,'ERROR: need at least two mesh files to morph')

n_feat = varargin{1};
nfile = nargin-1;
filename = cell(nfile,1);
for im = 1:nfile
    filename{im} = varargin{im+1};
end


%% Read Mesh
Vi = cell(nfile,1);
Fi = cell(nfile,1);
for im = 1:nfile
    [Vi{im},Fi{im}] = read_mesh(filename{im});
    if size(Vi{im},2)~=2 && size(Vi{im},2)~=3
        Vi{im} = Vi{im}'; Fi{im} = Fi{im}';
    end
end
v_id = zeros(n_feat,nfile); % Allocate n*nfile matrix to store corrsponding vertex id


%% Interactively choose feeature points in mesh (except last mesh)
fprintf('[Step 1/5] Choose Feature Points in Mesh\n')
fig = cell(nfile,1);
for im = 1:nfile-1
    fig{im} = figure(im);
    tsurf(Fi{im},Vi{im});
    hold on
    dcm_obj = datacursormode(fig{im});
    set(dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','on','Enable','on')
    for i = 1:n_feat
        title(['Click to choose feature point #' num2str(i) ', then press Return.'])
        % Wait while the user does this.
        pause
        c_info = getCursorInfo(dcm_obj);
        % Get position coordinates
        pos = c_info.Position;
        % Place marker
        text(pos(1),pos(2),pos(3),['Pt.' num2str(i)],'BackgroundColor',[.7 .7 .7])
        scatter3(pos(1),pos(2),pos(3),...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor','r')
        
        % Find vertex id as in array V
        a = intersect(find((pos(1)==Vi{im}(:,1))),find((pos(2)==Vi{im}(:,2))));
        b = intersect(a,find((pos(3)==Vi{im}(:,3))));
        if isempty(b)~=1
            vindex = b;
        else
            error('ERROR: failed to find vertex number')
        end
        v_id(i,im) = vindex;
    end
    hold off
    title(['Finished choosing feature points for mesh ' num2str(im) '. Click Enter to open mesh ' num2str(im+1)])
    pause
end

%% Interactively choose feeature points in final mesh
im = nfile;
fig{im} = figure(im);
tsurf(Fi{im},Vi{im});
hold on
dcm_obj = datacursormode(fig{im});
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','on','Enable','on')
for i = 1:n_feat
    title(['Click to choose feature point #' num2str(i) ', then press Return.'])
    % Wait while the user does this.
    pause
    c_info = getCursorInfo(dcm_obj);
    % Get position coordinates
    pos = c_info.Position;
    % Place marker
    text(pos(1),pos(2),pos(3),['Pt.' num2str(i)],'BackgroundColor',[.7 .7 .7])
    scatter3(pos(1),pos(2),pos(3),...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r')
    
    % Find vertex id as in array V
    a = intersect(find((pos(1)==Vi{im}(:,1))),find((pos(2)==Vi{im}(:,2))));
    b = intersect(a,find((pos(3)==Vi{im}(:,3))));
    if isempty(b)~=1
        vindex = b;
    else
        error('ERROR: failed to find vertex number')
    end
    v_id(i,im) = vindex;
end
hold off
title(['Finished choosing feature points for mesh' num2str(im) '. Click Enter to close all mesh.'])
pause
for im = 1:nfile
    close(fig{im})
end

%% Conformal Mean Curvature Flow to map to sphere
fprintf('[Step 2/5] Spherical Parametirization\n')
Ui = cell(nfile,1);
for im = 1:nfile
    Ui{im} = conformalized_mean_curvature_flow(Vi{im},Fi{im},'MaxIter',100,'delta',1e-2);
end

%% Check Parameterized Mesh
% Check mesh
for im = 1:nfile
    fig{im} = figure(im);
    tsurf(Fi{im},Ui{im})
    hold on
    scatter3(Ui{im}(v_id(:,im),1),Ui{im}(v_id(:,im),2),Ui{im}(v_id(:,im),3),...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r')
    hold off
    title('Check if curvature flow successful. Enter to check next mesh.')
    pause
end
for im = 1:nfile
    close(fig{im});
end

%% Normalize to Unit Sphere
for im = 1:nfile
    [Phi,Theta,R] = cart2sph(Ui{im}(:,1),Ui{im}(:,2),Ui{im}(:,3));
    R = ones(size(R)); % Normalize
    [Ui{im}(:,1),Ui{im}(:,2),Ui{im}(:,3)] = sph2cart(Phi,Theta,R);
end

%% Rotate to match feature point 1
b = Ui{1}(v_id(1,1),:);
for im = 2:nfile
    a = Ui{im}(v_id(1,im),:);
    M = SphRot(a',b');
    Ui{im} = M*Ui{im}'; Ui{im} = Ui{im}'; % Rotate Sphere
    Vi{im} = M*Vi{im}'; Vi{im} = Vi{im}'; % Rotate Original Mesh
end

%% Rotate to close-up feature point 2
a = Ui{1}(v_id(1,1),:); % Pivot Points
b = Ui{1}(v_id(2,1),:);
for im = 2:nfile
    c = Ui{im}(v_id(2,im),:);
    M = SphRot2(a',b',c');
    Ui{im} = M*Ui{im}'; Ui{im} = Ui{im}'; % Rotate Sphere
    Vi{im} = M*Vi{im}'; Vi{im} = Vi{im}'; % Rotate Original Mesh
end

%% Warp mesh using arap algorithm
fprintf('[Step 3/5] Mesh Warping\n')

% Choice 1, warp every other mesh to match mesh 1
% Choise 2, warp every single mesh to their middle ground
warp_choice = 2;
if warp_choice == 1
    for im = 2:nfile
        Ui{im} = iter_rbf_warp(Ui{im},v_id(:,im),Ui{1}(v_id(:,1),:),50);
    end
elseif warp_choice == 2
    a = zeros(1,3);
    for im = 1:nfile
        a = a+Ui{im}(v_id(:,im),:);
    end
    a = a/3;
    middle = SphNorm(a);
    for im = 1:nfile
        Ui{im} = iter_rbf_warp(Ui{im},v_id(:,im),middle,50);
    end
else
    error('warp_choice must be 1 or 2')
end
% Save 2

%% Check Mesh After Warping
% Check mesh
for im = 1:nfile
    fig{im} = figure(im);
    tsurf(Fi{im},Ui{im})
    hold on
    scatter3(Ui{im}(v_id(:,im),1),Ui{im}(v_id(:,im),2),Ui{im}(v_id(:,im),3),...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','r')
    hold off
    title('Check if curvature flow successful. Enter to check next mesh.')
    pause
end
for im = 1:nfile
    close(fig{im});
end

%% Construct Supermesh
fprintf('[Step 4/5] Construct Supermesh\n')
SuperMeshSize = 0;
ni = zeros(nfile,1);
for im = 1:nfile
    SuperMeshSize = SuperMeshSize + size(Vi{im},1);
    ni(im) = size(Vi{im},1);
end
U = zeros(SuperMeshSize,3);
U(1:ni(1),:) = Ui{1};
for im = 2:nfile
    nprev = sum(ni(1:im-1));
    nnow = sum(ni(1:im));
    U(nprev+1:nnow,:) = Ui{im};
end
F = convhulln(U);
VS = cell(nfile,1);
VS{1} = zeros(size(U));
VS{1}(1:ni(1),:) = Vi{1};
for im = 2:nfile
    VS{im} = zeros(size(U));
    nprev = sum(ni(1:im-1));
    nnow = sum(ni(1:im));
    VS{im}(nprev+1:nnow,:) = Vi{im};
end

o = [0,0,0];

% Interpolate Mesh (jm) based on Mesh (im)
count = 1;
unfound = cell(nfile,nfile);
adj = cell(nfile,nfile);
for im = 1:nfile
    for jm = 1:nfile
        if jm ~= im
            
            if jm>1
                nprev = sum(ni(1:(jm-1)));
            else
                nprev = 0;
            end
            
            for i = 1:ni(jm)
                d = Ui{jm}(i,:);
                [flag,~,lambda] = ray_mesh_intersect(o,d,Ui{im},Fi{im});
                fflag = find(flag);
                if isempty(intersect(i,v_id(:,jm)))~=1 % If this is a feature point, dont need to interpolate
                    vert = find(v_id(:,jm)==i);
                    VS{im}(i+nprev,:) = Vi{im}(v_id(vert,im),:);
                elseif length(fflag)==1
                    vert = Fi{im}(fflag,:);
                    VS{im}(i+nprev,:) = lambda(fflag,1).*Vi{im}(vert(1),:)+lambda(fflag,2).*Vi{im}(vert(2),:)...
                        +lambda(fflag,3).*Vi{im}(vert(3),:);
%                 elseif length(fflag)>1
%                     fflag = fflag(1);
%                     vert = Fi{im}(fflag,:);
%                     VS{im}(i+nprev,:) = lambda(fflag,1).*Vi{im}(vert(1),:)+lambda(fflag,2).*Vi{im}(vert(2),:)...
%                         +lambda(fflag,3).*Vi{im}(vert(3),:);
                else
                    std = lambda(:,1).^2+lambda(:,2).^2+lambda(:,3).^2;
                    crit = (std<1)+(lambda(:,1)>0)+(lambda(:,2)>0)+(lambda(:,3)>0)==4;
                    fflag = find(crit);
                    if length(fflag)>1
                        fflag = fflag(1);
                    end
                    if not(isempty(fflag))
                        vert = Fi{im}(fflag,:);
                        VS{im}(i+nprev,:) = lambda(fflag,1).*Vi{im}(vert(1),:)+lambda(fflag,2).*Vi{im}(vert(2),:)...
                            +lambda(fflag,3).*Vi{im}(vert(3),:);
                    else
                        unfound{im,jm} = [unfound{im,jm},i];
                        warning(['Unfound element when interpolating mesh ' num2str(jm) ' based on mesh ' num2str(im)])
                    end
                end
            end
            
            if isempty(unfound{im,jm})~=1
                adj{im,jm} = adjacency_list(Fi{jm},unfound{im,jm});
                for i = 1:length(unfound{im,jm})
                    if isempty(intersect(v_id(:,jm),adj{im,jm}{i}))~=1
                        vert = intersect(v_id(:,jm),adj{im,jm}{i});
                        vert = vert(1);
                        vert = find(v_id(:,jm)==vert);
                        VS{im}(unfound{im,jm}(i)+nprev,:) = Vi{im}(v_id(vert,im),:);
                    else
                        [~,nearest] = min(abs(adj{im,jm}{i}(:)-unfound{im,jm}(i)));
                        VS{im}(unfound{im,jm}(i)+nprev,:) = VS{im}(adj{im,jm}{i}(nearest),:);
                        warning(['Unresolved unfound vertex when interpolating mesh ' num2str(jm) ' based on mesh ' num2str(im)])
                    end
                end
            end
            
            fprintf(['[Step 4/5] Construct Supermesh ' num2str(count) ' of ' num2str(nchoosek(nfile,2)*2) ' Completed\n'])
            count = count+1;
        end
    end
end

%% Export Supermesh Matrix
varargout = cell(nfile,1);
for im = 1:nfile
    varargout{im} = VS{im};
end

