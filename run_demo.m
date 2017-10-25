% This is a simple script for running a demo testcase
% By Max Chiyu Jiang, July 13, 2016

[~,current_directory] = system('pwd'); % Only work on unix machines
current_directory(end) = [];

%% Load External Toolboxes, morphing code and mesh files
addpath(genpath([current_directory '/']))

%% Run Morphing Algorithm with cow and horse
n_feat = 9; % Number of feature points to match between models
filename1 = 'cow40k.ply';
filename2 = 'horse50k.ply';
[F,VS1,VS2] = morph_multi(n_feat,filename1,filename2);

%% Animate Morphed Shape
fprintf('Mesh Processing Complete. Starting Animation.\n')

% Preallocate movie structure
nFrames = 101; % Number of frames
mov(1:nFrames) = struct('cdata', [],...
    'colormap', []);
% Animate
alpha = 0; % Morphing ratio
d_alpha = 1/(nFrames-1); % Step size for alpha
for it = 1:nFrames
    V_morph = alpha*VS1+(1-alpha)*VS2;
    h = trisurf(F,V_morph(:,1),V_morph(:,2),V_morph(:,3),'EdgeColor','none','FaceColor','flat');
%     h = tsurf(F,V_morph);
    lighting gouraud
    axis equal;
    set(h,'FaceAlpha',1)
    title(['alpha = ' num2str(alpha)])
    view(2)
    mov(it) = getframe(gcf);
    alpha = alpha+d_alpha;
end