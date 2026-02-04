function R = PAGF2(Im, nlev, scale,SIZE,gamma)

if ~exist('nlev','var')
    nlev = 3;
end

if ~exist('scale','var')
    scale = 0.8;
end

if ~exist('gamma','var')
    gamma = 0.01;
end
if ~exist('SIZE','var')
    SIZE = 5;
end

[G,L]=pyramid(Im, nlev, scale);

R=G{size(G,1)};
for l=size(G,1)-1:-1:1
    % upsample
    R_up = imresize(R,[size(G{l},1),size(G{l},2)], 'bilinear');
    %R_hat = zeros([size(G{l}, 1), size(G{l}, 2), 3]);
    L{l}=anisgf(L{l},R_up,[],'scale',SIZE,'gamma',gamma);
    R = R_up+L{l};

    % laplacian
    
    
end



