function R = PAGF(Im, nlev, scale,SIZE,gamma)

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
r = cell(nlev-1,1);
u = cell(nlev-1,1);
g=cell(nlev-1,1);

R=G{size(G,1)};
for l=size(G,1)-1:-1:1
    % upsample
    R_up = imresize(R,[size(G{l},1),size(G{l},2)], 'bilinear');
    u{l}=R_up;

    %R_hat = zeros([size(G{l}, 1), size(G{l}, 2), 3]); 
    SIZE1=max(5,ceil(SIZE*0.5));
    L{l}=anisgf(L{l},R_up,[],'scale',SIZE1,'gamma',gamma);

    R_hat = R_up+L{l};
    g{l}=R_hat;


    % laplacian
    R=anisgf(G{l},R_hat,[],'scale',SIZE,'gamma',gamma);
    r{l}=R;

end



