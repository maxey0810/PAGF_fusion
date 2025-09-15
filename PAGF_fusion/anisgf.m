function [Y] = anisgf(varargin)
%ANISGF Anisotropic guided filtering of images
%
%   Usage with default parameters:
%       Y = anisgf(X)
%       Y = anisgf(X,G)
%       Y = anisgf(X,G,WG)
%       Y = anisgf(X,[],WG)
%
%   Usage with optional ordered parameters:
%       Y = anisgf(X,scale,gamma,alpha)
%       Y = anisgf(X,G,scale,gamma,alpha)
%       Y = anisgf(X,G,WG,scale,gamma,alpha)
%       Y = anisgf(X,[],WG,scale,gamma,alpha)
%
%   Usage with named parameters:
%       Y = anisgf(...,Name,Value)
%
%   Y = anisgf(X,...) applies the anisotropic guided filter on image X
%   using the same image as the guide and weight guide.
%
%   Y = anisgf(X,G,...) applies the anisotropic guided filter on image X
%   using the guide image G. The multiscale weights are calculated using
%   the same guide image G.
%
%   Y = anisgf(X,G,WG,...) applies the anisotropic guided filter on image X
%   using the guide image G. The multiscale weights are calculated using
%   the weight guide image WG.
%
%   Y = anisgf(X,[],WG,...) applies the anisotropic guided filter on
%   image X using the same input image X as the guide. The multiscale
%   weights are calculated using the weight guide image WG.
%
%
%   Parameters include:
%
%   'scale'     - Positive integer that dictates the window size of the
%                 guided filter. A larger window size will result in
%                 greater smoothing.
%                 Default value: 5 
%
%   'gamma'     - Positive scalar value that controls the smoothness of the
%                 guided filter. A small value preserves more detail while
%                 a large value promotes smoothing
%                 Default value: 0.01
%
%   'alpha'     - Positive scalar value that controls the effect of strong
%                 edges in calculating the anisotropic weights. A large
%                 value strongly avoids regions containing strong edges and
%                 details.
%                 Default value: max(log10(gamma) + 3,0)
%
%   'epsilon'   - Positive scalar value that regularizes the anisotropic
%                 weights. A large value will promote more isotropy while a
%                 small value will emphasize the anisotropic behavior.
%                 Default value: 2^(-8)
%
%   'adaptive'  - Boolean flag that controls the behavior of the guided
%                 filter. Setting this parameter to 'true' enables the
%                 adaptation of the guided filter regularizer (gamma) to
%                 better preserve details in the image.
%                 Default value: true
%
%   'mask'      - Matrix of logical values that excludes pixels from the
%                 guided filter calculations. The matrix should either
%                 contain a single channel or match the number of channels
%                 of the guide image.
%                 Default value: 1 (no mask)

    % Parse the inputs
    global param;
    global flags;
    [X,G,WG,M] = parseInputs(nargin,varargin);

    % Determine the dimensions of the image
    dims = size(X(:,:,1));
    ch = size(X,3);

    % Determine the mode depending on the supplied inputs and initialize
    % the cumulative sums
    if ~flags.useGuide && ~flags.useWeightGuide
        mode = 1;
    elseif flags.useGuide && ~flags.useWeightGuide
        mode = 2;
    elseif ~flags.useGuide && flags.useWeightGuide
        mode = 3;
    else
        mode = 4;
    end
    
    % Estimate the noise variance of the image
    if flags.useWeightGuide
        sigma2 = max((sqrt(pi/2) * sum(sum(abs(imfilter(WG,[1 -2 1; -2 4 -2; 1 -2 1],'replicate')),1),2) / 6 / prod(dims - 2)) .^ 2,[],3);
    elseif flags.useGuide
        sigma2 = max((sqrt(pi/2) * sum(sum(abs(imfilter(G,[1 -2 1; -2 4 -2; 1 -2 1],'replicate')),1),2) / 6 / prod(dims - 2)) .^ 2,[],3);
    else
        sigma2 = max((sqrt(pi/2) * sum(sum(abs(imfilter(X,[1 -2 1; -2 4 -2; 1 -2 1],'replicate')),1),2) / 6 / prod(dims - 2)) .^ 2,[],3);
    end
    
    % Determine the window size of the filter
    windowSize = param.scale;
    windowElem = windowSize .^ 2;

    if flags.useMask
        N = boxFilterFull(double(M),windowSize) + eps;
    else
        N = windowElem;
    end
    
    if flags.useGuide
        G1 = boxFilterFull(G,windowSize);
        G2 = boxFilterFull(G .^ 2,windowSize);
        X1 = boxFilterFull(M .* X,windowSize);
        XG = boxFilterFull(M .* X .* G,windowSize);
        MG1 = boxFilterFull(M .* G,windowSize);
    else
        X1 = boxFilterFull(M .* X,windowSize);
        X2 = boxFilterFull(M .* X .^ 2,windowSize);
    end
    
    if flags.useWeightGuide
        WG1 = boxFilterFull(WG,windowSize);
        WG2 = boxFilterFull(WG .^ 2,windowSize);
    end
    
    % Handle each mode
    if mode == 1 % anisgf(X,...)
        W = (X2 - (X1 .* X1) ./ N) ./ N;
        gamma = adaptiveGamma(W);
        
        A = W ./ (W + gamma);
        B = (1 - A) .* (X1 ./ N);
    elseif mode == 2 % anisgf(X,G,...)
        W = (G2 - (G1 .* G1) / windowElem) / windowElem;
        gamma = adaptiveGamma(W);
        
        A = ((XG - (X1 .* MG1) ./ N) ./ N) ./ (W + gamma);
        B = X1 ./ N - A .* G1 / windowElem;
    elseif mode == 3 % anisgf(X,[],WG,...)
        W = (WG2 - (WG1 .* WG1) / windowElem) / windowElem;
        A = (X2 - (X1 .* X1) ./ N) ./ N;
        gamma = adaptiveGamma(A);
        
        A = A ./ (A + gamma);
        B = (1 - A) .* (X1 ./ N);
    else % anisgf(X,G,WG,...)
        W = (WG2 - (WG1 .* WG1) / windowElem) / windowElem;
        A = (G2 - (G1 .* G1) / windowElem) ./ windowElem;
        gamma = adaptiveGamma(A);
        
        A = ((XG - (X1 .* G1) ./ N) ./ N) ./ (A + gamma);
        B = X1 ./ N - A .* G1 / windowElem;
    end

    % Calculate the weights
    W = max(max(W - sigma2,[],3),0);
    W = param.epsilon ./ ((windowElem * W) .^ param.alpha + param.epsilon);

    % Accumulate the weights to the parameters
    if ~flags.useGaussian
        A = boxFilterValid(W .* A,windowSize);
        B = boxFilterValid(W .* B,windowSize);
        W = boxFilterValid(W,windowSize);
    else
        A = gaussianFilterValid(W .* A,windowSize);
        B = gaussianFilterValid(W .* B,windowSize);
        W = gaussianFilterValid(W,windowSize);
    end

    % Generate the output image
    if flags.useGuide
        Y = (A .* G + B) ./ W;
    else
        Y = (A .* X + B) ./ W;
    end
end

function gamma = adaptiveGamma(V)
    % Calculate the adaptive gamma if necessary
    global flags;
    global param;
    if ~flags.useAdaptive
        gamma = param.gamma;
    else
        % Calculate the median variance from the current scale
        vm = 0.0002 * param.scale;
        gamma = vm * param.gamma ./ V;
    end
end

function Y = boxFilterFull(X,windowSize)
    if windowSize < 15
        X = padarray(X,[windowSize-1 windowSize-1],'replicate');
        Y = convn(X,ones(windowSize),'valid');
    else
        X = padarray(padarray(X,[windowSize windowSize],'pre','replicate'),[windowSize-1 windowSize-1],'post','replicate');
        X = cumsum(X,1);
        X = X(windowSize+1:end,:,:) - X(1:end-windowSize,:,:);
        
        X = cumsum(X,2);
        Y = X(:,windowSize+1:end,:) - X(:,1:end-windowSize,:);
    end
end

function Y = boxFilterValid(X,windowSize)
    if windowSize < 15
        Y = convn(X,ones(windowSize),'valid');
    else
        X = padarray(X,[1 1],'pre');
        
        X = cumsum(X,1);
        X = X(windowSize+1:end,:,:) - X(1:end-windowSize,:,:);
        
        X = cumsum(X,2);
        Y = X(:,windowSize+1:end,:) - X(:,1:end-windowSize,:);
    end
end

function Y = gaussianFilterValid(X,windowSize)
    Y = convn(X,fspecial('gaussian',[windowSize windowSize],windowSize/4),'valid');
end

function [X,G,WG,M] = parseInputs(narg,argin)
    % Initialize the input parser object
    p = inputParser;
    addRequired(p,'X');
    addParameter(p,'scale',5);
    addParameter(p,'gamma',0.01);
    addParameter(p,'alpha',1);
    addParameter(p,'epsilon',2^(-8));
    addParameter(p,'gaussian',false);
    addParameter(p,'adaptive',true);
    addParameter(p,'mask',1);

    % Initialize boolean flags
    global flags;
    flags.useGuide = false;
    flags.useWeightGuide = false;
    flags.useMask = false;
    flags.useAdaptive = false;
    flags.useGaussian = false;
    
    % Determine the arrangement of inputs
    G = [];
    WG = [];
    if narg == 0
        error('No input was found.');
    % anisgf(X)
    elseif narg == 1 
        parse(p,argin{:});
        X = p.Results.X;
        nparam = 0;
    elseif narg == 2
        % anisgf(X,G)
        if size(argin{2},1) > 1
            addRequired(p,'G');
            parse(p,argin{:});
            
            X = p.Results.X;
            G = p.Results.G;
            
            flags.useGuide = true;
            
            nparam = 0;
        % anisgf(X,scale)
        else 
            addOptional(p,'scaleParam',[]);
            addOptional(p,'gammaParam',[]);
            addOptional(p,'alphaParam',[]);
            parse(p,argin{:});
            
            X = p.Results.X;
            
            nparam = 1;
        end
    else
        % anisgf(X,[],WG,...)
        if isempty(argin{2}) 
            addRequired(p,'G');
            addRequired(p,'WG');
            addOptional(p,'scaleParam',[]);
            addOptional(p,'gammaParam',[]);
            addOptional(p,'alphaParam',[]);
            parse(p,argin{:});
            
            X = p.Results.X;
            WG = p.Results.WG;
            
            flags.useWeightGuide = true;
            
            nparam = narg - 3;
        elseif numel(argin{2}) > 1
            % anisgf(X,G,WG,...)
            if numel(argin{3}) > 1
                addRequired(p,'G');
                addRequired(p,'WG');
                addOptional(p,'scaleParam',[]);
                addOptional(p,'gammaParam',[]);
                addOptional(p,'alphaParam',[]);
                parse(p,argin{:});
                
                X = p.Results.X;
                G = p.Results.G;
                WG = p.Results.WG;
                
                flags.useGuide = true;
                flags.useWeightGuide = true;
                
                nparam = narg - 3;
            % anisgf(X,G,...)
            else
                addRequired(p,'G');
                addOptional(p,'scaleParam',[]);
                addOptional(p,'gammaParam',[]);
                addOptional(p,'alphaParam',[]);
                parse(p,argin{:});
                
                X = p.Results.X;
                G = p.Results.G;
                
                flags.useGuide = true;
                
                nparam = narg - 2;
            end
        % anisgf(X,...)    
        else
            addOptional(p,'scaleParam',[]);
            addOptional(p,'gammaParam',[]);
            addOptional(p,'alphaParam',[]);
            parse(p,argin{:});

            X = p.Results.X;

            nparam = narg - 1;
        end
    end
    
    % Assign default parameters
    global param;
    param.scale = p.Results.scale;
    param.gamma = p.Results.gamma;
    param.alpha = p.Results.alpha;
    param.epsilon = p.Results.epsilon;

    % Check for a mask input
    M = 1;
    if ~any(strcmpi(p.UsingDefaults,'mask'))
        M = logical(p.Results.mask);
        flags.useMask = true;
    end
    
    % Check if using gamma adaptation
    if p.Results.adaptive
        flags.useAdaptive = true;
    end
    
    % Check if using gamma adaptation
    if p.Results.gaussian
        flags.useGaussian = true;
    end
    
    % Parse the scale parameter
    parseFlag = true;
    if any(strcmpi(p.UsingDefaults,'scale')) && nparam >= 1 && parseFlag
        if ~isempty(p.Results.gammaParam)
            param.scale = p.Results.scaleParam;
        else
            parseFlag = false;
        end
    end
    validateattributes(param.scale,{'numeric'},{'scalar','positive','integer'},'anisgf','s');
    
    % Parse the gamma parameter
    if any(strcmpi(p.UsingDefaults,'gamma')) && nparam >= 2 && parseFlag
        if ~isempty(p.Results.gammaParam)
            param.gamma = p.Results.gammaParam;
        else
            parseFlag = false;
        end
    end
    validateattributes(param.gamma,{'numeric'},{'scalar','positive'},'anisgf','gamma');
    
    % Parse the alpha parameter
    if any(strcmpi(p.UsingDefaults,'alpha'))
        param.alpha = max(log10(param.gamma) + 3,0);
        if nparam >= 3 && parseFlag
            if ~isempty(p.Results.alphaParam)
                param.alpha = p.Results.alphaParam;
            end
        end
    end
    validateattributes(param.alpha,{'numeric'},{'scalar','nonnegative'},'anisgf','alpha');
    
    % Validate the epsilon parameter
    validateattributes(param.epsilon,{'numeric'},{'scalar','positive'},'anisgf','epsilon');

    % Validate the input images
    validateattributes(X,{'numeric'},{'nonempty','3d'},'anisgf','X');
    dims = size(X(:,:,1));

    if flags.useGuide
        validateattributes(G,{'numeric'},{'nonempty','3d','size',[dims NaN]},'anisgf','G');
        
        if size(G,3) ~= size(X,3)
            if size(G,3) == 3
                G = rgb2gray(G);
            else
                G = mean(G,3);
            end
        end
    end
    
    if flags.useWeightGuide
        validateattributes(WG,{'numeric'},{'nonempty','3d','size',[dims NaN]},'anisgf','G');
    end
    
    if flags.useMask
        validateattributes(M,{'logical'},{'nonempty','3d','size',[dims NaN]},'anisgf','M');
        assert(size(M,3) == size(X,3) || size(M,3) == 1,'Invalid number of mask channels');
    end
end
