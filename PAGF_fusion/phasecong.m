function [ResultPC] = phasecong(im)


nscale          = 4;     
norient         = 4;    
minWaveLength   = 6;    
mult            = 2;     
sigmaOnf        = 0.55; 
dThetaOnSigma   = 1.2;  
epsilon         = .0001; 

[rows, cols] = size(im);
imagefft = fft2(im);
zero = zeros(rows, cols);
EO = cell(nscale, norient);
ifftFilterArray = cell(1, nscale);

if mod(cols, 2)
    xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
else
    xrange = [-cols/2:(cols/2-1)]/cols;	
end

if mod(rows, 2)
    yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
else
    yrange = [-rows/2:(rows/2-1)]/rows;	
end

[x, y] = meshgrid(xrange, yrange);
radius = sqrt(x.^2 + y.^2);
theta = atan2(-y, x);
radius = ifftshift(radius);
theta = ifftshift(theta);
radius(1,1) = 1;

sintheta = sin(theta);
costheta = cos(theta);
clear x y theta;

thetaSigma = pi/norient/dThetaOnSigma;
lp = lowpassfilter([rows, cols], .45, 15);

logGabor = cell(1, nscale);
for s = 1:nscale
    wavelength = minWaveLength * mult^(s-1);
    fo = 1.0/wavelength;
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s} .* lp;
    logGabor{s}(1,1) = 0;
end

spread = cell(1, norient);
for o = 1:norient
    angl = (o-1) * pi / norient;
    ds = sintheta * cos(angl) - costheta * sin(angl);
    dc = costheta * cos(angl) + sintheta * sin(angl);
    dtheta = abs(atan2(ds, dc));
    spread{o} = exp((-dtheta.^2) / (2 * thetaSigma^2));
end

for o = 1:norient
    for s = 1:nscale
        filter = logGabor{s} .* spread{o};
        EO{s,o} = ifft2(imagefft .* filter);
    end
end
%% Adaptive K
    totalEnergyAllOrients = zeros(rows, cols);

    for x=1:norient
    E2_orient = abs(EO{1,x}).^2;

    totalEnergyAllOrients = totalEnergyAllOrients + E2_orient;
    end

    medianNoiseEnergy = 0.25*median(reshape(totalEnergyAllOrients, 1, rows*cols));


EnergyAll = zero;
AnAll = zero;

for o = 1:norient
    sumE_ThisOrient = zero;
    sumO_ThisOrient = zero;
    sumAn_ThisOrient = zero;
    Energy = zero;
    
    for s = 1:nscale
        ifftFilt = real(ifft2(filter)) * sqrt(rows*cols);
        ifftFilterArray{s} = ifftFilt;
        
        An = abs(EO{s,o});
        sumAn_ThisOrient = sumAn_ThisOrient + An;
        sumE_ThisOrient = sumE_ThisOrient + real(EO{s,o});
        sumO_ThisOrient = sumO_ThisOrient + imag(EO{s,o});
        
        if s == 1
            EM_n = sum(sum(filter.^2));
            maxAn = An;
        else
            maxAn = max(maxAn, An);
        end
    end
    
    XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;
    MeanE = sumE_ThisOrient ./ XEnergy;
    MeanO = sumO_ThisOrient ./ XEnergy;
    
    for s = 1:nscale
        E = real(EO{s,o});
        O = imag(EO{s,o});
        Energy = Energy + E.*MeanE + O.*MeanO - abs(E.*MeanO - O.*MeanE);
    end
    
    %% Noise power estimation   
    medianE2n = median(reshape(abs(EO{1,o}).^2, 1, rows*cols));
    meanE2n = -medianE2n / log(0.5);
    noisePower = meanE2n / EM_n; 
    
    
    EstSumAn2 = zero;
    for s = 1:nscale
        EstSumAn2 = EstSumAn2 + ifftFilterArray{s}.^2;
    end
    
    EstSumAiAj = zero;
    for si = 1:(nscale-1)
        for sj = (si+1):nscale
            EstSumAiAj = EstSumAiAj + ifftFilterArray{si} .* ifftFilterArray{sj};
        end
    end
    
    sumEstSumAn2 = sum(sum(EstSumAn2));
    sumEstSumAiAj = sum(sum(EstSumAiAj));
    EstNoiseEnergy2 = 2 * noisePower * sumEstSumAn2 + 4 * noisePower * sumEstSumAiAj;
    tau = sqrt(EstNoiseEnergy2/2);

    EstNoiseEnergy = tau * sqrt(pi/2); 
    EstNoiseEnergySigma = sqrt((2-pi/2) * tau^2); 


        a = 4;
        b = 10000;
        c = 2;
        k = a * log(b * medianNoiseEnergy + 1) + c;
   
        k = max(2, min(10, k))
        

    T = EstNoiseEnergy + k * EstNoiseEnergySigma;
    T = T / 1.7; 
    Energy = max(Energy - T, zero);
    
    EnergyAll = EnergyAll + Energy;
    AnAll = AnAll + sumAn_ThisOrient;
end

ResultPC = EnergyAll ./ AnAll;

return;




function f = lowpassfilter(sze, cutoff, n)

    if cutoff < 0 || cutoff > 0.5
	error('cutoff frequency must be between 0 and 0.5');
    end
    
    if rem(n,1) ~= 0 || n < 1
	error('n must be an integer >= 1');
    end

    if length(sze) == 1
	rows = sze; cols = sze;
    else
	rows = sze(1); cols = sze(2);
    end

    % Set up X and Y matrices with ranges normalised to +/- 0.5
    % The following code adjusts things appropriately for odd and even values
    % of rows and columns.
    if mod(cols,2)
	xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
    else
	xrange = [-cols/2:(cols/2-1)]/cols;	
    end

    if mod(rows,2)
	yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
    else
	yrange = [-rows/2:(rows/2-1)]/rows;	
    end
    
    [x,y] = meshgrid(xrange, yrange);
    radius = sqrt(x.^2 + y.^2);        % A matrix with every pixel = radius relative to centre.
    f = ifftshift( 1 ./ (1.0 + (radius ./ cutoff).^(2*n)) );   % The filter
    return;
