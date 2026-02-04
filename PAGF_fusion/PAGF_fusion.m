function F=PAGF_fusion(IR,VI,b_fixed)
%% base layer

B_IR = PAGF(IR,5,0.8,5,0.01);

B_VI=PAGF2(VI, 5, 0.5, 8, 1);
B_VI=PAGF2(B_VI, 5, 0.5, 8, 1);

D_IR=IR-B_IR;
D_VI=VI-B_VI;

gaussianfilter = fspecial('gaussian', 5, 1);
%% fuion of base layer
%saliency detection
S_IR=abs(0.5*(mean(IR(:))+median(IR(:)))-IR);
S_VI=abs(0.5*(mean(VI(:))+median(VI(:)))-VI);
%Normalization
S_IR=S_IR./max(S_IR(:));
S_VI=S_VI./max(S_VI(:));

W_BIR=0.5+0.5*(S_IR-S_VI);

if nargin < 3 || isempty(b_fixed)
SSIM=ssim(S_VI,S_IR);
THETA=zeros(size(IR));
THETA(S_IR>=S_VI)=0.5-SSIM(S_IR>=S_VI);
THETA(S_VI>S_IR)=SSIM(S_VI>S_IR)-1.5;
b=exp(THETA);
%b(S_IR>=S_VI)=2*log(2-SSIM(S_IR>=S_VI));
else
    b = b_fixed;
end

W_BIR=atan(b.*W_BIR)./atan(b);


F_B=W_BIR.*B_IR + (1-W_BIR).*B_VI;


%% fuion of detail layer
    W_DVI = zeros(size(IR));

    pir=phasecong(D_IR);
    pvi=phasecong(D_VI);

    W_DVI(pvi >= pir) = 1;
   
    W_DVI = imfilter(W_DVI,gaussianfilter,'symmetric');

    W_DIR=1-W_DVI;

    F_D = W_DIR .* D_IR + W_DVI .* D_VI;


F=F_B+F_D;
end
