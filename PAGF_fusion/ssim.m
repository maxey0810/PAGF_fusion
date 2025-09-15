function ssim_map = ssim(img1,fused)  


    w = fspecial('gaussian', 11, 1.5);  
    K(1) = 0.01;                      
    K(2) = 0.03;                      
    L = 255;       

    s=size(size(img1));
    if s(2)==3 
        img1=rgb2gray(img1);
    end 

    s1=size(size(fused));
    if s1(2)==3 
        fused=rgb2gray(fused);
    end 

    img1 = 255*double(img1);  
    fused = 255*double(fused);  

    C1 = (K(1)*L)^2;  
    C2 = (K(2)*L)^2;  
    w = w/sum(sum(w));  

    ua   = filter2(w, img1, 'same');  
    ub   = filter2(w, fused, 'same');
    ua_sq = ua.*ua;  
    ub_sq = ub.*ub;  
    ua_ub = ua.*ub;  
    siga_sq = filter2(w, img1.*img1, 'same') - ua_sq;  
    sigb_sq = filter2(w, fused.*fused, 'same') - ub_sq;  
    sigab = filter2(w, img1.*fused, 'same') - ua_ub;  

    ssim_map = ((2*ua_ub + C1).*(2*sigab + C2))./((ua_sq + ub_sq + C1).*(siga_sq + sigb_sq + C2));     
    
end