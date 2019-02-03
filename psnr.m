function PSNR = psnr(I1,I2);
% Compute the peak SNR given two images

[m,n]=size(I1);

I1=double(I1(:));
I2=double(I2(:));

MSE = 1/(m*n)*sum((I1-I2).^2);

MAXi = max(abs(I1));

PSNR = 	20*log10(MAXi/sqrt(MSE));

return;
