clear all; clc; close all;

m = 256;
n = m;
X = imresize(double(imread('Cameraman.png')),[m, n]);


%---strictly low-rank-----
[S,V,D] = svd(X);
v = diag(V);
v(round(0.15*m):end) = 0;
X = S*diag(v)*D';


%---salt-and-pepper noise-----
c_ratio = 0.2; % corruption ratio
J = randperm(m*n); 
J = J(1:round(c_ratio*m*n));    
y = X(:);
noise = randn(1,round(c_ratio*m*n));
noise(find(noise>=0))=255;
noise(find(noise<0))=0;
y(J) = noise;
Y = reshape(y, m, n);

figure(1);subplot(1,3,1); imshow(uint8(Y));
title(sprintf('Corrupted\n RelErr=%.3f, PSNR=%.2f dB',norm(Y - X,'fro')/norm(X,'fro'),psnr(X,Y)));


param.TOL =1e-16;
param.MAX_ITER = 2e3;
lamdas = logspace(-5, 2, 50);

% --- soft thresholding --------------------
q1 = 1; q2 = 1;
[X_L1, S_L1, out] = lq_lq_l2_bcd(Y, 0.071, q1, q2, param);

figure(1);subplot(1,3,2);
imshow(uint8(X_L1));
title(sprintf('Soft\n RelErr=%.3f, PSNR=%.2f dB',norm(X_L1(:)-X(:))/norm(X(:)),psnr(X,X_L1)));

% ----- Lq thresholding ----------
param.L0 = X_L1;
param.S0 = S_L1;
q1 = 0; q2 = 0;
[X_Lq, ~, ~] = lq_lq_l2_bcd(Y, 0.02, q1, q2, param);

figure(1);subplot(1,3,3);
imshow(uint8(X_Lq));
title(sprintf('Lq (q1=%.1f, q2=%.1f)\n RelErr=%.3e, PSNR=%.2f dB',q1,q2,norm(X_Lq(:)-X(:))/norm(X(:)),psnr(X, X_Lq)));
