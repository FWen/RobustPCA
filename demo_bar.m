clear all; clc; close all;

X = double(imread('barbara.png'));
[m,n] = size(X);

[S,V,D] = svd(X);
v = diag(V);

figure(1);
plot(1:min(m,n),v/max(v),'-'); xlim([1,min(m,n)]);
xlabel('Singular value index'); ylabel('Normalized amplitude');


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

figure(2);subplot(1,3,1); imshow(uint8(Y));
title(sprintf('Corrupted\n RelErr=%.3f, PSNR=%.2f dB',norm(Y - X,'fro')/norm(X,'fro'),psnr(X,Y)));


param.TOL =1e-6;
param.MAX_ITER = 1e3;

% --- soft thresholding --------------------
q1 = 1; q2 = 1;
[X_L1, S_L1, out] = lq_lq_l2_bcd(Y, 0.05, q1, q2, param);

figure(2);subplot(1,3,2);
imshow(uint8(X_L1));
title(sprintf('Soft\n RelErr=%.3f, PSNR=%.2f dB',norm(X_L1(:)-X(:))/norm(X(:)),psnr(X,X_L1)));

% ----- Lq thresholding ----------
param.L0 = X_L1;
param.S0 = S_L1;
q1 = 0.8; q2 = 0.6;
[X_Lq, ~, ~] = lq_lq_l2_bcd(Y, 0.07, q1, q2, param);

figure(2);subplot(1,3,3);
imshow(uint8(X_Lq));
title(sprintf('Lq (q1=%.1f, q2=%.1f)\n RelErr=%.3f, PSNR=%.2f dB',q1,q2,norm(X_Lq(:)-X(:))/norm(X(:)),psnr(X, X_Lq)));

