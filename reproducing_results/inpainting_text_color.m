clear all; %clc;  close all;

X = imread('turtle.png');
Y = imread('masked_turtle.png');%corrupted by text

[m,n,~] = size(X);
Xd = double(X);
Yd = double(Y);

figure(1);
subplot(1,4,1); imshow(uint8(Y));
title(sprintf('Corrupted\n RelErr=%.3f, PSNR=%.2f dB',norm(Yd(:) - Xd(:))/norm(Xd(:)),psnr(Xd(:),Yd(:))));


param.TOL =1e-6;
param.MAX_ITER = 1e3;
lamdas = logspace(-3, 0, 20);

% --- soft thresholding --------------------
for ch=1:3
    parfor k = 1:length(lamdas)
        [L, S, ~] = lq_lq_l2_bcd(Yd(:,:,ch), lamdas(k), 1, 1, param);         
        relerr(k)  = norm(L-Xd(:,:,ch),'fro')/norm(Xd(:,:,ch),'fro');
        X_l1(:,:,k) = L;
        S_l1(:,:,k) = S;
    end
    [RelErr_L1, mi] = min(relerr);

    X_L1(:,:,ch) = X_l1(:,:,mi); 
    S_L1(:,:,ch) = S_l1(:,:,mi); 
end

figure(1);subplot(1,4,2);
imshow(uint8(X_L1));
title(sprintf('Soft\n RelErr=%.3f, PSNR=%.2f dB',norm(X_L1(:) - Xd(:),'fro')/norm(Xd(:),'fro'),psnr(X(:),X_L1(:))));


% --- hard thresholding --------------------
for ch=1:3
    param.L0 = X_L1(:,:,ch);
    param.S0 = S_L1(:,:,ch);
    parfor k = 1:length(lamdas)
        [L, ~, ~] = lq_lq_l2_bcd(Yd(:,:,ch), lamdas(k), 0, 0, param);
        relerr(k)  = norm(L-Xd(:,:,ch),'fro')/norm(Xd(:,:,ch),'fro');
        X_l0(:,:,k) = L;
    end
    [~, mi] = min(relerr);
    X_L0(:,:,ch) = X_l0(:,:,mi); 
end

figure(1);subplot(1,4,3);
imshow(uint8(X_L0));
title(sprintf('Hard\n RelErr=%.3f, PSNR=%.2f dB',norm(X_L0(:) - Xd(:))/norm(Xd(:)),psnr(X(:),X_L0(:))));


% --- Lq thresholding --------------------
for ch=1:3
    param.L0 = X_L1(:,:,ch);
    param.S0 = S_L1(:,:,ch);
    parfor k = 1:length(lamdas)
        [L, ~, ~] = lq_lq_l2_bcd(Yd(:,:,ch), lamdas(k), 0.6, 0, param);
        relerr(k) = norm(L-Xd(:,:,ch),'fro')/norm(Xd(:,:,ch),'fro');
        X_lq(:,:,k) = L;
    end   
    [~, mi] = min(relerr);
    X_Lq(:,:,ch) = X_lq(:,:,mi); 
end

figure(1);subplot(1,4,4);
imshow(uint8(X_Lq));
title(sprintf('Lq (best, q1=0.6, q2=0.0)\n RelErr=%.3f, PSNR=%.2f dB',norm(X_Lq(:) - Xd(:))/norm(Xd(:)),psnr(X(:),X_Lq(:))));
