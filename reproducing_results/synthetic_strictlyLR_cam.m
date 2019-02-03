clear all; clc; % close all;

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

figure(1);subplot(1,4,1); imshow(uint8(Y));
title(sprintf('Corrupted\n RelErr=%.3f, PSNR=%.2f dB',norm(Y - X,'fro')/norm(X,'fro'),psnr(X,Y)));


param.TOL =1e-16;
param.MAX_ITER = 2e3;

lamdas = logspace(-5, 2, 50);

% --- soft thresholding --------------------
parfor k = 1:length(lamdas)
    [L, S, out] = lq_lq_l2_bcd(Y, lamdas(k), 1, 1, param);
    relerr(k) = norm(L-X,'fro')/norm(X,'fro');
    X_l1(:,:,k) = L;
    S_l1(:,:,k) = S;
end
[RelErr_L1, mi] = min(relerr);
X_L1 = X_l1(:,:,mi); 
S_L1 = S_l1(:,:,mi); 

figure(1);subplot(1,4,2);
imshow(uint8(X_L1));
title(sprintf('Soft\n RelErr=%.3f, PSNR=%.2f dB',RelErr_L1,psnr(X,X_L1)));


% ----- Lq thresholding ----------
param.L0 = X_L1;
param.S0 = S_L1;
qs = 0:0.2:1;
for l1=1:length(qs)
    for l2=1:length(qs)
        
        t0=tic;
        parfor k = 1:length(lamdas)           
            [L, ~, ~] = lq_lq_l2_bcd(Y, lamdas(k), qs(l1), qs(l2), param);
            relerr(k) = norm(L-X,'fro')/norm(X,'fro');
            xx(:,:,k) = L;
        end

        [RelErrs(l1,l2), mi] = min(relerr); 
        x_Lq(:,l1,l2)  = reshape(xx(:,:,mi), [m*n,1]); 
        PSNR(l1,l2) = psnr(X, xx(:,:,mi));

        sprintf('BCD with q1=%.1f and q2=%.1f completed, elapsed time: %.1f seconds',qs(l1),qs(l2),toc(t0))
    end
end

v0 = min(min(PSNR));v1 = max(max(PSNR));
figure(2);
contourf(qs,qs,PSNR,[v0:10:v1]); colorbar; xlabel('q_2'); ylabel('q_1');
set(gca, 'CLim', [v0, v1]);

[w1, e1] = max(PSNR); [~, lo] = max(w1); ko = e1(lo);
figure(2);hold on;
plot(qs(lo),qs(ko),'r*')
hold off;

figure(1);subplot(1,4,4);
imshow(uint8(reshape(x_Lq(:,ko,lo),[m,n])));
title(sprintf('Lq (best, q1=%.1f, q2=%.1f)\n RelErr=%.3e, PSNR=%.2f dB',qs(ko),qs(lo),RelErrs(ko,lo),PSNR(ko,lo)));

figure(1);subplot(1,4,3);
imshow(uint8(reshape(x_Lq(:,1,1),[m,n])));
title(sprintf('Hard\n RelErr=%.3e, PSNR=%.2f dB',RelErrs(1,1),PSNR(1,1)));
