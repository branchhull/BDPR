L = 3000;
N = 20;
K = 20;

rng(1)
% B = randn(L,K);
% B = B/sqrt(L);


idx2 = 1:L;
idx2 = idx2(1:K);
B = eye(L);
B = B(:,idx2); % sparse


C = randn(L,N);
C = C/sqrt(L);

ms = [ones(round(N/2),1);-ones(N-round(N/2),1)];
hs = [ones(round(K/2),1);-ones(K-round(K/2),1)];

MatA1 = fft(C);
MatA2 = fft(B);

yhat = abs(MatA1*ms).*abs(MatA2*hs);

[X1, X2] = PSD_BH(MatA1, MatA2, yhat.^2, 600, 50);

figure('position',[100 100 800 500])
subplot(211);
plot(ms);title('original m-signal');

subplot(212);
plot(real(X1(1,:)));title('recovered m-signal');
