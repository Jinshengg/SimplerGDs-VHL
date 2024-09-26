%% run a fxied 2D case via our algorithms, and plot the delay-doppler location via 2D MUSIC

clear all;
N1 = 17;
N2 = 13;
% r number of complex exponentials
r = 5;
% s dimension of subspace, which is less than N
s = 4;
maxit = 500;
t_end = 0.32;
% %% Genearate ofdm signal
% [fs1, fs2, cs, H, ~, ~, X0, B, ~] = getSignals_ofdm(r, s, N1, N2);
% 
% % load('fs.mat'); 
% % load('cs.mat');
% % load('H.mat');
% % load('B.mat');
% % load('y.mat');
% % load('X0.mat');
% 
% fs1 = sort(fs1);
% fs2 = sort(fs2);
% A1 = exp(-1i*2*pi*fs1.* (0:1:N1-1));
% A2 = exp(-1i*2*pi*fs2.* (0:1:N2-1));
% 
% X0 = zeros(s, N1*N2);
% 
% for nn = 1:N2
%     CsA2 = diag(cs.*A2(:,nn)); %r x r diag matrix
%     X0(:,(nn-1)*N1+1:nn*N1) = H * CsA2 * A1;
% end
% 
% y = zeros(N1*N2,1);
% 
% for i = 1:N1*N2
%     y(i,1) =  B(i,:) * X0(:,i);
% end

% load('./data/real2Dcase20231130T180639.mat');
load('real2Dcase.mat');

% 

% sigma_e = 0.15;
% e = randn(N1*N2,1);
% e = sigma_e*norm(y)*e/norm(e);
% y = y+e;

% solve by VGD-VHL (2D)

% tic;
% [~,~,X_pgd,~,~,~,~,time_pgd,cond1,err1] = solverPgd2d_fh(y,B,N1,N2,r,s,maxit,1/4,0,0,X0,1e-7,1e-5,1e-6,0);
% t_pgd = toc;
% err_pgd = norm(X_pgd(:)-X0(:))/norm(X0(:));

%% PGD-2D
tic;
[~,~,X_pgd,~,~,~,~,~,~,time_pgd,cond1,err1] = solverPgd2d(y,B,N1,N2,r,s,maxit,1/4,0,0,X0,1e-7,1e-5,1e-6,t_end,0);
t_pgd = toc;
err_pgd = norm(X_pgd(:)-X0(:))/norm(X0(:));

freq = music_2d(X_pgd,s,r,N1,N2);

% fs1_rec = sort(freq(1,:));
% fs2_rec = sort(freq(2,:));
fs1_rec = freq(1,:);
fs2_rec = freq(2,:);
%% recover rotations Hs
A_linSytem = zeros(N1*N2, s*length(fs1_rec));
for jj = 1:N2
    A_rec = diag(exp(-1i*2*pi*fs2_rec'.* (jj-1)))*exp(-1i*2*pi*fs1_rec'.* (0:1:N1-1));
    for kk = 1:N1
        A_linSytem((jj-1)*N1+kk,:) = kron(A_rec(:,kk).', B((jj-1)*N1+kk,:));
    end
end


Hs_ds = pinv(A_linSytem)*y;

amp = zeros(length(fs1_rec),1);
Hs_rec= zeros(s,length(fs1_rec));
for jj = 1:length(fs1_rec)
    idx = (jj-1)*s+1:jj*s;
    tmp = Hs_ds(idx);
    amp(jj) = norm(tmp);
    Hs_rec(:,jj) = tmp/norm( tmp);
%     fprintf('Inner : %f\n', abs(Hs_rec(:,jj)'*H(:,jj)));
end
[~,ind] = maxk(abs(amp),r);

figure
h1=stem3(fs1,fs2,abs(cs),'k--+','LineWidth',3,'MarkerSize',8);
hold on
h2=stem3(fs1_rec(ind)',fs2_rec(ind)',amp(ind),'r:o','LineWidth',3,'MarkerSize',11);
axis([0,1,0,1,0,5]);
grid on 
l0=legend('Original','PGD-VHL');
set(l0,'position',[0.45 0.80 0.25 0.14]);
% set(l0,'Location','north');

xlabel('$\phi$','interpreter','latex');
ylabel('$\psi$','interpreter','latex');
zlabel('magnitude');

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gca,'FontName','times new roman','FontSize',25 ,'Layer','top','LineWidth',3);
set(gca,'ztick', 0:2.5:5);
% print('-depsc','PGD_2D_freq.eps');
%% solve by VGD-VHL (2D)
% tic;
% [~,~,X_pgd,~,~,~,~,time_pgd,cond1,err1] = solverPgd2d_fh(y,B,N1,N2,r,s,maxit,1/4,0,0,X0,1e-7,1e-5,1e-6,0);
% t_pgd = toc;
% err_pgd = norm(X_pgd(:)-X0(:))/norm(X0(:));

tic;
[~,~,X_vgd,~,~,time_vgd,cond2,err2] = solverVgd2d(y,B,N1,N2,r,s,maxit,0,X0,1e-7,1e-5,1e-6,t_end,0);
t_vgd = toc;
err_vgd = norm(X_vgd(:)-X0(:))/norm(X0(:));

freq = music_2d(X_vgd,s,r,N1,N2);
fs1_rec = freq(1,:);
fs2_rec = freq(2,:);
%% recover rotations Hs
A_linSytem = zeros(N1*N2, s*length(fs1_rec));
for jj = 1:N2
    A_rec = diag(exp(-1i*2*pi*fs2_rec'.* (jj-1)))*exp(-1i*2*pi*fs1_rec'.* (0:1:N1-1));
    for kk = 1:N1
        A_linSytem((jj-1)*N1+kk,:) = kron(A_rec(:,kk).', B((jj-1)*N1+kk,:));
    end
end


Hs_ds = pinv(A_linSytem)*y;

amp = zeros(length(fs1_rec),1);
Hs_rec= zeros(s,length(fs1_rec));
for jj = 1:length(fs1_rec)
    idx = (jj-1)*s+1:jj*s;
    tmp = Hs_ds(idx);
    amp(jj) = norm(tmp);
    Hs_rec(:,jj) = tmp/norm( tmp);
%     fprintf('Inner : %f\n', abs(Hs_rec(:,jj)'*H(:,jj)));
end
[~,ind] = maxk(abs(amp),r);

figure
h1=stem3(fs1,fs2,abs(cs),'k--+','LineWidth',3,'MarkerSize',8);
hold on
h2=stem3(fs1_rec(ind)',fs2_rec(ind)',amp(ind),'r:o','LineWidth',3,'MarkerSize',11);
axis([0,1,0,1,0,5]);
grid on 
l0=legend('Original','VGD-VHL');
set(l0,'position',[0.45 0.80 0.25 0.14]);
% set(l0,'Location','north');

xlabel('$\phi$','interpreter','latex');
ylabel('$\psi$','interpreter','latex');
zlabel('magnitude');

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gca,'FontName','times new roman','FontSize',25 ,'Layer','top','LineWidth',3);
set(gca,'ztick', 0:2.5:5);

% print('-depsc','VGD_2D_freq.eps');


%% solve by ScalGD-VHL (2D)
% tic;
% [~,~,X_pgd,~,~,~,~,time_pgd,cond1,err1] = solverPgd2d_fh(y,B,N1,N2,r,s,maxit,1/4,0,0,X0,1e-7,1e-5,1e-6,0);
% t_pgd = toc;
% err_pgd = norm(X_pgd(:)-X0(:))/norm(X0(:));

tic;
[~,~,X_scal,~,~,time_scalgd,cond2,err3] = solverScaledgd2d(y,B,N1,N2,r,s,maxit,0,X0,1e-7,1e-5,1e-6,t_end,0);
t_scalgd = toc;
err_scalgd = norm(X_scal(:)-X0(:))/norm(X0(:));

freq = music_2d(X_scal,s,r,N1,N2);

fs1_rec = freq(1,:);
fs2_rec = freq(2,:);
%% recover rotations Hs
A_linSytem = zeros(N1*N2, s*length(fs1_rec));
for jj = 1:N2
    A_rec = diag(exp(-1i*2*pi*fs2_rec'.* (jj-1)))*exp(-1i*2*pi*fs1_rec'.* (0:1:N1-1));
    for kk = 1:N1
        A_linSytem((jj-1)*N1+kk,:) = kron(A_rec(:,kk).', B((jj-1)*N1+kk,:));
    end
end


Hs_ds = pinv(A_linSytem)*y;

amp = zeros(length(fs1_rec),1);
Hs_rec= zeros(s,length(fs1_rec));
for jj = 1:length(fs1_rec)
    idx = (jj-1)*s+1:jj*s;
    tmp = Hs_ds(idx);
    amp(jj) = norm(tmp);
    Hs_rec(:,jj) = tmp/norm( tmp);
%     fprintf('Inner : %f\n', abs(Hs_rec(:,jj)'*H(:,jj)));
end
[~,ind] = maxk(abs(amp),r);

figure
h1=stem3(fs1,fs2,abs(cs),'k--+','LineWidth',3,'MarkerSize',8);
hold on
h2=stem3(fs1_rec(ind)',fs2_rec(ind)',amp(ind),'r:o','LineWidth',3,'MarkerSize',11);
axis([0,1,0,1,0,5]);
grid on 
l0=legend('Original','ScalGD-VHL');
set(l0,'position',[0.45 0.80 0.25 0.14]);
% set(l0,'Location','north');

xlabel('$\phi$','interpreter','latex');
ylabel('$\psi$','interpreter','latex');
zlabel('magnitude');

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gca,'FontName','times new roman','FontSize',25 ,'Layer','top','LineWidth',3);
set(gca,'ztick', 0:2.5:5);
