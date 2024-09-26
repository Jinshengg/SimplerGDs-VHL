clear; close all; clc;
%% parameters 

n = 1024; % length of signal 
s = 4; %  dim of subspace 
r = 4; % number of point sources
maxit = 100; % maximum iteratioins
cond = 0;
% kappa = 0;
%% generate signal 
[fs, ~, ~, ~, X0, B, y] = getSignals_bdft_withsep(r, s, n,cond);

%% recovery 

[~,~,~,~,~,~,~,time_pgd,~,err1] = solverPgd_fh(y,B,n,r,s,maxit,1/4,0,0,X0,1e-7,1e-5,1e-6,1);

% tic
% [~,~,~,~,~,time_fiht,cond4,err4] = solverFIHT_fh(y,B,n,r,s,maxit,0,X0,1e-7,1e-5,1e-6,0);
% toc

[~,~,~,~,~,time_vgd,~,err2] = solverVgd_fh(y,B,n,r,s,maxit,0,X0,1e-7,1e-5,1e-6,1);


[~,~,~,~,~,time_scalgd,~,err3] = solverScaledgd_fh(y,B,n,r,s,maxit,0,X0,1e-7,1e-5,1e-6,1);


%% plot error versus time

fontsz = 24;

 clrs = {[.5,0,.5], [1,.5,0], [1,0,0], [0,.5,0], [0,0,1]};
mks = {'o', 'x', 'p', 's', 'd'};
figure('Position', [0,0,800,600], 'DefaultAxesFontSize', 24);
pace = 3;
ends = maxit; 

semilogy(1:pace:ends,err1(1:pace:ends),'Color',[1 0 1] , 'Marker', 'p', 'MarkerSize', 9,'LineWidth',1.5);
hold on;grid on;
semilogy(1:pace:ends,err2(1:pace:ends),'Color', 'b', 'Marker', 'o', 'MarkerSize', 11,'LineWidth',1.5);
hold on;grid on;
semilogy(1:pace:ends,err3(1:pace:ends),'Color', 'r', 'Marker', 'x', 'MarkerSize', 11,'LineWidth',1.5);
hold on;grid on;
% set(gca,'ytick', -32:4:0, 'xtick', 0:200:2000);
% set(gca,'ytick', -32:4:0, 'xtick', 0:20:200);
title(sprintf('$s = %d$, $r = %d$',s,r),'interpreter','latex','fontsize', fontsz);
xlabel('Iteration','interpreter','latex','fontsize', fontsz);
ylabel('Relative error','interpreter','latex','fontsize', fontsz);
% legend({'$pgd$','$vgd$'},'interpreter','latex','FontSize',14);
set(gca,'FontName','times new roman','FontSize',24,'Layer','top','linewidth',1.5);


legend({'PGD-VHL','VGD-VHL','ScalGD-VHL'},'interpreter','latex','FontSize',18);

