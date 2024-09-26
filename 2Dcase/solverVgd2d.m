function [si,iter,X,ratio,mg,timer,cond,err] = ...
    solverVgd2d(y,B,N1,N2,r,s,maxit,trace,X0,tol_1,tol_2,tol_3,t_end,test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs 
%     y: linear measurements.
%     B: low-dimensional matrix.
%     N1,N2: length of the signal to be recovered.
%     r: model order of the spectrally sparse signal.
%     s: the dimension of subspace.
%    mu: fixed stepsize (0.5).
% maxit: maximum number of iterations.

% trace: display relative change in signal and function value per iteraiton
%        if not zero.
%    X0: target matrix.
% tol_1: stopping criterion measuring relative change in the iterations. 
% tol_2: stopping criterion measuring relative change in gradient magnitude.
% tol_3: stopping criterion mearsuring recovery error
% t_end: stop when maximum run time is reached
%   test : 0, run with stopping criterion;1, run until the maximum
%          iteration is reached.

%Outputs
%    si: success indicator, 1 if convergence is achieved, 0 otherwise.
%  iter: iteration number at convergence.
%     X: recovered matrix.
% ratio: relative change in signal per iteration.
%    fv: function value per iteration.
%    mg: gradient magnitude per iteration. 
%  step: stepsize.
% ui,vi: indicator of projection.
%  cond: condition number of X0.
%   err: relative error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
tic
time = zeros(maxit,1);

si = 0;
ratio = zeros(maxit,1);
mg = zeros(maxit,1);
err = zeros(maxit,1);


%% Spectral initialization
if mod(N1,2) == 0
    N11 = N1/2;
    DD1 = [1:N11 N11 N11-1:-1:1]';
else
    N11 = (N1+1)/2;
    DD1 = [1:N11 N11-1:-1:1]';
end
N12 = N1+1-N11;
D1 = sqrt(DD1);
dy1 = repmat(D1',N2,1);
dy1 = dy1(:)';

if mod(N2,2) == 0
    N21 = N2/2;
    DD2 = [1:N21 N21 N21-1:-1:1]';
else
    N21 = (N2+1)/2;
    DD2 = [1:N21 N21-1:-1:1]';
end
N22 = N2+1-N21;
D2 = sqrt(DD2);
dy2 = repmat(D2',N1,1);
dy2 = dy2(:)';

L0 = B'*diag(y)*eye(N1*N2); %%A^*(y)(s*n^2)
L00 = zeros(s*N11,N12*N2); %%(sn1*n2n)
HL0 = zeros(s*N11*N21, N12*N22); %%GDA^*(y)
HL = zeros(s*N11,N12);
for i = 1:N2
    for j1 = 1:N11
        for j2 = 1:N12
            row_idx = (j1-1)*s+1:j1*s;
            HL(row_idx, j2) = L0(:, j1+j2-1);
        end 
    end
    L00(:,(i-1)*N12+1:i*N12) = HL;
end
for j1 = 1:N21
    for j2 = 1:N22
        row_idx = (j1-1)*s*N11+1:j1*s*N11;
        col_idx = (j2-1)*N12+1:j2*N12;
        HL0(row_idx, col_idx) = L00(:, (j1+j2-2)*N12+1:(j1+j2-1)*N12);
    end
end

%% best r-rank approximation
[U,sig,V] = svd(HL0);
sig = diag(sig);
sr = sig(1:r);
U0 = U(:,1:r);
V0 = V(:,1:r);

%% GD(X_0)
X00 = zeros(s*N11,N12*N2); %%(sn1*n2n)
HX0 = zeros(s*N11*N21, N12*N22); %%GDA^*(y)
HX = zeros(s*N11,N12);
for i = 1:N2
    for j1 = 1:N11
        for j2 = 1:N12
            row_idx = (j1-1)*s+1:j1*s;
            HX(row_idx, j2) = X0(:, j1+j2-1);
        end 
    end
    X00(:,(i-1)*N12+1:i*N12) = HX;
end
for j1 = 1:N21
    for j2 = 1:N22
        row_idx = (j1-1)*s*N11+1:j1*s*N11;
        col_idx = (j2-1)*N12+1:j2*N12;
        HX0(row_idx, col_idx) = X00(:, (j1+j2-2)*N12+1:(j1+j2-1)*N12);
    end
end

%% calculate the condition number of H(X_0)
[~,sigx,~] = svd(HX0);
sigx = diag(sigx);
srx = sigx(1:r);
%U0x = Ux(:,1:r);
%V0x = Vx(:,1:r);
cond = srx(1)/srx(r);


%% set the initialization
Zu = U0*diag(sqrt(sr));
Zv = V0*diag(sqrt(sr));


Z_hat = Zu*Zv';

Z0u = Zu;
Z0v = Zv;
s0 = max(norm(Z0v),norm(Z0u));
X = Gstar2d(Z_hat,s,N1,N2)*diag(1./dy2);
%% Successive Iteration
for iter = 1:maxit
   X_old = X;

%     W = X-Astar(Aop(X,B)-y,B);
% 
%     gu = Zu*(Zv'*Zv)-G(W*diag(dy2),s,N1,N2)*Zv;
%     
% %     gv1 = Zv*(Zu'*Zu) - (G(X*diag(dy2),s,N1,N2))'*Zu;
% %     gv2 = (G(Astar(Aop(X*diag(dy2),B),B),s,N1,N2))'*Zu;
% %     gv3 = (G(Astar(diag(dy2)*y,B),s,N1,N2))'*Zu;
% %     gv = (gv1 + (gv2 - gv3)) ;
%     gv =  Zv*(Zu'*Zu) -G(W*diag(dy2),s,N1,N2)'*Zu;
    
    gu1 = Zu*(Zv'*Zv) - G(Gstar2d(Zu*Zv',s,N1,N2),s,N1,N2)*Zv;
    gu2 = G(Astar(Aop(Gstar2d(Zu*Zv',s,N1,N2),B),B),s,N1,N2)*Zv;
    gu3 = G(Astar(diag(dy2)*y,B),s,N1,N2)*Zv;
    gu = (gu1 + (gu2 - gu3)); 
    
    gv1 = Zv*(Zu'*Zu) - (G(Gstar2d(Zu*Zv',s,N1,N2),s,N1,N2))'*Zu;
    gv2 = (G(Astar(Aop(Gstar2d(Zu*Zv',s,N1,N2),B),B),s,N1,N2))'*Zu;
    gv3 = (G(Astar(diag(dy2)*y,B),s,N1,N2))'*Zu;
    gv = (gv1 + (gv2 - gv3)) ;
    
    mg(iter) = sqrt(norm(gu(:))^2+norm(gv(:))^2);
%     fv(iter) = fvalue(Zu,Zv,s,y,N1,N2,dy2,lambda,B);
    
    opt = 0;
    if opt == 1 %% backtrack line search
        
        k = 1;
        
        mu = s0^2; %% initial stepsize for line search
        ls = 0; %% the number of line search
        
        while k
            
            Zu_update = Zu - mu/s0^2*gu;
            Zv_update = Zv - mu/s0^2*gv;
            
            value_update = fvalue(Zu_update,Zv_update,s,y,N1,N2,dy2,lambda,B);
            
            if value_update < fv(iter)-1/50*(mu/s0^2)*(mg(iter)^2)
                k = 0;
            else
                mu = mu/2;
                ls = ls+1;
            end
            
        end
        
        step(iter) = 1/(2^ls);
        
    elseif opt == 0
        
        mu = 0.5;
        step(iter) = mu/s0^2;
        
    end
   
        
        Zu = Zu - mu/s0^2*gu;
        Zv = Zv - mu/s0^2*gv;
        
    MZ = Zu * Zv';
    X = Gstar2d(MZ,s,N1,N2)*diag(1./dy2);
    ratio(iter) = norm(X-X_old)/norm(X_old);
    err(iter) = norm(X(:)-X0(:))/norm(X0(:));
    
    if trace
        fprintf('Iteration %4d: relative.change = %.10f, mg = %.10f, err = %.10f \n',iter,ratio(iter),mg(iter),err(iter))
    end
    time(iter) = toc;
    timer = time(1:iter);
    if test==0
    if ratio(iter) < tol_1 || mg(iter) < tol_2 || err(iter) < tol_3 || toc>t_end
        si = 1;
        ratio = ratio(1:iter);
%         fv = fv(1:iter);
        mg = mg(1:iter);
%         step = step(1:iter); 
        err = err(1:iter);
%         ui = ui(1:iter);
%         vi = vi(1:iter);
        %t = t(1:iter);
        return; 
    end
    end

end   
end
% 
% %% G operator
% function g = G(X)
% [s,n]  = size(X);
% if mod(n,2) == 0
%     n1 = n/2;
%     DD = [1:n1 n1 n1-1:-1:1].';
% else
%     n1 = (n+1)/2;
%     DD = [1:n1 n1-1:-1:1].';
% end
% n2 = n+1 - n1;
% D = sqrt(DD);
% 
% dx = X*diag(1./D);
% g = zeros(s*n1, n2);
% for j1 = 1:n1
%     for j2 = 1:n2
%         row_idx = (j1-1)*s+1:j1*s;
%         g(row_idx, j2) = dx(:, j1+j2-1);
%     end
% end
% end
% 
% %% G^* operator
% function  gstz = gstar(Z,s,D)
% [N,n2] = size(Z);
% n1 = N/s;
% gstz = zeros(s,n1+n2-1);
% for i = 1:(n1+n2-1)
%     if i <= n1
%         comp = zeros(s,1);
%         for j = 1:i
%             comp = comp + Z(s*(i-j)+1:s*(i+1-j),j);
%         end
%         gstz(:,i) = 1/D(i)*comp;
%     else
%         comp = zeros(s,1);
%         for j = i+1-n1:n2
%             comp = comp + Z(s*(i-j)+1:s*(i+1-j),j);
%         end
%         gstz(:,i) = 1/D(i)*comp;
%     end
% end
% end

%% G operator
function g = G(X,s,N1,N2)
%[s,nn] = size(X);
%n = sqrt(nn);
if mod(N1,2) == 0
    N11 = N1/2;
    DD1 = [1:N11 N11 N11-1:-1:1].';
else
    N11 = (N1+1)/2;
    DD1 = [1:N11 N11-1:-1:1].';
end
N12 = N1+1-N11;
D1 = sqrt(DD1);

if mod(N2,2) == 0
    N21 = N2/2;
    DD2 = [1:N21 N21 N21-1:-1:1].';
else
    N21 = (N2+1)/2;
    DD2 = [1:N21 N21-1:-1:1].';
end
N22 = N2+1-N21;
D2 = sqrt(DD2);
XX = zeros(s*N11,N12*N2);
for i = 1:N2
    GX = zeros(s*N11,N12);
    for j1 = 1:N11
        for j2 = 1:N12
            row_idx = (j1-1)*s+1:j1*s;
            GX(row_idx, j2) = X(:, (i-1)*N1 + (j1+j2-1));
        end 
    end
    XX(:,(i-1)*N12+1:i*N12) = 1/D2(i)*GX;
end
%dx = XX*diag(1./D);
g = zeros(s*N11*N21, N12*N22);
for j1 = 1:N21
    for j2 = 1:N22
        row_idx = (j1-1)*s*N11+1:j1*s*N11;
        col_idx = (j2-1)*N12+1:j2*N12;
        g(row_idx, col_idx) = XX(:, (j1+j2-2)*N12+1:(j1+j2-1)*N12);
    end
end
end


%% G^* operator
function gst = Gstar2d(Z,s,N1,N2) %Z sn11*n12n2
if mod(N1,2) == 0
    N11 = N1/2;
    DD1 = [1:N11 N11 N11-1:-1:1].';
else
    N11 = (N1+1)/2;
    DD1 = [1:N11 N11-1:-1:1].';
end
N12 = N1+1-N11;
D1 = sqrt(DD1);
if mod(N2,2) == 0
    N21 = N2/2;
    DD2 = [1:N21 N21 N21-1:-1:1].';
else
    N21 = (N2+1)/2;
    DD2 = [1:N21 N21-1:-1:1].';
end
D2 = sqrt(DD2);
gx = gstar_2d(Z,N1,N2,s);
gst = zeros(s,N1*N2);
for i = 1:N2
    gst(:,(i-1)*N1+1:i*N1) = D2(i)*gstar(gx(:,(i-1)*N12+1:i*N12),s);
end
end

function  gstz_2d = gstar_2d(Z,N1,N2,s) %Z sn11n21*n12n22
if mod(N1,2) == 0
    N11 = N1/2;
    DD1 = [1:N11 N11 N11-1:-1:1].';
else
    N11 = (N1+1)/2;
    DD1 = [1:N11 N11-1:-1:1].';
end
N12 = N1+1-N11;
D1 = sqrt(DD1);

if mod(N2,2) == 0
    N21 = N2/2;
    DD2 = [1:N21 N21 N21-1:-1:1].';
else
    N21 = (N2+1)/2;
    DD2 = [1:N21 N21-1:-1:1].';
end
N22 = N2+1-N21;
D2 = sqrt(DD2);
gstz_2d = zeros(s*N11,(N21+N22-1)*N12);
for i = 1:(N21+N22-1)
    if i <= N21
        comp = zeros(s*N11,N12);
        for j = 1:i
            row_idx = s*N11*(i-j)+1:s*(i+1-j)*N11;
            col_idx = (j-1)*N12+1:j*N12;
            comp = comp+Z(row_idx,col_idx);
        end
        gstz_2d(:,(i-1)*N12+1:i*N12) = (1/D2(i)^2)*comp;
    else
        comp = zeros(s*N11,N12);
        for j = i+1-N21:N22
            row_idx = s*N11*(i-j)+1:s*(i+1-j)*N11;
            col_idx = (j-1)*N12+1:j*N12;
            comp = comp+Z(row_idx,col_idx);
        end
        gstz_2d(:,(i-1)*N12+1:i*N12) = (1/D2(i)^2)*comp;
    end
end
end

function gstz = gstar(Z,s) %Z sn11*n12
[N,n2] = size(Z);
n1 = N/s;
n = n1 + n2 - 1;
if mod(n,2) == 0
    N11 = n/2;
    DD = [1:N11 N11 N11-1:-1:1].';
else
    N11 = (n+1)/2;
    DD = [1:N11 N11-1:-1:1].';
end
D = sqrt(DD);
gstz = zeros(s,n1+n2-1);
for i = 1:(n1+n2-1)
    if i <= n1
        comp = zeros(s,1);
        for j = 1:i
            comp = comp + Z(s*(i-j)+1:s*(i+1-j),j);
        end
        gstz(:,i) = (1/D(i)^2)*comp;
    else
        comp = zeros(s,1);
        for j = i+1-n1:n2
            comp = comp + Z(s*(i-j)+1:s*(i+1-j),j);
        end
        gstz(:,i) = (1/D(i)^2)*comp;
    end
end
end

%% A operator
function y = Aop(X,B)
[~,n] = size(X);
y=zeros(n,1);
I = eye(n);
for i = 1:n
    y(i) = B(i,:)*X*I(:,i);
end    
end

%% A^* operator
function X = Astar(y,B)
X = B'*diag(y);
end
