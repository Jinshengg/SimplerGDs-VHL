function [si,iter,X,ratio,fv,mg,step,timer,cond,err] = ...
    solverPgd_fh(y,B,n,r,s,maxit,lambda,opt,trace,X0,tol_1,tol_2,tol_3,test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs 
%     y: linear measurements.
%     B: low-dimensional matrix.
%     n: length of the signal to be recovered.
%     r: model order of the spectrally sparse signal.
%     s: the dimension of subspace.
% maxit: maximum number of iterations.
%lambda: the tuning parameter of balance term. (1/4)
%   opt: options for stepsize. 1 (line search), 0 (fixed).
% trace: display relative change in signal and function value per iteraiton
%        if not zero.
%    X0: target matrix.
% tol_1: stopping criterion measuring relative change in the iterations. 
% tol_2: stopping criterion measuring relative change in gradient magnitude.
% tol_3: stopping criterion measuring recovery error.
%   test : 0, run with stopping criterion;1, run until the maximum
%          iteration is reached.
% 

%
%Outputs
%    si: success indicator, 1 if convergence is achieved, 0 otherwise.
%  iter: iteration number at convergence.
%     X: recovered matrix.
% ratio: relative change in signal per iteration.
%    fv: function value per iteration.
%    mg: gradient magnitude per iteration. 
%  step: stepsize.
%  timer:run time
% ui,vi: indicator of projection.
%  cond: condition number of X0.
%   err: relative error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = zeros(maxit,1);
si = 0;
ratio = zeros(maxit,1);
fv = zeros(maxit,1);
mg = zeros(maxit,1);
step = zeros(maxit,1);
err = zeros(maxit,1);


ui = zeros(maxit,1);
vi = zeros(maxit,1);
%% Spectral initialization
if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1]';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1]';
end
n2 = n+1 - n1;
D = sqrt(DD);

L0 = B'*diag(y)*eye(n); %%A^*(y)

% HL0 = zeros(s*n1, n2); %%GDA^*(y)
% for j1 = 1:n1
%     for j2 = 1:n2
%         row_idx = (j1-1)*s+1:j1*s;
%         HL0(row_idx, j2) = L0(:, j1+j2-1);
%     end
% end

%% best r-rank approximation
[U, sig, V] = svds(@(v,tflag)vecHankelVecMul(L0, v, tflag, n, n1, n2, s), [s*n1, n2], r);
sig = diag(sig);
sr = sig(1:r);
U0 = U(:,1:r);
V0 = V(:,1:r);

%% GD(X_0)
% HX0 = zeros(s*n1, n2);
% for j1 = 1:n1
%     for j2 = 1:n2
%         row_idx = (j1-1)*s+1:j1*s;
%         HX0(row_idx, j2) = X0(:, j1+j2-1);
%     end
% end

%% calculate the condition number of H(X_0)
[~, sigx, ~] = svds(@(v,tflag)vecHankelVecMul(X0, v, tflag, n, n1, n2,s), [s*n1, n2], r);
sigx = diag(sigx);
srx = sigx(1:r);
cond = srx(1)/srx(r);

%% compute the projection parameter
UU = U0.*conj(U0);
UUb = reshape(UU',[r,s,n1]);
U0_inf = max(sqrt(sum(UUb,[1,2])));
V0_inf = max(sqrt(sum(V0.*conj(V0),2)));
param = max(U0_inf,V0_inf)*sqrt(sr(1));

%% set the initialization
Zu = U0*diag(sqrt(sr));
Zv = V0*diag(sqrt(sr));

%% projections
Zub = reshape(Zu',[r,s,n1]);
Zusb = reshape((Zu.*conj(Zu))',[r,s,n1]);
Zusb_sum = sqrt(sum(Zusb,[1,2]));
Zu_inf = max(Zusb_sum);

if Zu_inf > param %% Projection of Zu
    ub_ind = find(Zusb_sum > param);
    Zub(:,:,ub_ind) = Zub(:,:,ub_ind)./Zusb_sum(:,:,ub_ind)*param;
end

Zu = reshape(Zub,[r,s*n1])';
Zv_sum = sqrt(sum(Zv.*conj(Zv),2));
Zv_inf = max(Zv_sum);

if Zv_inf > param %% Projection of Zv
    vb_ind = find(Zv_sum > param);
    Zv(vb_ind,:) = Zv(vb_ind,:)./Zv_sum(vb_ind,:)*param;
end

% Z_hat = Zu*Zv';

Z0u = Zu;
Z0v = Zv;
s0 = max(norm(Z0v),norm(Z0u));
% X = gstar(Z_hat,s,D)*diag(1./D);
X = Gstar(Zu,Zv,n1,n2,s,D)*diag(1./D);

tic
%% Successive Iteration
for iter = 1:maxit
    err(iter) = norm(X(:)-X0(:))/norm(X0(:));
    X_old = X;
    
    W = Astar(Aop(X,B)-y,B)-X;
    
    gu1  = hankelmul_Zu(W,Zv,n1,n2,s,r);
    gu2 = Zu*(Zv'*Zv);
    
%     gu1 = (Zu*Zv')*Zv - G(gstar(Zu*Zv',s,D))*Zv;
%     gu2 = G(Astar(Aop(gstar(Zu*Zv',s,D),B),B))*Zv;
%     gu3 = G(Astar(diag(D)*y,B))*Zv;
    bgu = Zu*(Zu'*Zu - Zv'*Zv); %% the gradient of balance term
    gu = (gu1 + gu2 ) + lambda*bgu; 
    
    gv1 =  hankelmul_Zv(W,Zu,n1,n2,s,r);
    
    gv2 = Zv*(Zu'*Zu);
    
  
    bgv = Zv*(Zv'*Zv - Zu'*Zu); %% the gradient of balance term
    gv = (gv1 + gv2 ) + lambda*bgv;
    
    
%     mg(iter) = sqrt(sum(diag(gu'*gu+gv'*gv)));
%     fv(iter) = fvalue(Zu,Zv,s,y,B,D);
    mg(iter) =sqrt( norm(gv,'fro')^2+norm(gu,'fro')^2); 
    if opt == 1 %% backtrack line search
        
        k = 1;
        
        mu = s0^2; %% initial stepsize for line search
        ls = 0; %% the number of line search
        
        while k
            
            Zu_update = Zu - mu/s0^2*gu;
            Zv_update = Zv - mu/s0^2*gv;
            
            value_update = fvalue(Zu_update,Zv_update,s,y,B,D);
            
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
%         mu = 1/(s^2*r^2);
        step(iter) = mu/s0^2;
        
    end
   
        Zu = Zu - mu/s0^2*gu;
        Zv = Zv - mu/s0^2*gv;
        
    %% projections
    Zub = reshape(Zu',[r,s,n1]);
    Zusb = reshape((Zu.*conj(Zu))',[r,s,n1]);
    Zusb_sum = sqrt(sum(Zusb,[1,2]));
    Zu_inf = max(Zusb_sum);
    
    if Zu_inf > param %% Projection of Zu
        ui(iter) = 1;
        ub_ind = find(Zusb_sum > param);
        Zub(:,:,ub_ind) = Zub(:,:,ub_ind)./Zusb_sum(:,:,ub_ind)*param;
    end
    
    Zu = reshape(Zub,[r,s*n1])';
    Zv_sum = sqrt(sum(Zv.*conj(Zv),2));
    Zv_inf = max(Zv_sum);
    
    if Zv_inf > param %% Projection of Zv
        vi(iter) = 1;
        vb_ind = find(Zv_sum > param);
        Zv(vb_ind,:) = Zv(vb_ind,:)./Zv_sum(vb_ind,:)*param;
    end

    
%     MZ = Zu * Zv';
    X = Gstar(Zu,Zv,n1,n2,s,D)*diag(1./D);
%     X = gstar(MZ,s,D)*diag(1./D);
    ratio(iter) = norm(X-X_old)/norm(X_old);
    
    if trace
        fprintf('Iteration %4d: relative.change = %.10f, mg = %.10f, err = %.10f  t = %.10f \n',iter,ratio(iter),mg(iter),err(iter),toc)
    end
  time(iter) = toc;  
  timer = time (1:iter);      
  if (test ==0)
   if ratio(iter) < tol_1 || mg(iter) < tol_2 || err(iter) < tol_3
        si = 1;
        ratio = ratio(1:iter);

        mg = mg(1:iter);

        err = err(1:iter);
 
        return; 
   end
  end

end   
end

%% G operator
function g = G(X)
[s,n]  = size(X);
if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1].';
end
n2 = n+1 - n1;
D = sqrt(DD);

dx = X*diag(1./D);
g = zeros(s*n1, n2);
for j1 = 1:n1
    for j2 = 1:n2
        row_idx = (j1-1)*s+1:j1*s;
        g(row_idx, j2) = dx(:, j1+j2-1);
    end
end
end

%% G^* operator
function  gstz = gstar(Z,s,D)
[N,n2] = size(Z);
n1 = N/s;
gstz = zeros(s,n1+n2-1);
for i = 1:(n1+n2-1)
    if i <= n1
        comp = zeros(s,1);
        for j = 1:i
            comp = comp + Z(s*(i-j)+1:s*(i+1-j),j);
        end
        gstz(:,i) = 1/D(i)*comp;
    else
        comp = zeros(s,1);
        for j = i+1-n1:n2
            comp = comp + Z(s*(i-j)+1:s*(i+1-j),j);
        end
        gstz(:,i) = 1/D(i)*comp;
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

%% function value
function fv = fvalue(Zu,Zv,s,y,B,D)
Z = Zu*Zv';
% [N,n2] = size(Z);
% n1 = N/s;
% n = n1+n2-1;
fv = 1/2*norm(Z - G(gstar(Z,s,D)),'f')^2 + 1/2*norm(Aop(gstar(Z,s,D),B)-diag(D)*y)^2;
end

%% G^* fast dehankel operator
function  Gstz = Gstar(Zu,Zv,n1,n2,s,D)
    nd = n1+n2-1;
    Gstz = zeros(s,nd);
    n = 2^nextpow2(nd);
    zv = fft(conj(Zv),n);
    for i1 = 1:s
    row_ind = (i1-1)*n1+1:i1*n1;
    x = sum(ifft(fft(Zu(row_ind,:), n) .* zv), 2);
    Gstz(i1,:) = x(1:nd)./D;
   end
end

function Zu = hankelmul_Zu(W,Zv,n1,n2,s,r)
 nd = n1+n2 -1;
 n = 2^nextpow2(nd);
 Zu = zeros(n1*s,r);
 tempV = fft(flip(Zv),n);
 for i1 = 1:s
     row_ind = (i1-1)*n1+1:i1*n1;
     tempU = ifft(bsxfun(@times, tempV, fft(W(i1,:).', n)));
     Zu(row_ind,:) = tempU(n2:nd,:);
 end
end

function Zv = hankelmul_Zv(W,Zu,n1,n2,s,r)
 nd = n1+n2 -1;
 n = 2^nextpow2(nd);
 Zv= zeros(n2,r);
%  tempV = zeros(s,n2,r);
 for i1 = 1:s
     row_ind = (i1-1)*n1+1:i1*n1;
     tempV = ifft(bsxfun(@times, fft(flip(Zu(row_ind,:)),n), fft(W(i1,:)', n)));
     Zv = tempV(n1:nd,:)+Zv;    
 end
end


function z = vecHankelVecMul(L, v, tflag, nd, n1, n2,s)
% Hankel multiplication Z=H[L]v if tflag='notransp'; z=H[x]'v if tflag='transp'
    n = 2^nextpow2(nd);
   
    if strcmp(tflag, 'notransp')
         for i1 = 1:s
        row_ind = (i1-1)*n1+1:i1*n1;
        tempz = ifft(fft(flip(v), n) .* fft(L(i1,:).', n));
        z(row_ind,1) = tempz(n2:nd);
         end
    else
        z = zeros(n2,1);
          for i1 = 1:s
        row_ind = (i1-1)*n1+1:i1*n1;
        tempz = ifft(fft(flip(v(row_ind,:)), n) .* fft(L(i1,:)', n));
        z =tempz(n1:nd)+z;
          end
    end
end