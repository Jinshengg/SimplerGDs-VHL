function [fs, cs, H, A, X0, B, y] = getSignals_bdft_withsep(r, s, n,kappa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs 
%       r: number of point sources
%       s: dimension of subspace
%       n:length of signal
% kappa : conditon number of the target matrix
%           0, means random strength within dynamic range, resulting in a
%              the condition number of target Hankel lifted  matrix is random
%       others, the condition number of the target Hankel lifted matrix is  kappa
%        
% Outputs        
%       fs: locations of point sources
%       cs:amplititudes of point sources
%       H: unknow coefficients vector matrix for the subspace
%       A: Vandermonde matrix associated with the point sources
%       X0: the signal matrix
%       B:known subspace
%       y: measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = rand(r,1);

if r == 1
    sep = 1;
else
    fs_sort = sort(fs);
    sep = min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
end

while(sep<16/n) %primal 1/n
    fs = rand(r,1);
    if r == 1
        sep = 1;
    else
        fs_sort = sort(fs);
        sep = min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
    end
end

if (kappa==0)
dynamic_range= 20;
cs = exp(-1i*2*pi*rand(r,1)) .* (1 + 10.^(rand(r,1).*(dynamic_range/20)));
else
    
cs = exp(-1i*2*pi*rand(r,1)) .* linspace(1,kappa,r)';
% cs = linspace(1,kappa,r)'; % no random phase
end
% dynamic_range= 10;

Cs = diag(cs); %J x J diag matrix

H = zeros(s, r);
for kk = 1:r
    H(:,kk) = randn(s,1);
    H(:,kk) = H(:, kk)/norm(H(:,kk));
end

A = zeros(r, n);
for kk  = 1:r
    tmp = amplitude(fs(kk), n);
    A(kk,:) = tmp';
end
X0 = H * Cs * A;



% get the observation

% Fourier samples
F = fft(eye(s));
B = zeros(n,s);
for i = 1 : n
    B(i,:) = F(randi(s),:);
end


y = zeros(n,1);

for i = 1:n
    y(i,1) =  B(i,:) * X0(:,i);
end

end

%%
function a = amplitude(fs, n)
N = 0:1:n-1;
a = exp(1i*2*pi*fs.*N);
end
