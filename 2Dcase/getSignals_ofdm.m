function [fs1, fs2, cs, H, A1, A2, X0, B, y] = getSignals_ofdm(r, s, N1, N2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs 
%       r: number of point sources
%       s: dimension of subspace
%       N1, N2:length of 2D signal

%        
% Outputs        
%       fs1, fs2: locations of point sources
%       cs:amplititudes of point sources
%       H: unknow coefficients vector matrix for the subspace
%       A1,A2: Vandermonde matrix associated with the point sources
%       X0: the signal matrix
%       B: known subspace
%       y: measurements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs1 = rand(r,1);
fs2 = rand(r,1);


% dynamic_range= 20;
dynamic_range= 10;
cs = exp(-1i*2*pi*rand(r,1)) .* (1 + 10.^(rand(r,1).*(dynamic_range/20)));

H = zeros(s, r);
h = randn(s,1);
for kk = 1:r
    H(:,kk) = h/norm(h);
end

A1 = zeros(r, N1);
A2 = zeros(r, N2);
X0 = zeros(s, N1*N2);
for kk  = 1:r
    tmp1 = amplitude(fs1(kk), N1);
    A1(kk,:) = tmp1';
end

for nn = 1:N2
    for kk = 1:r   
        tmp2 = amplitude(fs2(kk), N2);
        A2(kk,:) = tmp2';
    end
    Cs = diag(cs.*A2(:,nn)); %r x r diag matrix
    X0(:,(nn-1)*N1+1:nn*N1) = H * Cs * A1;
end


% get the observation
B = zeros(N1*N2,s);
for i = 1 : N1*N2
    sigma = randn(1);
    B((i),:) = amplitude(sigma, s);
end


y = zeros(N1*N2,1);

for i = 1:N1*N2
    y(i,1) =  B(i,:) * X0(:,i);
end

end

%%
function a = amplitude(fs, n)
N = 0:1:n-1;
a = exp(1i*2*pi*fs.*N);
end
