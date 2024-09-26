function freq = music_2d(X,s,r,N1,N2)
if mod(N1,2) == 0
    N11 = N1/2;
else
    N11 = (N1+1)/2;
end
N12 = N1+1-N11;

if mod(N2,2) == 0
    N21 = N2/2;
else
    N21 = (N2+1)/2;
end
N22 = N2+1-N21;

XX = zeros(s*N11,N12*N2);
for i = 1:N2
    hX = zeros(s*N11,N12);
    for j1 = 1:N11
        for j2 = 1:N12
            row_idx = (j1-1)*s+1:j1*s;
            hX(row_idx, j2) = X(:, (i-1)*N1 + (j1+j2-1));
        end 
    end
    XX(:,(i-1)*N12+1:i*N12) = hX;
end
%dx = XX*diag(1./D);
HX = zeros(s*N11*N21, N12*N22);
for j1 = 1:N21
    for j2 = 1:N22
        row_idx = (j1-1)*s*N11+1:j1*s*N11;
        col_idx = (j2-1)*N12+1:j2*N12;
        HX(row_idx, col_idx) = XX(:, (j1+j2-2)*N12+1:(j1+j2-1)*N12);
    end
end

[U,D,V] = svd(HX);
d = diag(D);
Vn = V(:,r+1:end);
N11 = [0:1:N12-1]';
N21 = [0:1:N22-1]';

grid_size = 1e3;
f1 = [0:1/grid_size:1-1/grid_size];
f2 = [0:1/grid_size:1-1/grid_size];

for g1 = 1:grid_size 
    for g2 = 1:grid_size
        F = exp(1i*2*pi*f1(g1)*N11) * exp(1i*2*pi*f2(g2)*N21).';
        f = vec(F);
        %P(g1,g2) = abs((f'*Vn)*(f'*Vn)');
        Val(g1,g2) = norm(f'*Vn);
    end
end

[ind1,ind2] = find(abs(Val) < 7);
v1 = [];
v2 = [];
for t = 1:length(ind1)
    v = [ind1(t)-1,ind1(t)-1,ind1(t)-1,ind1(t),ind1(t),ind1(t)+1,ind1(t)+1,ind1(t)+1;
        ind2(t)-1,ind2(t),ind2(t)+1,ind2(t)-1,ind2(t)+1,ind2(t)-1,ind2(t),ind2(t)+1];
    
    for i = 1:length(v)
        if v(1,i) < 1
            v(1,i) = grid_size - v(1,i);
        end
        if v(1,i) > grid_size
            v(1,i) = 1;
        end
        if v(2,i)<1
            v(2,i) = grid_size - v(2,i);
        end
        if v(2,i) > grid_size
            v(2,i) = 1;
        end
    end
    
    for i=1:length(v)
        Val_edge(i) = Val(v(1,i), v(2,i));
    end
    
    if abs(Val(ind1(t),ind2(t))) < min(abs(Val_edge))
        v1 = [v1,ind1(t)];
        v2 = [v2,ind2(t)];
    end
end

for t = 1:length(v1)
    freq(1,t) = f1(v1(t));
    freq(2,t) = f2(v2(t));
end

end