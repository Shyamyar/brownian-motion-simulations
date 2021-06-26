clear;
clc;
n = 100;
% H = 0.9;
% n = [50,100,250,500,1000,1250,1500,1750,2000];
H = [0.1,0.3,0.5,0.7,0.9];
vec_no = 1; % no. of random vectors
size_tot = (n+1)^2;
t = zeros(size_tot,2);
pos1 = 0;

%% Time Field
range = [0.01,1:n];
for i = range
    for j = range
        t(pos1+1,1) = i/n;
        t(pos1+1,2) = j/n;
        pos1 = pos1 + 1;
    end
end

%% Time norm
t_norm = zeros(1,size_tot);
for i = 1:size_tot
    t_norm(i) = norm(t(i,:));
end

%% Time difference norm
tic;
t_diff = zeros(size_tot,size_tot);
for i = 1:size_tot
    for j = 1:size_tot
    t_diff(i,j) = norm(t(i,:)-t(j,:));
    end
end
toc
%% Covariance Matrix with Max n
tic;
cov_total = zeros(size_tot,size_tot,length(H));
for f = 1:length(H)
    i = t_norm;
    j = i.';
    cov_total(:,:,f) = bsxfun(@(i,j) 0.5.*(i.^(2*H(f))+j.^(2*H(f))), i, j);
    cov_total(:,:,f) = cov_total(:,:,f) - (0.5*(t_diff.^(2*H(f))));
%     for i=1:length(t)
%         for j=1:length(t)
%             cov_total(i,j,f)= 0.5*(t_norm(i)^(2*H(f))+t_norm(j)^(2*H(f))-norm(t(i,:)-t(j,:))^(2*H(f)));
%         end
%     end
end
toc
%% Computation for Different n and H
Z = normrnd(0,1,vec_no,size_tot);

for f = 1:length(H)
    cov = cov_total(:,:,f);
    
    %% Cholesky Method
    M = chol(cov,'lower');
    B_chol = M*Z';

    %% Direct Method for H=0.5
    if  H(f) == 0.5
        B_direct = sqrt(n+1).\cumsum(Z');
    end

    %% Built in Multivariate Package (mvnrnd)
    B_mv = mvnrnd(zeros(1,size_tot),cov,vec_no);

    %% Plots
    time_space = linspace(0,1,n+1);
    array = [0,0.2,0.4,0.6,0.8,1.0];
    str = strings(1,n+1);
    for i=1:length(array)
        c = find(time_space==array(i));
        str(c) = string(array(i));
    end
    
    figure(f);
    h1 = heatmap(reshape(B_chol,[n+1,n+1]),'Colormap', parula(20), 'ColorbarVisible', 'on', 'XLabel', 't', 'YLabel', 't');
    title = sprintf('Levy Brownian Field for H = %.2f', H(f));
    grid off;
    h1.Title = title;
    h1.XDisplayLabels = str;
    h1.YDisplayLabels = str;
    set(struct(h1).NodeChildren(3), 'XTickLabelRotation', 0);

end