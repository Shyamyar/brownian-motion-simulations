clear;
clc;
n = 1000;
H = 0.5;
vec_no = 1;
cov = zeros(n+1);
for i=1:length(cov)
    for j=1:length(cov)
        cov(i,j)= 0.5*(i^(2*H)+j^(2*H)-abs(j-i)^(2*H));
    end
end

Z = normrnd(0,1/sqrt(n),vec_no,n);
M = chol(cov,'lower');
G = M*[zeros(vec_no,1),Z]';
% 
figure(1);
hold on
 
plot((0:n)/n,G)

if  H == 0.5
    plot((0:n)/n,[zeros(vec_no,1),cumsum(Z')'],'--')
end

xlabel('t'); ylabel('B_t');
hold off;

figure(2)
X = mvnrnd(zeros(1,n+1),cov,vec_no);
plot((0:n)/n,X);