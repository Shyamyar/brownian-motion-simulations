clear;
clc;
n = 100;
H = [0.5,0.9];
% n = [50,100,250,500,1000,1250,1500,1750,2000];
% H = [0.1,0.3,0.5,0.7,0.9];
vec_no = n; % no. of random vectors
T = zeros(3,length(H),length(n)); % Time for Three Methods over different H's and n's

%% Covariance Matrix with Max n
cov_total = zeros(max(n)+1,max(n)+1,length(H));
for f = 1:length(H)
    c = (1/(max(n)))*[0.0001,1:(length(cov_total)-1)];
    y = c.';
    cov_total(:,:,f) = bsxfun(@(c,y) 0.5.*(c.^(2*H(f))+y.^(2*H(f))-abs(y-c).^(2*H(f))), c, y);
%     for i=c
%         for j=c
%             cov_total(find(c==i),find(c==j),f)= 0.5*(i^(2*H(f))+j^(2*H(f))-abs(j-i)^(2*H(f)));
%         end
%     end
end

%% Computation for Different n and H
for g = 1:length(n)
    Z = normrnd(0,1,vec_no,n(g));
    
    for f = 1:length(H)
        cov = cov_total(1:n(g)+1,1:n(g)+1,f);
        
        %% Cholesky Method
        tic;        
        M = chol(cov,'lower');
        B_chol = M*[zeros(vec_no,1),Z]';
        T(1,f,g) = toc;

        %% Direct Method for H=0.5
        tic;
        if  H(f) == 0.5
            B_direct = [zeros(vec_no,1),sqrt(n(g)).\cumsum(Z')'];
            T(2,f,g) = toc;
        end

        %% Built in Multivariate Package (mvnrnd)
        tic;
        B_mv = mvnrnd(zeros(1,n(g)+1),cov,vec_no);
        T(3,f,g) = toc;
        
        %% Plots
        figure(2*f-1);
        hold on
        plot((0:n(g))/n(g),B_chol)
        if  H(f) == 0.5 && vec_no == 1
            plot((0:n(g))/n(g),B_direct','--')
            legend('Chol','Direct(0.5)')
        elseif vec_no == 1
            legend('Chol')
        end
        xlabel('t'); ylabel('B_t');
        plot_title = sprintf('Fractional Brownian Motion, H = %.2f by Cholesky', H(f));
        title(plot_title)
        hold off;

        figure(2*f);
        plot((0:n(g))/n(g),B_mv);
        xlabel('t'); ylabel('B_t');
        plot_title = sprintf('Fractional Brownian Motion, H = %.2f by mvnrnd', H(f));
        title(plot_title)
        
    end
end

%% Plot times
%%
legends = {'Cholesky','Direct(H=0.5)','MV\_Package'};
figure(2*f+1);
plot(H,(reshape(T(:,:,1),[3,f])'),'o') % Time vs H
xlabel('H'); ylabel('Time');
legend(legends)

%%
slope = zeros(1,length(H));
interc = zeros(1,length(H));
start = 4;
for H_no=1:length(H)
%     H_no = 5;
    method_no = 1;    
    method = string(legends{method_no});
    
    P = polyfit(log(n(start:end)),log(reshape(T(method_no,H_no,start:end),[g-(start-1),1]))',1);
    slope(H_no) = P(1);
    interc(H_no) = P(2);
    yfit = P(1)*(log(n(start:end)))+P(2);
    yfit2 = 3*(log(n(start:end)))-(3*log(n(start))-yfit(1));
    
%     figure(2*H_no-1);
%     plot(n,smooth(reshape(T(method_no,H_no,:),[g,1])),'b-*')
%     xlabel('n'); ylabel('Time');    
%     plot_title = sprintf('Time Required for Method: %s with H: %.1f n: %d-%d', method, H(H_no), 50, n(g));
%     title(plot_title)
%     
%     figure(2*H_no);
%     scatter(log(n),log(reshape(T(method_no,H_no,:),[g,1])),'*')
%     hold on;
%     plot(log(n(start:end)),yfit,'r-.')
%     plot(log(n(start:end)),yfit2,'g-.')
%     xlabel('Log(n)'); ylabel('Log(Time)');
%     plot_title = sprintf('Time Required for Method: %s with H: %.1f n: %d-%d', method, H(H_no), 50, n(g));
%     legend('Actual Time Points','Trend Line','Line, Slope = 3')
%     title(plot_title)
end
slope_H = [H;slope]