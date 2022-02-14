clc
clear
close all

%% create samples:
%%
m = [1 1; 5 5; 9 1]';
S = [1 0.4; 0.4 1; 
    1 -0.6; -0.6 1; 
    1 0; 0 1];
N = 500;

randn('seed', 0);
rng('default') % For reproducibility
% generate data
[X, y] = generate_gauss_classes(m,S,N);
plot_data(X,y,m);
title('generated data');
legend('Gaussian 2', '', 'Gaussian 1', 'Gaussian 3');

data = X';

figure
plot(data(:,1),data(:,2),'k*');

%% q=2
% m = 2; % fuzzier

%% q=3
m = 3; % fuzzier

cluster_n = 3; %類別數量
iter = 100; % 疊代停止次數

num_data = size(data,1);%樣本個數
num_d = size(data,2);%樣本維度

% --初始化隸屬值uij，且含限制條件
U = rand(cluster_n,num_data);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n,1),:);

%%
for i = 1:iter
    % 更新 Ci
    for j = 1:cluster_n
        u_ij_m = U(j,:).^m;
        sum_u_ij = sum(u_ij_m);
        
        c(j,:) = u_ij_m*data./sum_u_ij;
    end
    
    % 更新 uij
    for j = 1:cluster_n
        for k = 1:num_data
            sum1 = 0;
            for j1 = 1:cluster_n
                temp = (norm(data(k,:)-c(j,:))/norm(data(k,:)-c(j1,:))).^(2/(m-1));
                sum1 = sum1 + temp;
            end
            U(j,k) = 1./sum1;
        end
    end
    
    % 計算目標函數 J
    temp1 = zeros(cluster_n,num_data);
    for j = 1:cluster_n
        for k = 1:num_data
            temp1(j,k) = U(j,k)^m*(norm(data(k,:)-c(j,:)))^2;
        end
    end
    J(i) = sum(sum(temp1));    
    
end

% confusion matrix: https://www.mathworks.com/help/stats/confusionmat.html
[~,label] = max(U);
figure
g1 = y';
g2 = label';
C = confusionmat(g1,g2)
confusionchart(C)

%% 繪圖
figure;
plot(J);

figure
hold on
[~,label] = max(U); % 找到所屬的類
gscatter(data(:,1),data(:,2),label)
plot(c(:,1),c(:,2),'kd')

%% function part
% Data set generation from Gaussian classes.
function [X, y] = generate_gauss_classes(m, S, N)
%     disp("In function ");
    [~ ,c]=size(m); % c-dimension
    X = [];
    y = [];
    n=0;
    while(n < N)
        % Generating the [p(j)*N] vectors from each distribution 
        % 1~25 to class1; 26~50 to class2; 51~75 to classes3; 76~100 to classes4
        for j=1:4
            switch j
                case {1,2}
                     t=mvnrnd(m(:,2), S(3:4,:),1);
                     X=[X;t];
                     y=[y ones(1,1)*2];
                case 3
                     t=mvnrnd(m(:,1), S(1:2,:),1);
                     X=[X;t];
                     y=[y ones(1,1)*1];
                case 4
                     t=mvnrnd(m(:,3), S(5:6,:),1);
                     X=[X;t];
                     y=[y ones(1,1)*3];
            end       
        end
        n=n+4;
    end
    X=X';
end

function plot_data(X, y, m)
    [l, N]=size(X);
    [l, c]=size(m);
    pale=["r."; "bo"; "g+"];
    figure();
    hold on
    for i=1:N
        plot(X(1,i), X(2,i), pale(y(i),:));
    end
    for j=1:c
        plot(m(1,j), m(2,j), "k*");
    end
    hold off
end

% Random initialization.
function w=rand_vec(X,m,sed)
    rand('seed',sed)
    mini=min(X');
    maxi=max(X');
    w=rand(size(X,1),m);
    for i=1:m
        w(:,i)=w(:,i).*(maxi'-mini')+mini';
    end
end