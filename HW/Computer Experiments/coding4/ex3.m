clc
clear

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

k=3;
sed=0;
w = rand_vec(X, k, sed);

%% q=2
q=2;
%% q=3
% q=3;

% reference: https://www.mathworks.com/help/fuzzy/fcm.html#buv_wo_-4
[centers,U] = fcm(X',3,q);

% Classify each data point into the cluster with the largest membership value.
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);

% Plot the clustered data and cluster centers.
plot(X(1,index1),X(2,index1),'ob');
hold on;
plot(X(1,index2),X(2,index2),'or')
plot(X(1,index3),X(2,index3),'og')
plot(centers(1,1),centers(1,2),'xk','MarkerSize',15,'LineWidth',3)
plot(centers(2,1),centers(2,2),'xk','MarkerSize',15,'LineWidth',3)
plot(centers(3,1),centers(3,2),'xk','MarkerSize',15,'LineWidth',3)
hold off

% confusion matrix: https://www.mathworks.com/help/stats/confusionmat.html
figure
% g1 = y';
[~,label] = max(U);
% g2 = label';
C = confusionmat(y',label');
confusionchart(C);

disp("mu1 = ")
disp(centers(1,:))
disp("mu2 = ")
disp(centers(2,:))
disp("mu3 = ")
disp(centers(3,:))

gm = fitgmdist(X',3);

% Visualize the fitting model gm by using pdf and fcontour
figure
scatter(X(1,:),X(2,:),10,'.') % Scatter plot with points of size 10
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
fcontour(gmPDF,[-1.5 11 -1.5 9])
c1 = colorbar;
ylabel(c1,'Probability Density Function')

[~,index1_size]=size(index1);
[~,index2_size]=size(index2);
[~,index3_size]=size(index3);

index1_size = index1_size/500;
index2_size = index2_size/500;
index3_size = index3_size/500;

disp("P1 = "+index1_size)
disp("P2 = "+index2_size)
disp("P3 = "+index3_size)

disp(" ")
disp("sigma1 = ")
disp(gm.Sigma(:,:,1))
disp("sigma2 = ")
disp(gm.Sigma(:,:,3))
disp("sigma3 = ")
disp(gm.Sigma(:,:,2))

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

% k_means algorithm
function [w, bel]=k_means(X,w)
    [l,N]=size(X);
    [l,m]=size(w);
    e=1;
    iter=0;
    while(e~=0)
        iter=iter+1;
        w_old=w;
        dist_all=[];
        for j=1:m
            dist=sum(((ones(N,1)*w(:,j)'-X').^2)');
            dist_all=[dist_all;dist];
        end
        [q1, bel]=min(dist_all);
        for j=1:m
            if(sum(bel==j)~=0)
                w(:,j)=sum(X'.*((bel==j)'*ones(1,l)))/sum(bel==j);
            end
        end
        e=sum(abs(w-w_old));
    end
end
    