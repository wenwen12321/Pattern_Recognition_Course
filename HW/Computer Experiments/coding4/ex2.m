clc
clear

%%
m = [1 1; 5 5; 9 1]';
S = [1 0.4; 0.4 1; 
    1 -0.6; -0.6 1; 
    1 0; 0 1];
N = 500;

% % initialization P0, M0, S0, e
% P0 = [0.25, 0.5, 0.25]';
% M0 = [0 0]';
% S0 = [0 0];
% e = 0;

randn('seed', 0);
rng('default') % For reproducibility
% generate data
[X, y] = generate_gauss_classes(m,S,N);
plot_data(X,y,m);
title('generated data');
legend('Gaussian 2', '', 'Gaussian 1', 'Gaussian 3');

%% k=2
% k=2;
% sed=0;
% w = rand_vec(X, k, sed);
% [w, bel]= k_means(X, w);
% 
% % clustering data
% plot_data(X,bel,w);
% title('k = 2,  clustering data');
% legend('cluster 2', '', 'cluster 1');
% 
% disp("mu1 = ")
% disp(w(:,1))
% disp("mu2 = ")
% disp(w(:,2))
% % disp("mu3 = ")
% % disp(w(:,3))
% 
% gm = fitgmdist(X',2);
% 
% % Visualize the fitting model gm by using pdf and fcontour
% figure
% scatter(X(1,:),X(2,:),10,'.') % Scatter plot with points of size 10
% hold on
% gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
% fcontour(gmPDF,[-1.5 11 -1.5 9])
% c1 = colorbar;
% ylabel(c1,'Probability Density Function')
% 
% disp("P1 = "+gm.ComponentProportion(:,2))
% disp("P2 = "+gm.ComponentProportion(:,1))
% 
% disp(" ")
% disp("sigma1 = ")
% disp(gm.Sigma(:,:,2))
% disp("sigma2 = ")
% disp(gm.Sigma(:,:,1))
% 
% [posterior_p, nlogL] = posterior(gm, X');
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,2))
% c2 = colorbar;
% ylabel(c2,'Posterior Probability of Component 1')
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,1))
% c3 = colorbar;
% ylabel(c3,'Posterior Probability of Component 2')

%% k=3
% k=3;
% sed=0;
% w = rand_vec(X, k, sed);
% [w, bel]= k_means(X, w);
% 
% % clustering data
% plot_data(X,bel,w);
% title('k = 3,  clustering data');
% legend('cluster 1', '', 'cluster 2', 'cluster 3');
% 
% disp("mu1 = ")
% disp(w(:,1))
% disp("mu2 = ")
% disp(w(:,2))
% disp("mu3 = ")
% disp(w(:,3))
% 
% gm = fitgmdist(X',3);
% 
% % Visualize the fitting model gm by using pdf and fcontour
% figure
% scatter(X(1,:),X(2,:),10,'.') % Scatter plot with points of size 10
% hold on
% gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
% fcontour(gmPDF,[-1.5 11 -1.5 9])
% c1 = colorbar;
% ylabel(c1,'Probability Density Function')
% 
% disp("P1 = "+gm.ComponentProportion(:,1))
% disp("P2 = "+gm.ComponentProportion(:,3))
% disp("P3 = "+gm.ComponentProportion(:,2))
% 
% disp(" ")
% disp("sigma1 = ")
% disp(gm.Sigma(:,:,1))
% disp("sigma2 = ")
% disp(gm.Sigma(:,:,3))
% disp("sigma3 = ")
% disp(gm.Sigma(:,:,2))
% 
% [posterior_p, nlogL] = posterior(gm, X');
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,1))
% c2 = colorbar;
% ylabel(c2,'Posterior Probability of Component 1')
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,3))
% c3 = colorbar;
% ylabel(c3,'Posterior Probability of Component 2')
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,2))
% c4 = colorbar;
% ylabel(c4,'Posterior Probability of Component 3')

%% k=4
k=4;
sed=0;
w = rand_vec(X, k, sed);
[w, bel]= k_means(X, w);

% clustering data
plot_data(X,bel,w);
title('k = 4,  clustering data');
legend('cluster 2', '', 'cluster 1', 'cluster 4', 'cluster 3');

disp("mu1 = ")
disp(w(:,1))
disp("mu2 = ")
disp(w(:,2))
disp("mu3 = ")
disp(w(:,3))
disp("mu4 = ")
disp(w(:,4))

gm = fitgmdist(X',4);

% Visualize the fitting model gm by using pdf and fcontour
figure
scatter(X(1,:),X(2,:),10,'.') % Scatter plot with points of size 10
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
fcontour(gmPDF,[-1.5 11 -1.5 9])
c1 = colorbar;
ylabel(c1,'Probability Density Function')

disp("P1 = "+gm.ComponentProportion(:,1))
disp("P2 = "+gm.ComponentProportion(:,2))
disp("P3 = "+gm.ComponentProportion(:,3))
disp("P4 = "+gm.ComponentProportion(:,4))

disp(" ")
disp("sigma1 = ")
disp(gm.Sigma(:,:,1))
disp("sigma2 = ")
disp(gm.Sigma(:,:,2))
disp("sigma3 = ")
disp(gm.Sigma(:,:,3))
disp("sigma4 = ")
disp(gm.Sigma(:,:,4))
% 
% [posterior_p, nlogL] = posterior(gm, X');
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,1))
% c2 = colorbar;
% ylabel(c2,'Posterior Probability of Component 1')
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,2))
% c3 = colorbar;
% ylabel(c3,'Posterior Probability of Component 2')
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,3))
% c4 = colorbar;
% ylabel(c4,'Posterior Probability of Component 3')
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,4))
% c5 = colorbar;
% ylabel(c5,'Posterior Probability of Component 4')

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
    pale=["r."; "bo"; "g+"; "cv"];
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
    