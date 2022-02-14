clc
clear

% from scipy.stats import multivariate_normal

%% a
m = [1 1; 5 5; 9 1]';
S = [1 0.4; 0.4 1; 
    1 -0.6; -0.6 1; 
    1 0; 0 1];
N = 500;

% initialization P0, M0, S0, e, iter
P0 = [0.25, 0.5, 0.25]';
M0 = [0 0]';
S0 = [0 0];
% e = 0;
% k =  0;% k個數
iter = 20;

randn('seed', 0);
rng('default') % For reproducibility
% generate data
[X, y] = generate_gauss_classes(m,S,N);
plot_data(X,y,m);
title('a. Generated data');
legend('Gaussian 2', '', 'Gaussian 1', 'Gaussian 3');

% Posterior probability of Gaussian mixture component
% reference: https://www.mathworks.com/help/stats/gmdistribution.posterior.html#d123e650844
gm = fitgmdist(X',3);

% Visualize the fitting model gm by using pdf and fcontour
figure
scatter(X(1,:),X(2,:),10,'.') % Scatter plot with points of size 10
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
fcontour(gmPDF,[-1.5 11 -1.5 9])
c1 = colorbar;
ylabel(c1,'Probability Density Function')

% Compute the posterior probabilities of the componments.
% Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
[posterior_p, nlogL] = posterior(gm, X');

% Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
figure
scatter(X(1,:),X(2,:),10,posterior_p(:,1))
c2 = colorbar;
ylabel(c2,'Posterior Probability of Component 1')

% Plot the posterior probabilities of Component 2.
figure
scatter(X(1,:),X(2,:),10,posterior_p(:,3))
c3 = colorbar;
ylabel(c3,'Posterior Probability of Component 2')

% Plot the posterior probabilities of Component 3.
figure
scatter(X(1,:),X(2,:),10,posterior_p(:,2))
c4 = colorbar;
ylabel(c4,'Posterior Probability of Component 3')

% disp("P = ")
disp("P1 = "+gm.ComponentProportion(:,1))
disp("P2 = "+gm.ComponentProportion(:,3))
disp("P3 = "+gm.ComponentProportion(:,2))

disp(" ")
disp("mu1 = ")
disp(gm.mu(1,:))
disp("mu2 = ")
disp(gm.mu(3,:))
disp("mu3 = ")
disp(gm.mu(2,:))

disp(" ")
disp("sigma1 = ")
disp(gm.Sigma(:,:,1))
disp("sigma2 = ")
disp(gm.Sigma(:,:,3))
disp("sigma3 = ")
disp(gm.Sigma(:,:,2))

%% b
m = [1 1; 3.5 3.5; 6 1]';

% randn('seed', 0);
% rng('default') % For reproducibility
% % generate data
% [X, y] = generate_gauss_classes(m,S,N);
% plot_data(X,y,m);
% title('b. Generated data');
% legend('Gaussian 2', '', 'Gaussian 1', 'Gaussian 3');
% 
% % Posterior probability of Gaussian mixture component
% % reference: https://www.mathworks.com/help/stats/gmdistribution.posterior.html#d123e650844
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
% % Compute the posterior probabilities of the componments.
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% [posterior_p, nlogL] = posterior(gm, X');
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,1))
% c2 = colorbar;
% ylabel(c2,'Posterior Probability of Component 1')
% 
% % Plot the posterior probabilities of Component 2.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,2))
% c3 = colorbar;
% ylabel(c3,'Posterior Probability of Component 2')
% 
% % Plot the posterior probabilities of Component 3.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,3))
% c4 = colorbar;
% ylabel(c4,'Posterior Probability of Component 3')
% 
% % disp("P = ")
% disp("P1 = "+gm.ComponentProportion(:,1))
% disp("P2 = "+gm.ComponentProportion(:,2))
% disp("P3 = "+gm.ComponentProportion(:,3))
% 
% disp(" ")
% disp("mu1 = ")
% disp(gm.mu(1,:))
% disp("mu2 = ")
% disp(gm.mu(2,:))
% disp("mu3 = ")
% disp(gm.mu(3,:))
% 
% disp(" ")
% disp("sigma1 = ")
% disp(gm.Sigma(:,:,1))
% disp("sigma2 = ")
% disp(gm.Sigma(:,:,2))
% disp("sigma3 = ")
% disp(gm.Sigma(:,:,3))

%% c
% m = [1 1; 2 2; 3 1]';
% 
% randn('seed', 0);
% rng('default') % For reproducibility
% % generate data
% [X, y] = generate_gauss_classes(m,S,N);
% plot_data(X,y,m);
% title('c. Generated data');
% legend('Gaussian 2', '', 'Gaussian 1', 'Gaussian 3');
% 
% % Posterior probability of Gaussian mixture component
% % reference: https://www.mathworks.com/help/stats/gmdistribution.posterior.html#d123e650844
% gm = fitgmdist(X',3);
% 
% % Visualize the fitting model gm by using pdf and fcontour
% figure
% scatter(X(1,:),X(2,:),10,'.') % Scatter plot with points of size 10
% hold on
% gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
% fcontour(gmPDF,[-1.5 7 -1.5 6])
% c1 = colorbar;
% ylabel(c1,'Probability Density Function')
% 
% % Compute the posterior probabilities of the componments.
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% [posterior_p, nlogL] = posterior(gm, X');
% 
% % Plot the posterior probabilities of Component 1 by using the scatter function. Use the circle colors to visualize the posterior probability values.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,3))
% c2 = colorbar;
% ylabel(c2,'Posterior Probability of Component 1')
% 
% % Plot the posterior probabilities of Component 2.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,1))
% c3 = colorbar;
% ylabel(c3,'Posterior Probability of Component 2')
% 
% % Plot the posterior probabilities of Component 3.
% figure
% scatter(X(1,:),X(2,:),10,posterior_p(:,2))
% c4 = colorbar;
% ylabel(c4,'Posterior Probability of Component 3')
% 
% % disp("P = ")
% disp("P1 = "+gm.ComponentProportion(:,3))
% disp("P2 = "+gm.ComponentProportion(:,1))
% disp("P3 = "+gm.ComponentProportion(:,2))
% 
% disp(" ")
% disp("mu1 = ")
% disp(gm.mu(3,:))
% disp("mu2 = ")
% disp(gm.mu(1,:))
% disp("mu3 = ")
% disp(gm.mu(2,:))
% 
% disp(" ")
% disp("sigma1 = ")
% disp(gm.Sigma(:,:,3))
% disp("sigma2 = ")
% disp(gm.Sigma(:,:,1))
% disp("sigma3 = ")
% disp(gm.Sigma(:,:,2))

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

% function [X, y]=generate_gauss_classes2(m, S, P, N)
%     [~, c]=size(m); % c-dimension
%     X=[];
%     y=[];
%     for j=1:c
%     % Generating the [p(j)*N] vectors from each distribution 
%     % 1~125 to class1; 126~375 to class2; 376~500 to classes3;
%         t = mvnrnd(m(:,j),S(j*2-1:j*2,:),fix(P(j)*N));
%         X = [X ;t];
%         y = [y ones(1,fix(P(j)*N))*j];
%     end
%     X=X';
% end

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