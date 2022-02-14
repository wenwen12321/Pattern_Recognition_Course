clc;
clear;

%% (a) generate data set 

m = [1 1; 12 8; 16 1]'; % mean vectors
S = [4 0; 0 4];         % covariance matrices
P = 1/3;
N = 1000;

% for the reproducibility of the results.
randn('seed',0);
% generate data set by gaussian distribution
[X, y]=generate_gauss_classes(m, S, P, N);
% plot data set
plot_data(X,y,m);

%% (b)、(c)
% Bayesian classifier
z_B = bayes_classifier(m, S, P, X);
bayes_error=compute_error(y,z_B, X);
disp(['Bayesian classifier error: ', num2str(bayes_error)]);

% Euclidean classifier
z_E = euclidean_classifier(m,X);
euclidean_error=compute_error(y,z_E, X);
disp(['Euclidean classifier error: ', num2str(euclidean_error)]);

% Mahalanobis classifiers
z_M = mahalanobis_classifier(m,S,X);
mahalanobis_error=compute_error(y,z_M, X);
disp(['Mahalanobis classifier error: ', num2str(mahalanobis_error)]);

%% function part

% Data set generation from Gaussian classes.
function [X, y]=generate_gauss_classes(m, S, P, N)
    [~, c]=size(m);
    X=[];
    y=[];
    for j=1:c
    % Generating the [p(j)*N] vectors from each distribution 
    % 1~333 to class1; 334~666 to class2; 667~1000 to classes3
        if(j==c)
            t = mvnrnd(m(:,j),S,(fix(P*N)+1));
            X = [X ;t];
            y = [y ones(1,fix(P*N)+1)*j];
        else
            t = mvnrnd(m(:,j),S,fix(P*N));
            X = [X ;t];
            y = [y ones(1,fix(P*N))*j];
        end
    end
    X=X';
end

function plot_data(X,y,m)
    [l,N]=size(X); % N=no. of data vectors, l=dimensionality
    [l,c]=size(m); % c=no. of classes
    pale=["r.";"bo";"g+"]; 
    figure(1);
    hold on
    % Plot of the data vectors
    for i=1:N
        plot(X(1,i), X(2,i), pale(y(i),:));
    end
    % Plot of the class means
    for j=1:c 
        plot(m(1,j), m(2, j), 'k*');
    end
    %hold off
end

% Gaussian function evaluation.
function z=comp_gauss_dens_val(m, S, x)
    [l, ~]=size(m);
    z=(1/((2*pi)^(l/2)*det(S)^0.5))*exp(-0.5*(x-m)'*inv(S)*(x-m));
end

% Bayesian classifier ( for Gaussian Processes).
function z = bayes_classifier(m, S, P, X)
    [~, c]=size(m);% c=no. of classes
    [~, N]=size(X);% N=no. of vectors
    for i=1:N
        for j=1:c
            t(j)=P*comp_gauss_dens_val(m(:,j),S,X(:,i));
        end
        % Determining the maximum quantity Pi*p(x|wi)
        %disp(['all t: ', num2str(t)]);
        [num,z(i)]=max(t); % num -> max data, z(i) -> class 
        %disp(['max posterior: ', num2str(z(i))]);
    end
end

% Euclidean distance classifier.
function z = euclidean_classifier(m,X)
    [l, c]=size(m);
    [l, N]=size(X);
    for i=1:N
        for j=1:c
            t(j)=sqrt((X(:,i)-m(:,j))'*(X(:,i)-m(:,j)));
        end
        % Determining the maximum quantity Pi*p(x|wi)
        [num,z(i)]=min(t);
        %disp(['all t: ', num2str(t)]);
        %disp(['min class: ', num2str(z(i))]);
    end
end

% Mahalanobis distance classifier.
function z=mahalanobis_classifier(m,S,X)
    [~, c]=size(m);
    [l, N]=size(X);
    for i=1:N
        for j=1:c
            t(j)=sqrt((X(:,i)-m(:,j))'*inv(S)*(X(:,i)-m(:,j)));
        end
        % Determining the maximum quantity Pi*p(x|wi)
        [num,z(i)]=min(t);
    end
end

 % Classification error evaluation.
function clas_error=compute_error(y,y_est,X)
    [q,N]=size(y);
    %c=max(y);
    clas_error=0;
    disp(['error point: '])
    for i=1:N
        if(y(i)~=y_est(i))
            clas_error=clas_error+1;
            %plot(X(1,i), X(2,i), 'mh');
            disp([num2str(i)]);
        end
    end
    clas_error=clas_error/N;
end