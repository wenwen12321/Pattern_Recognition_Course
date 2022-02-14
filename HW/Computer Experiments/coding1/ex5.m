clc;
clear;

%% (a)

m = [1 1; 4 4; 8 1]';
S = [2 0; 0 2];
P = [1/3, 1/3, 1/3]';
P2= [0.8 0.1 0.1]';
N = 1000;

%for the reproducibility of the results.
randn('seed',0);
[X5, y]=generate_gauss_classes(m, S, P, N);
[X5_, y_]=generate_gauss_classes(m, S, P2, N);
%plot_data(X5,y,m);
%plot_data(X5_,y_,m);

%% (b)、(c)
% Bayesian classifier
z = bayes_classifier(m, S, P, X5);
z_ = bayes_classifier(m, S, P2, X5_);
bayes_error=compute_error(y,z);
bayes_error_=compute_error(y_,z_);
disp(['X5 Bayesian classifier error: ', num2str(bayes_error)]);
disp(['X5^(T) Bayesian classifier error: ', num2str(bayes_error_)]);
plot_result(X5, z);
plot_result(X5_, z_);

% Euclidean classifier
z = euclidean_classifier(m,X5);
z_ = euclidean_classifier(m,X5_);
euclidean_error=compute_error(y,z);
euclidean_error_=compute_error(y_,z_);
fprintf('\n');
disp(['X5 Euclidean classifier error: ', num2str(euclidean_error)]);
disp(['X5^(T) Euclidean classifier error: ', num2str(euclidean_error_)]);
plot_result(X5, z);
plot_result(X5_, z_);

%% function part
function [X, y]=generate_gauss_classes(m, S, P, N)
    [~, c]=size(m);
    X=[];
    y=[];
    for j=1:c
        % Generating the [p(j)*N] vectors from each distribution
        if(j==c && P(j)==1/3)
            t = mvnrnd(m(:,j),S,(fix(P(j)*N)+1));
            X = [X ;t];
            y = [y ones(1,fix(P(j)*N)+1)*j];
        else
            t = mvnrnd(m(:,j),S,fix(P(j)*N));
            X = [X ;t];
            y = [y ones(1,fix(P(j)*N))*j];
        end
    end
    X=X';
end

function plot_data(X,y,m)
    [l,N]=size(X);
    [l,c]=size(m);
    pale=["r.";"bo";"g+"];
    figure();
    hold on
    for i=1:N
        plot(X(1,i), X(2,i), pale(y(i),:));
    end
    for j=1:c 
        plot(m(1,j), m(2, j), 'k*');
    end
    hold off
end

function plot_result(X,y)
    [l,N]=size(X);
    pale=["r.";"bo";"g+"];
    figure();
    hold on
    for i=1:N
        plot(X(1,i), X(2,i), pale(y(i),:));
    end
    hold off
end

function z=comp_gauss_dens_val(m, S, x)
    [l, ~]=size(m);
    z=(1/((2*pi)^(l/2)*det(S)^0.5))*exp(-0.5*(x-m)'*inv(S)*(x-m));
end

function z = bayes_classifier(m, S, P, X)
    [~, c]=size(m);% c=no. of classes
    [~, N]=size(X);% N=no. of vectors
    for i=1:N
        for j=1:c
            t(j)=P(j)*comp_gauss_dens_val(m(:,j),S,X(:,i));
        end
        % Determining the maximum quantity Pi*p(x|wi)
        %disp(['all t: ', num2str(t)]);
        [num,z(i)]=max(t); % num -> max data, z(i) -> class 
        %disp(['max posterior: ', num2str(z(i))]);
    end
end

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

function clas_error=compute_error(y,y_est)
    [q,N]=size(y);
    %c=max(y);
    clas_error=0;
    for i=1:N
        if(y(i)~=y_est(i))
            clas_error=clas_error+1;
        end
    end
    clas_error=clas_error/N;
end