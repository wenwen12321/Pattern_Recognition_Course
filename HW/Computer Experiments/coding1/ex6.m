clc;
clear;

%% (a)

m = [1 1; 8 6; 13 1]';
S = [6 0; 0 6];
P = 1/3;
N = 1000;
% k = 1 and 11

%for the reproducibility of the results.
randn('seed',0);
[X, y]=generate_gauss_classes(m, S, P, N);
[Z, y1]=generate_gauss_classes(m, S, P, N);

k=1;
z1 = k_nn_classifier(Z,y1,k,X);

k=11;
z11 = k_nn_classifier(Z,y1,k,X);

plot_data(X,z1,m);
plot_data(X, z11,m);

%% function part
function [X, y]=generate_gauss_classes(m, S, P, N)
    [~, c]=size(m);
    X=[];
    y=[];
    for j=1:c
    % Generating the [p(j)*N] vectors from each distribution
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

function z = k_nn_classifier(Z,v,k,X)
    [l,N1]=size(Z);
    [l,N]=size(X);
    c=max(v) % the number of classes
    
    % Computaion of the (squared) Eculidean distance
    % of a point from each reference vector
    for i=1:N
        dist=sum((X(:,i)*ones(1,N1)-Z).^2);
        % Sorting the above distances in ascending order
        [sorted,nearest]=sort(dist);    
        refe=zeros(1,c)
        for q=1:k
            class=v(nearest(q));
            refe(class)=refe(class)+1;
        end
        [val,z(i)]=max(refe);
        %disp(['z: ', num2str(val)]);
    end
end