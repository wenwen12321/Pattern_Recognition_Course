clc
clear

%% a. (i) Generate two sets, each one consisting of 100 two-dimensional vectors
m = [2 4; 2.5 10]';
S = [1 0; 0 1];
P = [0.5 0.5]';
N = 200;

randn('seed',0);
[X, y]=generate_gauss_classes(m, S, P, N);
plot_data(X,y,m);
title(' (a) Generate two sets');

% a. (ii) Compute the value of the FDR index for both features.
FDR = [];
for i=1:2
    FDR = [FDR ;FDR_comp(X, y, i)];
end

%% b. (i)
m = [2 4; 2.5 10]';
S = [0.25 0; 0 0.25];
P = [0.5 0.5]';
N = 200;

randn('seed',0);
[X, y]=generate_gauss_classes(m, S, P, N);
plot_data(X,y,m);
title(' (b) Generate two sets');

% b. (ii) Compute the value of the FDR index for both features.
FDR = [];
for i=1:2
    FDR = [FDR; FDR_comp(X, y, i)];
end

%% c. Discuss the result (write in report)

%% function part
% Data set generation from Gaussian classes.
function [X, y]=generate_gauss_classes(m, S, P, N)
    [~, c]=size(m); % c-dimension
    X=[];
    y=[];
    for j=1:c
    % Generating the [p(j)*N] vectors from each distribution 
    % 1~25 to class1; 26~50 to class2; 51~75 to classes3; 76~100 to classes4
        t = mvnrnd(m(:,j),S,fix(P(j)*N));
        X = [X ;t];
        y = [y ones(1,fix(P(j)*N))*j];
    end
    X=X';
end


function plot_data(X,y,m)
    [l,N]=size(X);
    [l,c]=size(m);
    pale=["r.";"bo";"g+"; "y*"];
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

function FDR=FDR_comp(X,y,ind)
    [l,N]=size(X);
    c=max(y);
    for i=1:c
        y_temp=(y==i);
        X_temp=X(ind,y_temp);
        m(i)=mean(X_temp);
        vari(i)=var(X_temp);
    end
    a=nchoosek(1:c,2);
    q=(m(a(:,1))-m(a(:,2))).^ 2 ./ (vari(a(:,1))+vari(a(:,2)))';
    FDR=sum(q);
end