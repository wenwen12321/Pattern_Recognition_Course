clc
clear

%% a, (i) Generate four sets, each one consisting of 100 two-dimensional vectors.
m = [-10 -10; -10 10; 10 -10; 10 10]';
S = [0.2 0; 0 0.2];
P = [0.25 0.25 0.25 0.25]';
N = 400;

randn('seed',0);
[X, y]=generate_gauss_classes(m, S, P, N);
plot_data(X,y,m);
title(' (a) Generate four sets');

% a. (ii) Compute the Sw, Sb, and Sm scatter matrices.
[Sw, Sb, Sm]=scatter_mat(X,y); 
disp("Sw = ");
disp(Sw);
disp("Sb = ");
disp(Sb);
disp("Sm = ");
disp(Sm);

% %a. (iii) Compute the value for the criterion J 3.
J3 = trace(inv(Sw)*Sm);
disp("J3 = " + J3);

%% b. (i) Repeat(a) when mean = [1, 1]T , [1, 1]T , [1, 1]T , [1, 1]T .
m = [-1 -1; -1 1; 1 -1; 1 1]';
S = [0.2 0; 0 0.2];
P = [0.25 0.25 0.25 0.25]';
N = 400;

randn('seed',0);
[X, y]=generate_gauss_classes(m, S, P, N);
plot_data(X,y,m);
title('(b) Generate four sets');

% b. (ii) Compute the Sw, Sb, and Sm scatter matrices.
[Sw, Sb, Sm]=scatter_mat(X,y); 
disp("Sw = ");
disp(Sw);
disp("Sb = ");
disp(Sb);
disp("Sm = ");
disp(Sm);

% b. (iii) Compute the value for the criterion J 3.
J3 = trace(inv(Sw)*Sm);
disp("J3 = " + J3);

%% c. (i)
m = [-10 -10; -10 10; 10 -10; 10 10]';
S = [3 0; 0 3];
P = [0.25 0.25 0.25 0.25]';
N = 400;

randn('seed',0);
[X, y]=generate_gauss_classes(m, S, P, N);
plot_data(X,y,m);
title('(c) Generate four sets');

% c. (ii) Compute the Sw, Sb, and Sm scatter matrices.
[Sw, Sb, Sm]=scatter_mat(X,y); 
disp("Sw = ");
disp(Sw);
disp("Sb = ");
disp(Sb);
disp("Sm = ");
disp(Sm);

% c. (iii) Compute the value for the criterion J 3.
J3 = trace(inv(Sw)*Sm);
disp("J3 = " + J3);

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

function [Sw, Sb, Sm] = scatter_mat(X, y)
    [l, N]=size(X);
    c=max(y);
    %Computation of class mean vectors, a priori prob, and Sw
    m=[];
    Sw=zeros(l);
    for i=1:c
        y_temp=(y==i);
        X_temp=X(:, y_temp);
        P(i)=sum(y_temp) / N;
        m(:, i)=(mean(X_temp'))';
        Sw=Sw+P(i)*cov(X_temp');
    end
    %Computation of Sb
    m0=(sum(((ones(l,1)*P).*m)'))';
    Sb=zeros(l);
    for i=1:c
        Sb=Sb+P(i)*((m(:,i)-m0)*(m(:,i)-m0)');
    end
    %Computation of Sm
    Sm=Sw+Sb;
end

function J3=J3_comp(Sw, Sm)
    J3=trace(inv(Sw)*Sm)
end