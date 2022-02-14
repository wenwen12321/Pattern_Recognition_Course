clc
clear

%% a.
w = [1 1]';
w0 = 0;
a =10;
e = 1;
N = 1000;
sed = 0;

X = generate_hyper(w, w0, a, e, N, sed);

%% b. PCA analysis
[pc,variances]=pcacov(cov(X'))


%% function part
function X=generate_hyper(w,w0,a,e,N,sed)
    l=length(w);
    t=(rand(l-1,N)-.5)*2*a;
    t_last=-(w(1:l-1)/w(l))'*t + 2*e*(rand(1,N)-.5)-(w0/w(l));
    X=[t; t_last];
    %Plots for the 2d and 3d case
    if(l==2)
        figure(1), plot(X(1,:),X(2,:),'.b')
    elseif(l==3)
        figure(1), plot3(X(1,:),X(2,:),X(3,:),'.b')
    end
    figure(1), axis equal
end