clc
clear;
randn('seed',0);

%% (a)
m = [1, 1]';
S = [5 3; 3 4];
N = 1000;
% Data set generation from Gaussian classes.
Xa = generate_gauss_classes(m, S, 1, N);

disp('(a) estimate 𝛍 is:')
Mua_ML = estimate_mu(Xa, N);

disp('(a) estimate 𝚺 is:');
Sigmaa_ML = estimate_sigma(Xa, Mua_ML, N);


%% (b)
m = [10, 5]';
S = [7 4; 4 5];
N =1000;
% Data set generation from Gaussian classes.
Xb = generate_gauss_classes(m, S, 1, N);

disp('(b) estimate 𝛍 is:');
Mub_ML = estimate_mu(Xb, N);

disp('(b) estimate 𝚺 is:');
Sigmab_ML = estimate_sigma(Xb, Mub_ML, N);

%% function part
% Data set generation from Gaussian classes.
function X=generate_gauss_classes(m, S, P, N)
    X=[];
    t = mvnrnd(m(:,1), S, N);
    X = [X ;t];
    figure();
    plot(X(:,1), X(:,2), '+');
    %X=X';
end

% estimate mu => sample mean
% ( 𝚺 Xi) / N
function mu_ML = estimate_mu(X,N)
    mu_ML = [];
    tempX = 0;
    tempY = 0;
    
    for i=1:N
        tempX = tempX + X(i, 1);
        tempY = tempY + X(i, 2);
    end
    
    tempX = tempX/N;
    tempY = tempY/N;
    
    mu_ML = [tempX; tempY];
    disp(mu_ML);
end

% ( 𝚺 (Xi-𝛍)(Yi-𝛍) ) / N
function Sigma_ML =  estimate_sigma(X, Mu_ML, N)
    Sigma_ML = [];
    sumx = 0;
    sumy = 0;
    sumxy = 0;
    for i = 1:N
        sumx = sumx + ((X(i, 1) - Mu_ML(1, 1)) * (X(i, 1) - Mu_ML(1, 1))');
        sumy = sumy + ((X(i, 2) - Mu_ML(2, 1)) * (X(i, 2) - Mu_ML(2, 1))');
        sumxy = sumxy + ((X(i, 1) - Mu_ML(1, 1)) * (X(i, 2) - Mu_ML(2, 1))');
    end    
    sumx = sumx / N;
    sumy = sumy / N;
    sumxy = sumxy / N;

    Sigma_ML = [sumx sumxy; sumxy sumy];
    disp(Sigma_ML);
end