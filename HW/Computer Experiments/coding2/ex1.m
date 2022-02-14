clc;
clear;

%% (a)
Pa=0.25;
N=1000;
% generates the Bernoulli samples 𝑋
Xa = rand(N, 1) < Pa;

disp('(a) estimate p is:');
P_a= estimate_p(Xa, N);

%% (b)
Pb=0.5;
N=1000;
% generates the Bernoulli samples 𝑋
Xb = randi([0, 1], N, 1);

disp('(b) estimate p is:');
P_b = estimate_p(Xb, N);

%% function part
function p = estimate_p(X, N)
    p = 0;    
    for i = 1:N    
        p = p + X(i);
    end
    p = p / N;
    disp(p);
end