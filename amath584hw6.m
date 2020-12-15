clc;clear;close all;
%% load data
[X,n,Height,Width] = loadMNISTImages('train-images.idx3-ubyte');
L = loadMNISTLabels('train-labels.idx1-ubyte');
L(L==0) = 10;
SX = sparse(X);
SB = sparse(L',(1:n)',ones(n,1),10,n);
B = full(SB);
clearvars SB SX
%% /
A1 = B/X; modename1 = 'B/X';
%% pinv
A2 = B*pinv(X); modename2 = 'B*pinv';
%% lasso 1
A3 = zeros(10,Width*Height);
for i = 1:10
    b = (B(i,:))';
    a = lasso(X',b,'Lambda',0.1,'alpha',0.2);
    A3(i,:) = a';
end
modename3 = 'lasso (\lambda = 0.1 \alpha = 0.2)';
%% lasso 2
A4 = zeros(10,Width*Height);
for i = 1:10
    b = (B(i,:))';
    a = lasso(X',b,'Lambda',0.1,'alpha',0.5);
    A4(i,:) = a';
end
modename4 = 'lasso (\lambda = 0.1 \alpha = 0.5)';
%%
SolAnalyze(A1,modename1,0,X,L,Width,Height);
%%
SolAnalyze(A2,modename2,0,X,L,Width,Height);
%%
SolAnalyze(A3,modename3,0.01,X,L,Width,Height);
%%
SolAnalyze(A4,modename4,0.01,X,L,Width,Height);
%%
SA = A1;
SA(abs(A1) > 0.2) = 0;
SolAnalyze(SA,modename1,0,X,L,Width,Height);
%%

