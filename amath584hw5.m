clc; clear; close all;
%% I(a)
%	Generate a random, symmetric matrix A which is m by m where m = 10. 
%	Use the EIGS command in MATLAB (or the equivalent in Python) 
%	to give you the ground truth eigenvalues and eigenvectors.
m = 10;
A = rand(m,m);
A = (A+A')./2;
[V,D] = eigs(A,10);
EigenV = sort(diag(D),'descend'); 
figure(1);plot(1:10,EigenV);ylabel('eigenvalue');
%% (b) Find the largest eigenvalue with the power iteration method. 
%   Compare the accuracy of the method as a function of iterations.
max_iter = 20;
err = zeros(1,max_iter); %error (inversely proportional to accuracy)
v0 = rand(m,1); v0 = v0/norm(v0);
for i = 1 : max_iter
   [lambda,v] = PowerIter(A,v0,i);
    err(i) = abs(lambda - D(1,1));
end
figure(2);
subplot(1,2,1);plot(1:max_iter,err);xlabel('iteration') ;ylabel('error of PowerIter');
subplot(1,2,2);semilogy(1:max_iter,err);xlabel('iteration') ;ylabel('error (log) of PowerIter');

disp('largest eigenvalue:');
disp(lambda);
%% (c) Find all ten eigenvalues by Rayleigh Quotient iteration and guessing initial "eigenvectors". 
%Compare the accuracy of the method as a function of iterations and discuss your initial guesses 
%to and all eigen-value/eigenvector pairs.

%   Compare the accuracy of RQIter as a function of iterations.
max_iter = 20;
err = zeros(1,max_iter); %accuracy
v0 = rand(m,1);
v0 = v0/norm(v0);
for i = 1 : max_iter
   [lambda,v] = RQIter(A,v0,i,0);
    err(i) = abs(lambda - D(1,1));
end
figure(3);
subplot(1,2,1);plot(1:max_iter,err);xlabel('iteration') ;ylabel('error of RQIter');
subplot(1,2,2);semilogy(1:max_iter,err);xlabel('iteration') ;ylabel('error (log) of RQIter');
%%
% Try to get all the eigenvalues
GuessNum = 10;
EigenV1 = zeros(1,GuessNum);
EigenVec1 = zeros(m, GuessNum);
v0 = rand(m,1);
for i= 1:GuessNum
    [EigenV1(i),EigenVec1(:,i)] = RQIter(A,v0,100,1E-15);
    v0 = rand(m,1);v0 = v0/norm(v0);
    v0 = v0 - EigenVec1(:,1:i)*((EigenVec1(:,1:i))'*v0);
    % the next initial guess and the previous eigenvectors are orthogonal.
end
EigenV1 = sort(EigenV1,'descend');

figure(4);plot(1:10,EigenV,'o');ylabel('eigenvalue');
hold on 
plot(1:10,EigenV1,'*'); hold off;
legend('ground truth','RQIter');
disp('norm of EigenV - EigenV1:');
disp(norm(EigenV' - EigenV1));
%%
clc;clear;
%% (d) Repeat with a random matrix that is not symmetric. 
% Be sure to plot the eigenvalue in the complex plane.
m = 10;
A = rand(m,m);
[V,D] = eigs(A,10);
EigenV = sort(diag(D),'descend');
%%
figure(5);plot(EigenV,'*');xlabel('Re(\lambda)');ylabel('Im(\lambda)');
%% Find the largest magnitude eigenvalue with the power iteration method.
max_iter = 30;
err = zeros(1,max_iter); %accuracy
v0 = rand(m,1); v0 = v0/norm(v0);
for i = 1 : max_iter
   [lambda,v] = PowerIter(A,v0,i);
    err(i) = abs(lambda - D(1,1));
end
figure(6);
subplot(1,2,1);plot(1:max_iter,err);xlabel('iteration') ;ylabel('error of PowerIter');
subplot(1,2,2);semilogy(1:max_iter,err);xlabel('iteration') ;ylabel('error (log) of PowerIter');

disp('largest magnitude eigenvalue:');
disp(lambda);
%%
clearvars i;
warning('off');
% Try to get all the eigenvalues
GuessNum = 1000;
EigenV3 = zeros(1,GuessNum);
EigenVec3 = zeros(m, GuessNum);
v0 = 2*rand(m,1)-ones(m,1)+ (2*rand(m,1)-ones(m,1))*1i;
for j= 1:GuessNum
    [EigenV3(j),EigenVec3(:,j)] = RQIter(A,v0,15,1E-12);
    v0 = 2*rand(m,1)-ones(m,1)+ (2*rand(m,1)-ones(m,1))*1i;
    v0 = v0/norm(v0);
end
EigenV3 = sort(EigenV3,'descend');
figure(7);plot(EigenV,'o');ylabel('eigenvalue');
hold on
plot(EigenV3,'*'); hold off;
legend('groud truth', 'results of RQIter');
%%
clc;clearvars;
%% II import data
CDirNum = 37; % the images in the last folder are used as test images
CDirIndex = [1:13 15:39]; % I didn't notice in the first place that yaleB14 does not exist
CImgNum = 64; %number of images in each folder
CImgTest = double(imread('CroppedYale\YaleB39\yaleB39_P00A+000E+00.pgm')); %a sample image
[CImgHeight, CImgWidth]= size(CImgTest);
CRawImg = zeros(CImgHeight, CImgWidth, CImgNum*CDirNum);
Ccount = 0; % counting the number of images
CImgAvg = zeros(CImgHeight, CImgWidth); % store average face
for i = 1:CDirNum 
    ImgDir = ['CroppedYale\yaleB',num2str(CDirIndex(i),'%02d')];
    ImgDirName = dir(ImgDir);
    for j = 1: length( ImgDirName )
        if(isequal(ImgDirName(j).name, '.')||...
           isequal(ImgDirName(j).name, '..'))
           continue;
        else
            Ccount=Ccount+1;
            CRawImg(:,:,Ccount) = imread([ImgDir,'\',ImgDirName(j).name]);
            CImgAvg = CImgAvg + CRawImg(:,:,Ccount);
        end  
    end
end
CImgAvg = CImgAvg./Ccount; %caculate the average image
CImgAvg = reshape(CImgAvg, [CImgHeight*CImgWidth, 1])*ones(1,CImgNum*CDirNum);
CImgData = squeeze(reshape(CRawImg,[CImgHeight*CImgWidth,1,CImgNum*CDirNum]));
%% compute the square correlation matrix of cropped matrix
CImgData = CImgData - CImgAvg;
clearvars -except CImgData Ccount CImgHeight CImgWidth
C = CImgData' * CImgData; %square correlation matrix
%% (a) Power iterate on the matrix of images 
% and get the dominant eigenvector and eigenvalue. 
% Compare it to the leading order SVD mode.
[CU,CS,CV] = svd(CImgData, 'econ'); %The svd decomp of data
u1 = CU(:,1);
M1 = reshape(u1,CImgHeight,CImgWidth); % The first mode using SVD
%%
[lambda,v] = PowerIter(C,rand(Ccount,1),100);
u2 = -abs(CImgData*v./sqrt(lambda)); % -abs is for make u1 and u2 into the same direction
M2 = reshape(u2, CImgHeight,CImgWidth); % The first mode using power iteration

figure(8);
subplot(1,2,1);pcolor(flipud(M1)), shading interp, colormap(gray),axis off; 
title('first mode by SVD (U matrix)');
subplot(1,2,2);pcolor(flipud(M2)), shading interp, colormap(gray),axis off;
title('first mode by PowerIter');

disp('the gap between sigma and sqrt(lambda):')
disp(abs(CS(1,1) - sqrt(lambda)));
%% (b) Use randomized sampling to reproduce the SVD matrices
[OU,OS,OV] = RSVD(CImgData,500);
M3 = reshape(OU(:,1), CImgHeight,CImgWidth);
figure(9);
pcolor(flipud(M3)), shading interp, colormap(gray),axis off;
title('first mode by Randomized Sampling');

%% (c)leading mode
err = zeros(1,11);
k = 550:-50:50;
figure(10);
subplot(3,4,1);pcolor(flipud(M1)), shading interp, colormap(gray),axis off;
title('SVD');
hold on
for i = 1:11
   [U,~,~] = RSVD(CImgData,k(i));
   err(i) = norm(u1+abs(U(:,1)));
   M = reshape(-abs(U(:,1)), CImgHeight,CImgWidth);
   subplot(3,4,i+1);pcolor(flipud(M)), shading interp, colormap(gray),axis off;
   title(['k = ' num2str(k(i))]);
end
hold off
figure(11);
plot(k, err);
xlabel('k');
ylabel('error of the leading u');
%% singular value decay
figure(12);
CSigma = diag(CS);
CSigma = CSigma(CSigma > 1);
semilogy(1:600,CSigma(1:600),'ob');
hold on
for k = 50:50:550
   [~,S,~] = RSVD(CImgData,k);
   S = diag(S);
   S = S(S > 1);
   semilogy(1:size(S),S);
end
legend('svd','k=50','k=100','k=150','k=200','k=250','k=300','k=350','k=400','k=450','k=500','k=550');
ylabel('\sigma');
