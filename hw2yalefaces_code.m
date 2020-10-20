clc,clear,close all;
%MainDir = 'H:\UW course work\amath 584\hw2\code';
%% import cropped images
CDirNum = 37; % the images in the last folder are used as test images
CDirIndex = [1:13 15:39]; % I didn't notice in the first place that yaleB14 does not exist
CImgNum = 64;
CImgTest = double(imread('CroppedYale\YaleB39\yaleB39_P00A+000E+00.pgm')); % one of the test images
[CImgHeight, CImgWidth]= size(CImgTest);
CRawImg = zeros(CImgHeight, CImgWidth, CImgNum*CDirNum);
Ccount = 0;
CImgAvg = zeros(CImgHeight, CImgWidth); % average face
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
            %pcolor(flipud(RawImg(:,:,k))), shading interp, colormap(gray)
            %pause(0.1); % showing all the analysis images
        end  
    end
end
CImgAvg = CImgAvg./Ccount;
%% import uncropped data
UImgDir = ['yalefaces\'];
UImgNum = 164;
UImgTest = double(imread([UImgDir 'subject01.centerlight']));
[UImgHeight, UImgWidth]= size(UImgTest);
URawImg = zeros(UImgHeight, UImgWidth, UImgNum);
Ucount = 0;
UImgAvg = zeros(UImgHeight, UImgWidth); % average face
    UImgDirName = dir(UImgDir);
for j = 4: length( UImgDirName )
    if(isequal(UImgDirName(j).name, '.')||...
       isequal(UImgDirName(j).name, '..'))
       continue;
    else
        Ucount=Ucount+1;
        URawImg(:,:,Ucount) = imread([UImgDir,'\',UImgDirName(j).name]);
        UImgAvg = UImgAvg + URawImg(:,:,Ucount);
        %pcolor(flipud(RawImg(:,:,k))), shading interp, colormap(gray)
        %pause(0.1);
    end  
end
UImgAvg = UImgAvg./Ucount;
%% Show 20 uncropped & cropped images
ShowNum = 20;
figure(1);
for i=1:ShowNum
   subplot(4,5,i);
   pcolor(flipud(CRawImg(:,:,randi(CImgNum*CDirNum)))), shading interp, colormap(gray),axis off;
end
suptitle('Some cropped images');
figure(2);
for i = 1:ShowNum
   subplot(4,5,i);
   pcolor(flipud(URawImg(:,:,randi(UImgNum)))), shading interp, colormap(gray),axis off;
end
suptitle('Some uncropped images');
%% plot the average face
figure(3);
subplot(1,2,1),pcolor(flipud(CImgAvg)), shading interp, colormap(gray);
title('The average face of all the cropped images(2368 faces)');

subplot(1,2,2),pcolor(flipud(UImgAvg)), shading interp, colormap(gray);
title('The average face of all the uncropped images(164 faces)');
%% Center the sample pictures at the "origin"
%cropped
CImgAvg = reshape(CImgAvg, [CImgHeight*CImgWidth, 1])*ones(1,CImgNum*CDirNum);
CImgData = squeeze(reshape(CRawImg,[CImgHeight*CImgWidth,1,CImgNum*CDirNum]));
CImgData = CImgData - CImgAvg;
%uncropped
UImgAvg = reshape(UImgAvg, [UImgHeight*UImgWidth, 1])*ones(1,UImgNum);
UImgData = squeeze(reshape(URawImg,[UImgHeight*UImgWidth,1,UImgNum]));
UImgData = UImgData - UImgAvg;
%% SVD decomposition
[CU,CS,CV] = svd(CImgData, 'econ');
[UU,US,UV] = svd(UImgData, 'econ');

%% plot singular value spectrum
%cropped
figure(4);
CSigma = diag(CS);
CSigma = CSigma(CSigma > 1);
semilogy(1:size(CSigma),CSigma,'ob-');
title("singular value spectrum of cropped images");
xlabel('p') ;
ylabel('singular value \sigma_p');
%uncropped
figure(5);
USigma = diag(US);
USigma = USigma(USigma > 1);
semilogy(1:size(USigma),USigma,'or-');
title("singular value spectrum of uncropped images");
xlabel('p') ;
ylabel('singular value \sigma_p');
%% compare of the singular values of cropped and uncropped images
figure(6);
semilogy(1:size(USigma),CSigma(1:size(USigma),:),'b-');
hold on
semilogy(1:size(USigma),USigma,'r-');
hold off
title("singular value of cropped&uncropped images");
xlabel('p') ;
ylabel('singular value \sigma_p');
legend('cropped images','uncropped images');
%% Show some eigenfaces (modes)
num = 20;
CEigFaces= zeros(CImgHeight*CImgWidth,num,1);
CEigFaces(:,:,1) = CU(:,1:num);
CEigFaces = reshape(permute(CEigFaces,[1,3,2]),[CImgHeight,CImgWidth,num]);
figure(7);
for i = 1:num
    subplot(4,5,i);pcolor(flipud(CEigFaces(:,:,i))), shading interp, colormap(gray),axis off;
end
suptitle("The first 20 eigenfaces of cropped images");

UEigFaces= zeros(UImgHeight*UImgWidth,num,1);
UEigFaces(:,:,1) = UU(:,1:num);
UEigFaces = reshape(permute(UEigFaces,[1,3,2]),[UImgHeight,UImgWidth,num]);
figure(8);
for i = 1:num
    subplot(4,5,i);pcolor(flipud(UEigFaces(:,:,i))), shading interp, colormap(gray),axis off;
end
suptitle("The first 20 eigenfaces of uncropped images");
%% compare the first eigen face
figure(9);
subplot(1,2,1);pcolor(flipud(CEigFaces(:,:,1))), shading interp, colormap(gray),axis off
title('the first cropped eigenface');
subplot(1,2,2);pcolor(flipud(UEigFaces(:,:,1))), shading interp, colormap(gray),axis off;
title('the first uncropped eigenface');
%suptitle("The first principal component of cropped&uncropped faces");
%% using eigenfaces to reconstruct an image
% decompose&reconstruct a random cropped image
%CRandImg = CRawImg(:,:,randi(CImgNum*CDirNum,1));
figure(10);
subplot(2,5,1),pcolor(flipud(CImgTest)), shading interp, colormap(gray),title('test cropped image'),axis off;
CImgTest = double(reshape(CImgTest,[CImgHeight*CImgWidth, 1])) - CImgAvg(:,1);
subplotIndex = 2;
for r=[10 20 40 80 160 320 640 1280 2000]
    reconFace = CImgAvg(1) + (CU(:,1:r)*(CU(:,1:r)'*CImgTest));
    subplot(2,5,subplotIndex),pcolor(flipud(reshape(reconFace,CImgHeight,CImgWidth))),shading interp, colormap(gray),title(['r=' num2str(r)]),axis off;
    subplotIndex=subplotIndex+1;
end
% decompose&reconstruct a random uncropped image
%UImgTest = URawImg(:,:,randi(UImgNum,1));
figure(11);
subplot(2,5,1),pcolor(flipud(UImgTest)), shading interp, colormap(gray),title('test uncropped image'),axis off;
UImgTest = double(reshape(UImgTest,[UImgHeight*UImgWidth, 1])) - UImgAvg(:,1);
subplotIndex = 2;
for r=[30 45 60 75 90 105 120 135 150]
    reconFace = UImgAvg(1) + (UU(:,1:r)*(UU(:,1:r)'*UImgTest));
    subplot(2,5,subplotIndex),pcolor(flipud(reshape(reconFace,UImgHeight,UImgWidth))),shading interp, colormap(gray),title(['r=' num2str(r)]),axis off;
    subplotIndex=subplotIndex+1;
end
%% Images from 2 different people in the coordinates of 2 principle components
% a single person has 64 cropped images
CP = CImgData(:,1:192);
%RandIndex = randi(320,2,12); [RandIndex(1,i) RandIndex(2,i)]
figure(12);
for i=1:12
    subplot(3,4,i);
    CProjP = CU(:,(i-1)*10+1:(i-1)*10+2)'*(CP - CImgAvg(:,1:192));
    plot(CProjP(1,1:64),CProjP(2,1:64),'ob'); hold on
    plot(CProjP(1,65:128),CProjP(2,65:128),'xr'); 
    %plot(CProjP(1,129:192),CProjP(2,129:192),'>k'); 
    title(['r=' num2str((i-1)*10+1) '&' num2str((i-1)*10+2)]);
    %legend('person 1','person 2','Location','northwest');
    hold off
    %pause(1);
end
suptitle('o and x represent sets of images of 2 different people (cropped) ');
%% uncropped images
% a single person has 21 cropped images
UPNum = 21;
UP = UImgData(:,1:3*UPNum);
%RandIndex = randi(320,2,12); [RandIndex(1,i) RandIndex(2,i)]
figure(13);
for i=1:12
    subplot(3,4,i);
    UProjP = UU(:,(i-1)*10+1:(i-1)*10+2)'*(UP - UImgAvg(:,1:3*UPNum));
    plot(UProjP(1,1:UPNum),UProjP(2,1:UPNum),'ob'); hold on
    plot(UProjP(1,UPNum+1:2*UPNum),UProjP(2,UPNum+1:2*UPNum),'xr'); 
    %plot(CProjP(1,129:192),CProjP(2,129:192),'>k'); 
    title(['r=' num2str((i-1)*10+1) '&' num2str((i-1)*10+2)]);
    %legend('person 1','person 2','Location','northwest');
    hold off
    %pause(1);
end
suptitle('o and x represent sets of images of 2 different people (uncropped) ');


% %% single person experiment
% % do the SVD analysis using a single person's images
% SImgData = CImgData(:,1:64);
% [SU,SS,SV] = svd(SImgData ,'econ');
% figure(4);
% SSigma = diag(SS);
% SSigma = SSigma(SSigma > 1);
% semilogy(1:size(SSigma),SSigma,'*-');
% title("singular value spectrum of the first person");
% xlabel('p') ;
% ylabel('singular value \sigma_p');
