function  SolAnalyze(A,modename,sparseTol,X,L,Width,Height)
%  
figure(1);
for i = 1:10
    subplot(5,2,i);bar(A(i,:),5);
end
suptitle([modename ' solution'])

[numDigit,numCorrect] = ACU(A,X,L);
figure(2);
subplot(1,2,1);
xaxis = categorical({'1','2','3','4','5','6','7','8','9','0'});
xaxis = reordercats(xaxis,{'1','2','3','4','5','6','7','8','9','0'});
bar(xaxis, [numDigit(1:end-1),numCorrect(1:end-1)],'Group');
legend('number of digit i','correct number of solution');
title([modename ' prediction']);
subplot(1,2,2);
xaxis = categorical({'1','2','3','4','5','6','7','8','9','0','total'});
xaxis = reordercats(xaxis,{'1','2','3','4','5','6','7','8','9','0','total'});
acu = numCorrect./numDigit;
bar(xaxis,acu);
ylim([0 1]);
title(['Accuracy of ' modename]);
for i = 1:length(xaxis)
    text(xaxis(i),acu(i)+0.02,num2str(roundn(acu(i),-2)));
end
%% nonzero entries
j = [0,0.01,0.1,1];
nonzero = zeros(1,length(j));
for i = 1:length(j)
    nonzero(i) = sum(sum(abs(A) > j(i)));
end
figure(3);
J = categorical({'>0','>0.01','>0.1','>1'});
J = reordercats(J,{'>0','>0.01','>0.1','>1'});
bar(J,nonzero);title(['distribution of nonzero entries of ' modename]);
if sparseTol~= 0
    ylim([0,800]);
end
for i = 1:length(J)
    text(J(i),nonzero(i)+20,num2str(nonzero(i)));
end
%% nonzero pixels
j = [0,0.01,0.1,1];
AA = A;
AA1 = zeros(10,Height*Width,3);
spn = zeros(1,length(j));
for i = 1:length(j)
    Asp = A;
    Asp(abs(A)<j(i)) = 0;
    Asp = sum(Asp); %
    spn(i) = sum(Asp~=0);
    AA(:,Asp==0) = 0; 
    AA1(:,:,i)= AA; 
end
figure(4);
J = categorical({'>0','>0.01','>0.1','>1'});
J = reordercats(J,{'>0','>0.01','>0.1','>1'});
pb = bar(J,spn);title(['nonzero pixels of different sparse ' modename]);
ylim([0,800]);
for i = 1:length(J)
    text(J(i),spn(i)+20,num2str(spn(i)));
end
%% total accuracy of sparse sol
acu = zeros(1,length(j));
for i = 1:length(j)
    [numDigit,numCorrect] = ACU(AA1(:,:,i),X,L);
    acu(i) = roundn(numCorrect(11)/numDigit(11),-2);
end
figure(5);
%spn = spn + [1,2,3,4];%%%
J = categorical(split(num2str(spn)));
J = reordercats(J,split(num2str(spn)));
bar(J,acu);
title(['total accuracy of sparse ' modename]);
xlabel('number of pixels');
ylim([0 1]);
for i = 1:length(J)
    text(J(i),acu(i)+0.03,num2str(acu(i)));
end
%% per digit number of important pixels
figure(10);
psp = sum(abs(A')>sparseTol);
xaxis = categorical({'1','2','3','4','5','6','7','8','9','0'});
xaxis = reordercats(xaxis,{'1','2','3','4','5','6','7','8','9','0'});
bar(xaxis,psp);
ylim([0 800]);
xlabel('digit');ylabel('number of important pixels');
title(['number of most important pixels for each digit of ' modename]);
for i = 1:length(xaxis)
    text(xaxis(i),psp(i)+20,num2str(psp(i)));
end
%% per digit sol
SA = A;
if sparseTol~=0
    SA(abs(A)<=sparseTol) = 0;
end
figure(11);
for i = 1:10
    subplot(2,5,i);pcolor(reshape(SA(i,:),Width,Height)); colormap(hot),axis off;
    title(['digit ' num2str(i- (i ==10)*(10))]);
end
suptitle(['sparse (for digit) ' modename]);
%% per digit accuracy
[numDigit,numCorrect] = ACU(SA,X,L);
figure(12);
subplot(1,2,1);
xaxis = categorical({'1','2','3','4','5','6','7','8','9','0'});
xaxis = reordercats(xaxis,{'1','2','3','4','5','6','7','8','9','0'});
bar(xaxis,[numDigit(1:end-1),numCorrect(1:end-1)],'Group');
legend('number of digit i','correct number of solution');
title([modename ' prediction']);
subplot(1,2,2);
xaxis = categorical({'1','2','3','4','5','6','7','8','9','0','total'});
xaxis = reordercats(xaxis,{'1','2','3','4','5','6','7','8','9','0','total'});
acu = numCorrect./numDigit;
bar(xaxis,acu);
ylim([0 1]);
title(['Accuracy of sparse (for digit) ' modename]);
for i = 1:length(xaxis)
    text(xaxis(i),acu(i)+0.02,num2str(roundn(acu(i),-2)));
end
end

