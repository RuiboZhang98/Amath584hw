function [numDigit,numCorrect] = ACU(A,X,L)
% numDigit: number of training samples of digit i
% numCorrect: number of correct predictions of A
[m,I] = max(A*X,[],1);
I(m==0) = 0;
Sol = I';
numDigit = zeros(10,1);
numCorrect = zeros(10,1);
for i = 1:10
    numDigit(i) = sum(L == i);
    numCorrect(i) = sum(min(L==Sol,L==i));    
end
numDigit = [numDigit ; sum(numDigit)];
numCorrect = [numCorrect ; sum(numCorrect)];
end

