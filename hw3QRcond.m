%%
clc;close all;clear;
%% 1 modified gram schmit vs qrfactor vs QR
% Let's fisrt show the correction of our modified gram schmidt algorithm
TestMatrix = [magic(5);ones(1, 5)];
[Q,R] = MGramSchmidt(TestMatrix);
% Q is orthogonal
eq = norm(eye(5) - Q'*Q);
disp('the orthogonality of Q of the test matrix:'); disp(eq);
% QR = A
er = norm(Q*R - TestMatrix);
disp('the error of reconstructing A:'); disp(er);

%% experient on a Hilbert Matrix
H = hilb(15); % A Hillbert matrix is a typical ill-conditioned matrix
H = H(:,end-9:end); % make sure that this is a over-determined system
[m,n] = size(H);
%Show that Matrix H is ill-conditions
disp('condition number of H:');disp(cond(H));
disp('rank of H:'); disp(rank(H));

%QR decomp using 3 methods (reduced QR)
[Q1,R1] = MGramSchmidt(H);
[Q2,R2] = qrfactor(H);
[Q3,R3] = qr(H,0);
Q2 = Q2(1:m,1:n);
R2 = R2(1:n,:);
%How the algorithm works
% orthogonality of Q
eq1 = norm(eye(n) - Q1'*Q1);
eq2 = norm(eye(n) - Q2'*Q2);
eq3 = norm(eye(n) - Q3'*Q3);
disp('The 2 norm of I - Q''Q');
disp(strcat('Modified Gram Schmidt:',num2str(eq1)));
disp(strcat('qrfactor:',num2str(eq2)));
disp(strcat('reduced qr of Matlab:',num2str(eq3)));
% recontructing ability
er1 = norm(H - Q1*R1);
er2 = norm(H - Q2*R2);
er3 = norm(H - Q3*R3);
disp('The Frobenius norm of H - QR');
disp(strcat('Modified Gram Schmidt:',num2str(er1)));
disp(strcat('qrfactor:',num2str(er2)));
disp(strcat('reduced qr of Matlab:',num2str(er3)));
%% experiment on some random ill-conditioned matrix
Num = 20; noise = 0.01;
CondN = zeros(1, Num);
ON = zeros(3, Num); 
RN = zeros(3, Num);
for j = 1:Num
    S = rand(200, 99);  
    S = [S noise*S(:,1)];
    CondN(j) = cond(S);
    [Q1,R1] = MGramSchmidt(S);
    [Q2,R2] = qrfactor(S);Q2 = Q2(:,1:100);R2 = R2(1:100,:);
    [Q3,R3] = qr(S,0);
    ON(:,j) = [norm(eye(100) - Q1'*Q1);norm(eye(100) - Q2'*Q2)...
        ;norm(eye(100) - Q3'*Q3)];
    RN(:,j) = [norm(S - Q1*R1);norm(S - Q2*R2);norm(S - Q3*R3)];
end
figure(1);
plot(CondN,ON(1,:),'b*'); hold on
plot(CondN,ON(2,:),'r>'); 
plot(CondN,ON(3,:),'ko'); hold off
xlabel('condition number') ;
ylabel('norm of I - Q^*Q');

legend('Modified GS','qrfactor.m','qr of Matlab');
figure(2);
plot(CondN,RN(1,:),'b*'); hold on
plot(CondN,RN(2,:),'r>');
plot(CondN,RN(3,:),'ko'); hold off
legend('Modified GS','qrfactor.m','qr of Matlab');
xlabel('condition number') ;
ylabel('norm of A - QR');

OErrMGS = mean(ON(1,:)); RErrMGS = mean(RN(1,:));
OErrQRF = mean(ON(2,:)); RErrQRF = mean(RN(2,:));
OErrQRM = mean(ON(3,:)); RErrQRM = mean(RN(3,:));
disp('avg error of I-Q^*Q');disp([OErrMGS,OErrQRF,OErrQRM]);
disp('avg error of A-QR');disp([RErrMGS,RErrQRF,RErrQRM]);
%% 2 ploynomial
p1 = @(x)(x-2).^9;
p2 = @(x) x.^9 - 18*x.^8 + 144*x.^7 - 672*x.^6 + 2016*x.^5 ...
    - 4032*x.^4 + 5376*x.^3 - 4608*x.^2 + 2304*x - 512;
x = 1.920:0.001:2.080; % the interval we wanna show
figure(3);
plot(x,p1(x),'b');
hold on
plot(x,p2(x),'r');
%axis equal;
%xlim([-0.01 0.01]);
legend('left side', 'right side');
hold off


%% 3 conditioning of a matrix
%(a)
Num = 100; %generate 100 random matrices
condN = zeros(1,Num); % to store all the condition numbers 
for n = 1 : Num
    m = 2 * n; % set m = 2n to let A be a overdetermined system
    A = rand(m, n);
    condN(n) = cond(A); 
end
figure(4);
plot(1:Num, condN);
xlabel('column numbers n (m = 2n)') ;
ylabel('condition number');
%(b)
c1 = cond(A);    % the condition number of A
a1 = A(:,1);     % the fisrt colume of A
A1 = [A a1];     % append a1 as the (n+1)th column of A
c2 = cond(A1);   % condition number afterwards
disp('the condition number of the original matrix and the appended matrix A');
disp(c1);disp(c2);

B = rand(200, 199);
b1 = B(:,1);
B1 = [B, b1];
cb1 = cond(B); cb2 = cond(B1); detB1 = det(B1); 
disp('the condition number of the original matrix and the appended matrix B');
disp(cb1);disp(cb2);
disp('the determinant of B1');
disp(detB1);
%(c)
c = [ ];
for epsilon = 0.001: 0.001:0.1
    a1 = A(:,1);
    a1 = a1 + epsilon * rand(m, 1);
    A1 = [A a1];
    c = [c , cond(A1)];
end
figure(5);
plot(0.001: 0.001:0.1, c);
xlabel('\epsilon') ;
ylabel('conditional number');