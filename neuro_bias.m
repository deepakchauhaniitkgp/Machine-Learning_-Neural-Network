clear all;
clc;
%learning rate 
ita=0.015;
J=20;
K=20;
L=1;
I=1;
N=1;
lambda=1;
counter=0;

% num = xlsread('sn.xlsx');

num1 = 0.01:0.001:0.9;
num = num1.^2;

length1 = size(num);
length = length1(2);

error = ones(1,length);
u=rand(I,J); v= rand(J,K); w=rand(K,L);
% u = [0.5407 0.8699];
% v = [0.2648 0.1192 ; 0.3181 0.9398];
% w = [0.6465; 0.4795];
counter = 1;
pp=1;
for pp=1:100
    
% while sum1(error)>0.01
    plot( counter,sum1(error), '.');
    
    hold on ;
    counter = counter+1;
%     sum1(error)
for ii=1:length

input = num1;

x=input;
% u
% v
% w
% calculating net1
for i=1:J
    net1(i)=0;
    for j=1:I
    net1(i)=net1(i) + u(j,i)*input(ii);
    end
end
% displaying net1
% net1
for i=1:J
y(i,ii)=(1-exp(-lambda*net1(i)))/(1+exp(-lambda*net1(i)));
der_y(i,ii) = 0.5*(1-y(i,ii)^2);
end
% displaying y
% y

for i=1:K
    net2(i)=0;
    for j=1:J
        net2(i)=net2(i)+v(j,i)*y(j,ii);
    end
end

%Displaying net2
% net2
for i=1:K
z(i,ii)=(1-exp(-lambda*net2(i)))/(1+exp(-lambda*net2(i)));
der_z(i,ii) = 0.5*(1-z(i,ii)^2);
end

% Displaying z
% z

for i=1:L
    net3(i)=0;
    for j=1:K
        net3(i)=net3(i)+w(j,i)*z(j,ii);
    end
end

%Displaying net2
% net3
for i=1:L
o(i,ii)=(1-exp(-lambda*net3(i)))/(1+exp(-lambda*net3(i)));
der_o(i,ii) = 0.5*(1-o(i,ii)^2);
end

% Displaying o
% o

% output = tanh(x);
% disp( output); 

error(ii) = num(ii) - o(1,ii);
end

dww = zeros(K,L);
% Calculation of derivatives and change in error 
for ii=1:length
for j=1:L
   
    dell(j,ii) = error(ii) * der_o(j,ii);
   
    for i=1:K        
        dww(i,j) = dww(i,j) + 1/length * ita * z(i,ii)*dell(j,ii) ; 
    end
end
end


dvv=zeros(J,K);
for ii=1:length
for i=1:K
    delk(i,ii) =  dell(1,ii) * w (i,1) * der_z (i,ii);
    for j=1:J
        dvv(j,i) = dvv(j,i) + 1/length * ita * delk(i,ii)* y (j,ii);
    end
end
end

du = zeros(N,J);
for ii=1:length
for j=1:J
    sum(j)=0;
    for xx=1:K
        sum(j) = sum(j) + v(j,xx)*delk(xx,ii);
    end
    delj(j,ii)=der_y (j,ii)*sum(j);
    for i=1:N
        du(i,j) = du(i,j) + 1/length * x (ii) * delj (j,ii)*ita;
    end
end

% Display
% disp('INPUT IS '); disp(input);
% disp('initial weights u v w ');disp(u);disp(v);disp(w);
% disp('The value of y');disp(y);
% disp('The value of z'); disp(z);
% disp('Output and error'); disp(o); disp(error); 


% updated weights are follows

for i=1:I
    for j=1:J
    u(i,j) = u(i,j)+du(i,j);
    end
end
for j=1:J
    for k=1:K
    v(j,k) = v(j,k)+dvv(j,k);
    end
end
for k=1:K
    for l=1:L
    w(k,l) = w(k,l)+dww(k,l);
    end
end
end
end
