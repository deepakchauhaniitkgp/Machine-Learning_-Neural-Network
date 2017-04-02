%learning rate 
ita=0.15;
J=2;
K=2;
L=1;
I=1;
N=1;
error = 1;
lambda=1;
counter=0;
num = xlsread('sn.xlsx');
length = size(num);
u=rand(I,J); v= rand(J,K); w=rand(K,L);
% u = [0.5407 0.8699];
% v = [0.2648 0.1192 ; 0.3181 0.9398];
% w = [0.6465; 0.4795];
while error>0.00001
counter = counter+1;
input =[0.65];
x=input;
% u
% v
% w
% calculating net1
for i=1:J
    net1(i)=0;
    for j=1:I
    net1(i)=net1(i)+u(j,i)*input(j);
    end
end
% displaying net1
% net1
for i=1:J
y(i)=(1-exp(-lambda*net1(i)))/(1+exp(-lambda*net1(i)));
der_y(i) = 0.5*(1-y(i)^2);
end
% displaying y
% y

for i=1:K
    net2(i)=0;
    for j=1:J
        net2(i)=net2(i)+v(j,i)*y(j);
    end
end

%Displaying net2
% net2
for i=1:K
z(i)=(1-exp(-lambda*net2(i)))/(1+exp(-lambda*net2(i)));
der_z(i) = 0.5*(1-z(i)^2);
end

% Displaying z
% z

for i=1:L
    net3(i)=0;
    for j=1:K
        net3(i)=net3(i)+w(j,i)*z(j);
    end
end

%Displaying net2
% net3
for i=1:L
o(i)=(1-exp(-lambda*net3(i)))/(1+exp(-lambda*net3(i)));
der_o(i) = 0.5*(1-o(i)^2);
end

% Displaying o
% o

% output = tanh(x);
% disp( output); 

error = input - o;
% Calculation of derivatives and change in error 
for j=1:L
    dell(j) = error * der_o(j);
    for i=1:K        
        dw(i,j) = ita * z (i)*dell(j) ; 
    end
end

for i=1:K
    delk(i) =  dell * w (i,1) * der_z (i);
    for j=1:J
        dv(j,i) = ita * delk(i)* y (j);
    end
end

for j=1:J
    sum(j)=0;
    for xx=1:K
        sum(j) = sum(j) + v(j,xx)*delk(xx);
    end
    delj(j)=der_y (j)*sum(j);
    for i=1:N
        du(i,j) = x (i) * delj (j)*ita;
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
    v(j,k) = v(j,k)+dv(j,k);
    end
end
for k=1:K
    for l=1:L
    w(k,l) = w(k,l)+dw(k,l);
    end
end
end
