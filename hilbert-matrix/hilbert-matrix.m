H2 = zeros(2);
for i = 1:2
    for j = 1:2
        H2(i,j) = 1/(i+j-1);
    end
end

H3 = zeros(3);
for i = 1:3
    for j = 1:3
        H3(i,j) = 1/(i+j-1);
    end
end
H4 = zeros(4);
for i = 1:4
    for j = 1:4
        H4(i,j) = 1/(i+j-1);
    end
end

H5 = zeros(5);
for i = 1:5
    for j = 1:5
        H5(i,j) = 1/(i+j-1);
    end
end

H6 = zeros(6);
for i = 1:6
    for j = 1:6
        H6(i,j) = 1/(i+j-1);
    end
end

H7 = zeros(7);
for i = 1:7
    for j = 1:7
        H7(i,j) = 1/(i+j-1);
    end
end

H8 = zeros(8);
for i = 1:8
    for j = 1:8
        H8(i,j) = 1/(i+j-1);
    end
end

H9 = zeros(9);
for i = 1:9
    for j = 1:9
        H9(i,j) = 1/(i+j-1);
    end
end

H10 = zeros(10);
for i = 1:10
    for j = 1:10
        H10(i,j) = 1/(i+j-1);
    end
end
%% II. Decomposition of H to upper triangular R and lower triangular R
%% Compute for the upper triangular
R2 = chol(H2);
R3 = chol(H3);
R4 = chol(H4);
R5 = chol(H5);
R6 = chol(H6);
R7 = chol(H7);
R8 = chol(H8);
R9 = chol(H9);
R10 = chol(H10);


%% Compute for the lower triangular using the existing upper triangular
S2 = H2*(inv(R2));
S3 = H3*(inv(R3));
S4 = H4*(inv(R4));
S5 = H5*(inv(R5));
S6 = H6*(inv(R6));
S7 = H7*(inv(R7));
S8 = H8*(inv(R8));
S9 = H9*(inv(R9));
S10 = H10*(inv(R10));

%% III. Create a one-dimensional matrix b from Hx. Solve for x approx using b.
x2 = ones(2,1);
x3 = ones(3,1);
x4 = ones(4,1);
x5 = ones(5,1);
x6 = ones(6,1);
x7 = ones(7,1);
x8 = ones(8,1);
x9 = ones(9,1);
x10 = ones(10,1);

b2 = H2*(x2);
b3 = H3*(x3);
b4 = H4*(x4);
b5 = H5*(x5);
b6 = H6*(x6);
b7 = H7*(x7);
b8 = H8*(x8);
b9 = H9*(x9);
b10 = H10*(x10);


%% Solving for b using R and S %%
%% forward sub %%
l = S2 %define l
b = b2
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y2(j);
    end
    y2(i)=b(i)/l(i,i);
end

l = S3 %define l
b = b3
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y3(j);
    end
    y3(i)=b(i)/l(i,i);
end
disp(y3);

l = S4 %define l
b = b4
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y4(j);
    end
    y4(i)=b(i)/l(i,i);
end

l = S5 %define l
b = b5
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y5(j);
    end
    y5(i)=b(i)/l(i,i);
end

l = S6 %define l
b = b6
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y6(j);
    end
    y6(i)=b(i)/l(i,i);
end

l = S7 %define l
b = b7
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y7(j);
    end
    y7(i)=b(i)/l(i,i);
end

l = S8 %define l
b = b8
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y8(j);
    end
    y8(i)=b(i)/l(i,i);
end

l = S9 %define l
b = b9
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y9(j);
    end
    y9(i)=b(i)/l(i,i);
end

l = S10 %define l
b = b10
n = length(b);

for i=1:n
    for j=1:i-1 
        b(i)=b(i)-l(i,j)*y10(j);
    end
    y10(i)=b(i)/l(i,i);
end



%% backward sub %%

u = R2; %define u
c = y2; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k2(j);
    end
    k2(i)=c(i)/u(i,i);
end

u = R3; %define u
c = y3; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k3(j);
    end
    k3(i)=c(i)/u(i,i);
end

u = R4; %define u
c = y4; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k4(j);
    end
    k4(i)=c(i)/u(i,i);
end

u = R5; %define u
c = y5; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k5(j);
    end
    k5(i)=c(i)/u(i,i);
end

u = R6; %define u
c = y6; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k6(j);
    end
    k6(i)=c(i)/u(i,i);
end

u = R7; %define u
c = y7; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k7(j);
    end
    k7(i)=c(i)/u(i,i);
end

u = R8; %define u
c = y8; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k8(j);
    end
    k8(i)=c(i)/u(i,i);
end

u = R9; %define u
c = y9; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k9(j);
    end
    k9(i)=c(i)/u(i,i);
end

u = R10; %define u
c = y10; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u(i,j)*k10(j);
    end
    k10(i)=c(i)/u(i,i);
end

disp(k10)

%% IV. Computation of the infinity-norm of the residual r.
%% Computation of the residual 
r2 = b2 - H2*k2';
r3 = b3 - H3*k3';
r4 = b4 - H4*k4';
r5 = b5 - H5*k5';
r6 = b6 - H6*k6';
r7 = b7 - H7*k7';
r8 = b8 - H8*k8';
r9 = b9 - H9*k9';
r10 = b10 - H10*k10';

%% Solve for the inf-norms of residual %%
r2norm = norm(r2, inf);
r3norm = norm(r3, inf);
r4norm = norm(r4, inf);
r5norm = norm(r5, inf);
r6norm = norm(r6, inf);
r7norm = norm(r7, inf);
r8norm = norm(r8, inf);
r9norm = norm(r9, inf);
r10norm = norm(r10, inf);

%% Solve for the forward error %%
xc2 = x2 - k2;
xc3 = x3 - k3;
xc4 = x4 - k4;
xc5 = x5 - k5;
xc6 = x6 - k6;
xc7 = x7 - k7;
xc8 = x8 - k8;
xc9 = x9 - k9;
xc10 = x10 - k10;

%% Solve for the forward error inf-norm%%
xc2norm = norm(xc2, inf);
xc3norm = norm(xc3, inf);
xc4norm = norm(xc4, inf);
xc5norm = norm(xc5, inf);
xc6norm = norm(xc6, inf);
xc7norm = norm(xc7, inf);
xc8norm = norm(xc8, inf);
xc9norm = norm(xc9, inf);
xc10norm = norm(xc10, inf);

%% Condition number of H %%
H2cond = cond(H2, inf);
H3cond = cond(H3, inf);
H4cond = cond(H4, inf);
H5cond = cond(H5, inf);
H6cond = cond(H6, inf);
H7cond = cond(H7, inf);
H8cond = cond(H8, inf);
H9cond = cond(H9, inf);
H10cond = cond(H10, inf);

%% V. Plot of the errors
%% Creation of array %%
s = 2:1:10;
t = [r2norm, r3norm, r4norm, r5norm, r6norm, r7norm, r8norm, r9norm, r10norm];
v = [xc2norm, xc3norm, xc4norm, xc5norm, xc6norm, xc7norm, xc8norm, xc9norm, xc10norm];
w = [H2cond, H3cond, H4cond, H5cond, H6cond, H7cond, H8cond, H9cond, H10cond]

%% Plotting %%
subplot(3,1,1);
plot(s,v);
xlabel('n');
ylabel('forward error');

subplot(3,1,2);
plot(s,t);
xlabel('n');
ylabel('residual');

subplot(3,1,3);
semilogy(s,w);
xlabel('n');
ylabel('condition number');