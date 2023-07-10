%% I.
%% 1 Calculate the norm %%
b = [0.254; 0.127];
a = [0.913, 0.659; 0.457, 0.330];
x1 = [0.6391; -0.5];
x2 = [0.999; -1.001];
disp(b);
disp(a);
r1 = b - (a*x1);
r2 = b - (a*x2);
k1 = norm(r1,1);
k2 = norm(r2,1);
disp(k1);
disp(k2);

%% 2 Calculate the condition number %%
acond = cond(a,1)

%% 3 Comparison %%
anorm = norm(a,1);
x1norm = norm(x1,1);
x2norm = norm(x2,1);
M1 = acond * k1/(anorm * x1norm);
M2 = acond * k2/(anorm * x2norm);
disp(M1);
disp(M2);
fprintf('No comparison can be made for x1 and x2 given that the condition number with value %d is relatively high which makes the matrix A ill-conditioned.\n' , acond);

%% II. 
%% 1 Determine the orthogonal Q and upper triangular R %%
A = [3, 1, 2; 6, 3, 4; 3, 1, 5];
[Q,R] = qr(A);
disp(Q);
disp(R);

%% 2 Determine y using gaussian elim w/ back substitution %%

a = Q; %define matrix a
b = [0, 1, 3]'; %define RHS vector b

n = length(b);
eps = 1e-8;

%elimination

for j = 1 : n-1
    if abs(a(j,j))<eps; error('zero pivot encountered'); end
    for i = j+1 : n
        mult = a(i,j)/a(j,j);
        for k = j:n
            a(i,k) = a(i,k) - mult*a(j,k);
        end
        b(i) = b(i) - mult*b(j);
    end
end

u = a; % backward susbtitution
disp(a);
disp(b);
for i=n:-1:1
    for j=i+1:n 
        b(i)=b(i)-u(i,j)*y(j);
    end
    y(i)=b(i)/u(i,i);
end
disp(y)

%% 3 Calculate x using back substitution %%

u2 = R ; %define u2
c = y'; %define c
n = length(c);

for i=n:-1:1
    for j=i+1:n 
        c(i)=c(i)-u2(i,j)*x(j);
    end
    x(i)=c(i)/u2(i,i);
end
disp(x);
