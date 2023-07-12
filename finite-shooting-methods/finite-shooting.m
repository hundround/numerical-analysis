%% I. Finite Differences
% Using pen and paper, we solve for the finite differences formula for
% derivatives. We obtain the following expression 
% (1-h)*w_(i-1) + (-2-3h^2)*w_i + (1+h)*w_(i+1).

clear
close all
n = [9 19 39];
for i = 1:3
    h = 1/(n(i)+1);
    a = -2-3*h^2;
    b = 1+h;
    c = 1-h;
    f1 = exp(3);
    f2 = 1;

    w = g(a,b,c,n(i),f1,f2);
    t = linspace(0,1,n(i)+2);
    plot(t,w,'-o');
    hold on
end

hold on;
f=@(t) exp(3-3*t);
fplot(f,[0 1],'b--')
legend('n=9','n=19','n=39','y(t)')
ylim([0 20])
hold off;


%% I. Errors
% We need to yield the error vector whose entries are the absolute diff
% of y(t) and w.
n = 9;
h = 1/(n+1);
a = -2-3*h^2;
b = 1+h;
c = 1-h;
f1 = exp(3);
f2 = 1;
w1 = g(a,b,c,n,f1,f2);
M = linspace(0,1,n+2);
for i = 1:(n+2)
    error(i) = abs(w1(i) -f(M(i)));
end
semilogy(M,error,'b--');

hold on;
n = 19;
h = 1/(n+1);
a = -2-3*h^2;
b = 1+h;
c = 1-h;
w2 = g(a,b,c,n,f1,f2);
M = linspace(0,1,n+2);
for i = 1:(n+2)
    error(i) = abs(w2(i) -f(M(i)));
end
semilogy(M,error,'r--');

n = 39;
h = 1/(n+1);
a = -2-3*h^2;
b = 1+h;
c = 1-h;
w3 = g(a,b,c,n,f1,f2);
M = linspace(0,1,n+2);
for i = 1:(n+2)
    error(i) = abs(w3(i)-f(M(i)));
end
semilogy(M,error,'g--');
legend('n=9','n=19','n=39')
%% II. Shooting Method
clear
close all
% interval for t in y'(t)=f(t,y)
inter=[0 1]; 
ybounds=[1 -2/3];

% initial guess for s in fsolve
s0=2;
[ss,fval]=fsolve(@(s)F(s,inter),s0);
disp(ss)
disp(fval)

% plot the numerical solution y=y1 (first column) and the solution using initial guess s0
ydot=@(t,y) [y(1) - 3*y(1)*y(2);-6*(t*y(2)+log(y(1)))]; 

[t,y]=ode45(ydot,inter,[ybounds(1) ss]);
[t0,y0]=ode45(ydot,inter,[ybounds(1) s0]);
plot(t,y(:,1),'-',t,y(:,2),'-');
legend('y1','y2');
% ax=gca;
% ax.FontSize = 15;
hold on;
scatter(inter,ybounds,60,'filled')

%% Functions

function tridiag = A(a,b,c,N)
% a,b,c are diagonal entries and N is the length of the square matrix.
    tridiag = diag(a*ones(1,N)) + diag(b*ones(1,N-1),1) + diag(c*ones(1,N-1),-1);
end    

function solvector = B(b,c,N,f1,f2)
% a,b,c are diagonal entries, N is the number of rows of 1-column vector
% and f1,f2 are the initial values for y0 and y1, respectively.
    for i=1:N
        if i == 1
            solvector(i,1) = -c*f1;
        elseif i == N
            solvector(i,1) = -b*f2;
        else
            continue
        end
    end
end

function w = g(a,b,c,N,f1,f2)
% the function will solve for the linear equation 
    tridiag = A(a,b,c,N);
    solvector = B(b,c,N,f1,f2);
    w0 = tridiag\solvector;
    w = [f1;w0;f2]';
end

function z=F(s,inter)
% Example 1
init=[1 s];
yb=-2/3; % right end point in BVP
ydot=@(t,y) [y(1) - 3*y(1)*y(2);-6*(t*y(2)+log(y(1)))]; 
[t,y]=ode45(ydot,inter,init);
z=y(end,2)-yb; % end means last entry of solution y
end