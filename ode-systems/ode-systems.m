%% I.
close
clear all
format long
% a. Show that y_exact is an exact solution to the ODE.
% Written on the paper

% b. Convert to a system of ODEs
% Written on the paper

% c. Solve the ODE using RK Method.
syms y y1b1(t) y1b2(t) y1b3(t)
y1b1 = diff(y,t,1);
y1b2 = diff(y1b1,t,1);

f=@(t,y) y1b1 + t;
y0=0;
a=0;
b=10;
N=10;
h=(b-a)/N;
t=a:h:b;

N=length(t)-1;
w(1)=y0;
for i=1:N
    s1 = f(t(i),w(i));
    s2 = f(t(i)+h/2,w(i)+(h/2)*s1);
    s3 = f(t(i)+h/2,w(i)+(h/2)*s2);
    s4 = f(t(i)+h,w(i)+h*s3);
    w(i+1) = w(i) + (h/6)*(s1+2*s2+2*s3+s4);
end

plot(t,w,'--o');
hold on

[t,y] = ode45(@ode1b,[0 10],[0;0;0]);
plot(t,y,'r-')
title('y_{exact} vs. Solution using RK Method');
xlabel('t');
ylabel('y');
legend('y_{exact}','RK Method')
% % Plot of the exact solution y_exact
% fplot(@(t) y_exact1(t),[0 10],'b-');
% ylim([0 50]);
% legend('y_{exact}');
% hold off;

%% II
close
clear all

% a. Solve the exact solution to the ODE
% Written on the paper

% b. Error plot
y_euler = zeros(6,1);
error = zeros(6,1);
k = 0:5;
h=0.1.*2.^(-k);
for i=1:6
    len=0:h(i):1;
    [t,y]=euler1([0 1],0,length(len),h(i));
    y_euler(i) = y(10*2^(i-1)+1);
    error(i) = abs(y_exact2(1)-y_euler(i));
end

loglog(h,error,'b-')
xlabel('h')
ylabel('error')

%% III.
close
clear all

% % a. Plot the graph of the population
a1 = 1;
a2 = 0.5;
b1 = 0.1;
b2 = 0.02;
tspan = [0 50];
[t,y] = ode45(@(t,y) predpray(t,y,a1,a2,b1,b2),tspan,[100;10]);
plot(t,y(:,1),'b-',t,y(:,2),'r-');
legend('y_1 (Prey)','y_2 (Predator)')
xlabel('t');
ylabel('Population')

figure;
% b. Plot (y_1(t),y_2(t))
plot(y(:,1),y(:,2),'r-');
hold on;
xlabel('Prey');
ylabel('Predator');

% Explanation is written on the paper.
% c. What initial population would extinct the prey/predator, if possible?
% Written on the paper
%% Function Dump

function y = y_exact1(t)
    y = (exp(t) + exp(-t)-t^2)/2;
end

function y = y_exact2(t)
    y = (t - 1 + exp(-t));
end
function dydt = predpray(t,y,a1,a2,b1,b2)
    dydt = zeros(2,1);
    dydt(1) = y(1)*(a1-b1*y(2));
    dydt(2) = y(2)*(-a2+b2*y(1));
end 

function dydt = ode1b(t,y)
dydt = [y(1); y(2); y(1) + t];
end

function [t,y]=euler1(inter,y0,n,h)
t(1)=inter(1); y(1)=y0;
for i=1:n
t(i+1)=t(i)+h;
y(i+1)=eulerstep(t(i),y(i),h);
end
end

function y=eulerstep(t,y,h)
y=y+h*ydot(t,y);
end

function z=ydot(t,y)
%right-hand side of differential equation
z=t-y;
end
