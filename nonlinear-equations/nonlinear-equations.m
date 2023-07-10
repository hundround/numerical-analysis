%% Fixed-Point Iteration
clear
close all

%To use fpi.m, first define a function g such that x^2 - 4 = atan(2x)
g=@(x) sqrt(4 + atan(2*x));

%set initial guess and number of steps
init=2;
steps=100;

xc=fpi(g,init,steps);
fprintf('One of the solution is %d, found using Fixed-Point Iteration. \n', xc);
%% Secant Method
%Set the initial values x0 and x1
x0 = -2;
x1 = -1.5;
tol = 0.005;
N = 15;
f=@(x) atan(2*x) + 4 -x^2;

format long
[xp, i] = secant(f,x0,x1,tol,N);
fprintf('One of the solution is %d, found using Secant Method with %d iterations. \n', xp, i);
%% Plot the graph
x = -3:0.0001:3;
h=@(x) x.^2 - 4;
j=@(x) atan(2.*x);

plot(x,h(x),'LineStyle','-', 'Color', "#EDB120");
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

hold on
plot(x,j(x),'LineStyle','-','Color',"#D95319");
plot(xp,h(xp),'.', 'MarkerSize', 10, 'Color','b');
plot(xc,h(xc),'.', 'MarkerSize', 10, 'Color', 'c');

hold off
legend('$y = x^2 - 4$', '$y = \arctan(2x)$', 'Solution to the Secant Method', 'Solution to the Fixed-Point Iteration', 'Interpreter', 'latex')
