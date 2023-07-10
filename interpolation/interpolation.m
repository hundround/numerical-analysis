%% I. Interpolant using Monomial Basis
% x and y inputs
a = [-2, 0, 1];
b = [-27, -1, 0];
m = monocoeff(a,b);

%% II. Interpolants using a. mono basis, b. linear spline, c. cubic spline
format long 

f=@(k) 1./(1 + 25*k.^2);
x = linspace(-1,1,100);
y = f(x);


% Data points
x1 = linspace(-1,1,12);
y1 = f(x1);
plot(x1,y1,'o','Color','black')
ylim([-0.2 1])
xlim([-1 1]);
xticks(-1:0.5:1);

hold on 

% Monomial basis
MonoM = monocoeff(x1,y1);
disp(MonoM);
plot(x,polyval(MonoM,x),'b-')

% Linear Spline
%xq1 = linspace(-1,1,12);
%ls1 = interp1(xq1,f(xq1),x,'linear');
%plot(xq1,f(xq1),'-',x,ls1,'Color','black');

% Cubic Spline
%xq2 = linspace(-1,1,12);
%ls2 = interp1(xq2,f(xq2),x,'v5cubic');
%plot(xq2,f(xq2),'--',x,ls2,'Color','red');
cub = splineplot(x1,y1,12);

hold off
legend('data','monomial','linear','cubic')
