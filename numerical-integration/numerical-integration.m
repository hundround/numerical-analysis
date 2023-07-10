%% I. Composite Simpsons Rule
a = 1;
b = 2;
m = 4;
f1 = @(x) log(x);

format long

k = compsimp(f1,m,a,b)

%% II. Bessel Function
a = 0;
b = pi;

A = zeros(101,1);
B = zeros(101,1);
C = zeros(101,1);

xaxis = linspace(0,20,101);
for i = 0:100
    j = i*0.2;
    f =@(x) (1/pi)* cos(j*sin(x));
    A(i+1) = compsimp(f,10,a,b);
end
plot(xaxis,A,'b-')

hold on

for i = 0:100
    j = i*0.2;
    g =@(x) (1/pi)* cos(j*sin(x) - x);
    B(i+1) = compsimp(g,10,a,b);
end
plot(xaxis,B,'r-')

for i = 0:100
    j = i*0.2;
    h =@(x) (1/pi)* cos(j*sin(x) - 2*x);
    C(i+1) = compsimp(h,10,a,b);
end
plot(xaxis,C,'g-')

hold off
legend('$J_0(x)$','$J_1(x)$','$J_2(x)$','Interpreter','latex')

%% III. Three-point center difference formula
fb =@(x) (x-1)^4 + 8;
x0 = 0;
tol = 10e-7;
max_iter = 1000;
h = 10e-4;

x = x0;
iter = 0;

while iter < max_iter
    fpr1 = (fb(x+h)-fb(x-h))/(2*h);
    fpr2 = (fb(x+h)-2*fb(x)+fb(x-h))/(h^2);
    if (abs(fpr1) < tol) && (fpr2 > 0)
        break;
    else
        xdelta = -(fpr1/fpr2);
        x = x + xdelta;
        iter = iter + 1;
    end
end
x_output = x
