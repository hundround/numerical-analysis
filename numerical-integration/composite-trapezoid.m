%input: function f
%       number of panels m
%       [a,b] interval of integration

function y=compositetrapezoid(f,m,a,b)
    h=(b-a)/m;
    x=linspace(a,b,m+1);
    y=(h/2)*(f(a)+f(b)+2*sum(f(x(2:m))));