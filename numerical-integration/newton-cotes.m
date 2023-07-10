%input: function f
%       n=1 for (open) midpoint rule
%       n=2 for (closed) trapezoid rule 
%       n=3 for (closed) simpson's rule
%       [a,b] interval of integration

function y=newtoncotes(f,n,a,b)

if n==1 %midpoint
    h=b-a;
    y=h*f(a+h);
elseif n==2 %trapezoid rule
    h=b-a;
    y=(h/2)*(f(a)+f(b));
elseif n==3 %simpson's rule
    h=(b-a)/2;
    y=(h/3)*(f(a)+4*f(a+h)+f(b));
end