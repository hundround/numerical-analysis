function c=compsimp(f,m,a,b)
    k1 = linspace(a,b,m+1);
    h = (b-a)/(2*m);
    k2 = linspace(a+h,b-h,m);
    c = (h/3).*(2.*sum(f(k1(2:m)))+f(a)+f(b)+4.*sum(f(k2(1:m))));
end