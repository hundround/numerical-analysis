function M = monocoeff(x,y)
    A = vander(x);
    B = y';
    M = mldivide(A,B)';
end 