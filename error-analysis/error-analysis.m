%% Initial values
x = 1;
k = 0:16;
%% h inputs
h = 10.^(-k);
%% finite-diff derivative fnct subtracted by the true sec^2(x) value
m = mag(x,h);
m = abs(m);
%% log scale plot
loglog(h,m)
title('Magnitude Error of Finite Difference Method')
xlabel('h values')
ylabel('Magnitude Error')