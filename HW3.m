A = [1.1,.2,-.2,.5;
     .2,.9,.5,.3;
     .1,0.,1.,.4;
     .1,.1,.1,1.2];
b = [1;0;1;0];
y = zeros(length(b),1);
T_vector = zeros(2,4);
error_vector = zeros(2,4);
errors = [0.01 0.0001 0.000001 0.00000001];
%%{
for k = 1:length(errors)
    for i = 1:25
        y = A*jacobi(A,b,i)-b;
        temp = 0;
        for j=1:length(y)
            if abs(y(j)) > temp
                temp = abs(y(j));
            end
        end
        if temp < errors(k)
            T_vector(1,k) = i;
            error_vector(1,k) = temp;
            break
        end
    end
end
%}
for k = 1:length(errors)
    for i = 1:25
        y = A*gauss_seidel(A,b,i)-b;
        temp = 0;
        for j=1:length(y)
            if abs(y(j)) > temp
                temp = abs(y(j));
            end
        end
        if temp < errors(k)
            T_vector(2,k) = i;
            error_vector(2,k) = temp;
            break
        end
    end
end

A1 = T_vector(1,1);
A2 = error_vector(1,1);
A3 = T_vector(1,2);
A4 = error_vector(1,2);
A5 = T_vector(1,3);
A6 = error_vector(1,3);
A7 = T_vector(1,4);
A8 = error_vector(1,4);
A9 = T_vector(2,1);
A10 = error_vector(2,1);
A11 = T_vector(2,2);
A12 = error_vector(2,2);
A13 = T_vector(2,3);
A14 = error_vector(2,3);
A15 = T_vector(2,4);
A16 = error_vector(2,4);

%%

x0 = [ .9;    % S
       .09;   % I
       .01 ]; % R
x = x0;
p = 0;
M = [1-1/200-p, 0, 1/10000; 1/200, 1-1/1000, 0; p, 1/1000, 1-1/10000];
for i=1:1000
    x = M*x;
    if x(2) > .5
        D0 = i;
        break
    end
end
x = x0;
prev = 0;
for i=1:100000
    x = M*x;
    if abs(x(2)-prev)<10^(-8)
        F0 = x(2);
        break
    end
    prev = x(2);
end

p = 2/1000;
M = [1-1/200-p, 0, 1/10000; 1/200, 1-1/1000, 0; p, 1/1000, 1-1/10000];
x = x0;
for i=1:1000
    x = M*x;
    if x(2) > .5
        D1 = i;
        break
    end
end
x = x0;
prev = 0;
for i=1:100000
    x = M*x;
    if abs(x(2)-prev)<10^(-8)
        F1 = x(2);
        break
    end
    prev = x(2);
end

A17 = D0;
A18 = F0;
A19 = D1;
A20 = F1;

%%


function y = jacobi(A,b,T)
    n = length(b);
    y = zeros(n,1);
    L = tril(A,-1);
    D = diag(A);
    U = triu(A,1);
    for i=1:T
        y = -(L+U)*y./D+b./D;
    end
end

function y = gauss_seidel(A,b,T)
    n = length(b);
    y = zeros(n,1);
    L = tril(A,-1);
    D = diag(A);
    U = triu(A,1);
    for i=1:T
        temp1 = A-U;
        temp2 = temp1;
        temp3 = temp1;
        temp2(:,length(temp1)+1) = -(U)*y;
        temp3(:,length(temp1)+1) = b;
        y = Forsub(temp2)+Forsub(temp3);
    end
end