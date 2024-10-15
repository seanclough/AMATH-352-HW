r = RandStream('mt19937ar','Seed',1234);
A = r.randn(6,6);
[A1,A2] = problem1(A);
L = tril(A2) - diag(diag(A2)) + eye(length(A2));
U = triu(A2);
error = A(A1,:)-L*U;

b = [1;0;1;0;1;0];
A3 = b(A1,:);
temp1 = L;
temp1(:,length(temp1)+1) = A3;
A4 = Forsub(temp1);
temp1 = U;
temp1(:,length(temp1)+1) = A4;
A5 = Backsub(temp1);

A = [1.1,0.2,-0.2,0.5;
     0.2,0.9,0.5,0.3;
     0.1,0.0,1.0,0.4;
     0.1,0.1,0.1,1.2];
b = [1;0;1;0];
ANS = [1,1,1];
ERR = [1,1,1];
max_error = [0.01, 0.0001, 0.000001];
for i = 1:3
    for T = 1:25
        y = problem3(A,b,T);
        error = A*y - b;
        temp1 = true;
        temp2 = 0;
        for j = 1:length(error)
            if abs(error(j)) > temp2
                temp2 = abs(error(j));
            end
        end
        if temp2 < max_error(i)
            ERR(i) = temp2;
        else
            temp1 = false;
        end
        if temp1 == true
            ANS(i)=T;
            break
        end
    end
end
A6 = ANS(1);
A7 = ERR(1);
A8 = ANS(2);
A9 = ERR(2);
A10 = ANS(3);
A11 = ERR(3);

function y = problem3(A,b,T)
    n = length(A);
    y = zeros(n,1);
    M = eye(n)-A;
    for i = 1:T
        y = M*y + b;
    end
end

function [p, M] = problem1(A_)
    n = length(A_);
    p = (1:n)';
    M = A_;
    for i = 1:n-1
        for m = i:n
            for k = m+1:n
                if abs(M(k,i)) >= abs(M(m,i))
                    temp = M(i,:);
                    M(i,:) = M(k,:);
                    M(k,:) = temp;
                    temp = p(i);
                    p(i) = p(k);
                    p(k) = temp;
                end
            end
        end
        for j = i+1:n
            a = M(j,i)/M(i,i);
            M(j,i:n) = M(j,i:n) - a*(M(i,i:n));
            M(j,i) = a;
        end
    end
end