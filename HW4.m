r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,6);
Q = eye(length(A(:,1)));
R = A;
for i = 1:length(A(1,:))
    w = R(i:end,i);
    w(1) = w(1) + sign(w(1))*norm(w);
    w = w/norm(w);
    R(i:end,i:end) = R(i:end,i:end) - 2*w*(w'*R(i:end,i:end));
    Q(:,i:end) = Q(:,i:end) - 2*(Q(:,i:end)*w)*w';
end
A1 = Q;
A2 = R;

r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,6);
b = r.randn(10,1);
A = [0,2,1;1,1,-1;2,1,0;1,1,1;0,2,-1];
b = [1;0;1;-1;0]
A3 = A'*A;
A4 = A'*b;
A5 = A3\A4
A6 = abs(b-A*A5);

r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,10)+1i*r.randn(10,10);
[V,D] = eig(A);
for i = 1:length(D(1,:))
    for j = i+1:length(D(1,:))
        if real(D(j,j)) < real(D(i,i))
            temp = D(i,i);
            D(i,i) = D(j,j);
            D(j,j) = temp;
            temp = V(:,i);
            V(:,i) = V(:,j);
            V(:,j) = temp;
        end
    end
end
%V*D*inv(V)-A
A7 = real(V);
A8 = imag(V);
A9 = real(D);
A10 = imag(D);