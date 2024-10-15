function x = Forsub(L)
    x = L(:,end);
    for i = 1:length(x)
        SUM = 0;
        for j=1:i-1
           SUM = SUM + L(i,j)*x(j);
        end
        x(i) = (x(i) - SUM)/L(i,i);
    end
end