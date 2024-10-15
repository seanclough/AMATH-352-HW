function x = Backsub(U)
    n = size(U); n = n(1); % get the number of rows of U
    x = U(:,n+1);
    if U(n,n) == 0
        disp('Method failed: singular matrix');
        x = Inf; % return infinity if there is no solution
    else
        x(n) = U(n,n+1)/U(n,n);
        for k = 1:n-1
            i = n-k;
            if U(i,i) == 0
                disp('Method failed: singular matrix');
                x = Inf;
                break; % return infinity if there is no solution
            else
               SUM = 0;  % compute the sum as j varies from i+1 to n
               for j = i+1:n
                   SUM = SUM + U(i,j)*x(j);
               end
               x(i) = U(i,n+1)-SUM;
               x(i) = x(i)/U(i,i);
            end
            
        end
    end     
end