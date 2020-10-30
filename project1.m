% Assignment:
% Generates matrix A and b for any N and solves each iteratively
    % write uij as a vector and then Au=b is a system of equations
    % We find A and b for any N
% Jacobi (Kartik)
% Gauss-Seidel with over-relaxation parameter (Saket)
% Conjugate Gradient method (Heather)

% Test from lecture to see if it works
A = [7 3 1; 3 10 2; 1 2 15];
b = [28 31 22]';
x = [0 0 0]';
x = conjGradMethod(A, b, x);

% Conjugate Gradient Method
% Takes in matrix A, matrix b, and solves for our solution x
function x = conjGradMethod(A, b, x)
    % use a variable for tolerance to easily change it
    tolerance = 1e-12;
    
    % initialize r to -1/2(gradient)
    r = b - A*x;
    
    % initialize p
    p = r;
    % we need to use r'r for k and k+1, so keep track of both
    old = r' * r;
    new = old;
   
    % n iterations
    while(sqrt(old) >= tolerance)
        % calculate old/p'Ap to avoid multiple calculations
        oldpAp = old/(p' * A * p);
        % x = x + ((p'r)/(p'Ap))p = x + ((r'r)/p'Ap)p
        x = x + oldpAp * p;
        r = r - oldpAp * A * p;
        
        % new r'r, to use while we still need the old one
        new = r' * r;
         % return if norm is less than tolerance already
        if(sqrt(new) < tolerance)
            % print out our answer
            fprintf('Result of conjugate gradient method: ');
            % x as our result (may switch to return value as needed)
            display(x);
            return;
        end
        % p = r - ((r'Ap/p'Ap))p = r + new(r'r)/old(r'r) * p
        p = r + (new/old) * p;
        
        % make sure to update our r'r
        old = new;

    end
end
