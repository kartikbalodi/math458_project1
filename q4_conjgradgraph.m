% Assignment:
% Generates matrix A and b for any N and solves each iteratively
    % write uij as a vector and then Au=b is a system of equations
    % We find A and b for any N
% Jacobi (Kartik)
% Gauss-Seidel with over-relaxation parameter (Saket)
% Conjugate Gradient method (Heather)

N = 149; % N

% initialize b to g(xij,  yij)
h = 1/(N+1);
b = zeros(N^2,1);
for i = 1:N
    for j = 1:N
        if i/(N+1) > 1/5 && i/(N+1) < 3/5 && j/(N+1) > 1/4 && j/(N+1)<1/2
            b((j-1)*N+i,1) = -h^2; % bij = -h^2*g(xij,yij
        end
    end
end

% filling in A
A = zeros(N^2,N^2);
for j = 1:N
    for i = 1:N
        % 
        %    B I 0           4 1 0
        % A= I B I where B = 1 4 1 , I is 3x3 negative identity matrix and
        %    0 I B           0 1 4   0 is 03x3 zero matrix
        % we use this to get 
        % 4ui;j - ui+1;j - ui-1;j - ui;j+1 - ui;j-1
        % fill in identity matrices on the side diagonals
        c = (j-1)*N;
        if j > 1
            A(c+i-N, c+i) = -1;
        end
        if j < N
            A(c+i+N, c+i) = -1;
        end
        % then fill in B
        A(c+i,c+i)= 4;
        if i ~= 1
            A(c+i, c+i-1) = -1;
        end
        if i ~= N
            A(c+i, c+i+1) = -1;
        end
        
    end
end

x0 = zeros(N^2, 1); % initial guess
conjGradMethod(A, b, x0)

% Conjugate Gradient Method
% Takes in matrix A, matrix b, and solves for our solution x
function x = conjGradMethod(A, b, x)
    % use a variable for tolerance to easily change it
    programTimeLimit = 150;
    tolerance = 1e-3;
    figure(1);
    ax = axes();
    hold(ax);
    xlabel(ax, 'iterations');
    ylabel (ax, 'residual');
    iterCount = 0;
    tic
    % initialize r to -1/2(gradient)
    r = b - A*x;
    
    % initialize p
    p = r;

    % n iterations
    while(norm(r) >= tolerance)
        norm1 = norm(r);
        plot(ax,iterCount,norm1,'kx');
        time = toc;
        title(time);
        if time > programTimeLimit
            disp('Program time limit '+programTimeLimit+' reached');
            break
        end
        drawnow()
        
        % calculate Ap and p' * Ap to avoid multiple calculations
        Ap = A * p;
        pAp = p' * Ap;
        x = x + (r' * r)/(pAp) * p;
        r = r - (p' * r)/(pAp)*Ap;
        
         % return if norm is less than tolerance already
        if(sqrt(r' * r) < tolerance)
            % print out our answer
            return;
        end
        % calculate our new p
        p = r - (r'*Ap)/(pAp) * p;
        iterCount = iterCount + 1;

    end
end



