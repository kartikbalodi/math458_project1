%Test matrix for function
A = [7 3 1; 
    3 10 2; 
    1 2 15];
b = [28 31 22]';
x = [0 0 0]';
%Set value of omega, can change for tests
omega = 0.1;
x = SORMethod(A,b,x, omega);

function x = SORMethod(A,b,x_init, omega)
    %Stopping Conditions: a tolerance limit and max number of iterations
    %before quitting
    tolerance = 1e-12;
    max_iterations = 100000;
  
    %LDU Decomoposition of A
    %Isolates lower diagonal (below the main diagonal)
    L = -tril(A, -1);
    %Isolates main diagonal
    D = diag(diag(A));
    %Isolates upper diagonal (above the main diagonal)
    U = -triu(A, 1);
    
    %Creates initial solution array of zeros
    x = zeros(length(A));
    
    %Stored in reference since it is used repeatedly
    %A in the equation Ax_i = b where x_i represents the currrent
    %estimation of x
    ref = (D - omega * L);
    
    %Repeats until stopping criterion of iterations is met
    for i = 1:max_iterations
        %Solving the iteration step
        x = mldivide(ref,((1-omega)*D+omega*U)*x_init)+omega*mldivide(ref,b);
        
        %Ends the loop if stopping criterion of tolerance is met
        %Tolerance is the max. acceptable distance between x and x_init
        %(previous estimate of x)
        if norm(x - x_init) < tolerance
            break;
        end
        %x_init stores previous value of x
        x_init = x;
    end
end