N = 40;

% initialize b to the value of g(xij,  yij)
h = 1/(N+1);
b = zeros(N^2,1);
for i = 1:N
    for j = 1:N
        if i/(N+1) > 1/5 && i/(N+1) < 3/5 && j/(N+1) > 1/4 && j/(N+1)<1/2
            b((j-1)*N+i,1) = -h^2;
        end
    end
end

A = zeros(N^2,N^2);
for j = 1:N
    for i = 1:N
        % First take care of the tridiagonal matrices
        c = (j-1)*N;
        A(c+i,c+i)= 4;
        if i ~= 1
            A(c+i, c+i-1) = -1;
        end
        if i ~= N
            A(c+i, c+i+1) = -1;
        end
        % Then fill in the identity matrices
        if j > 1
            A(c+i-N, c+i) = -1;
        end
        if j < N
            A(c+i+N, c+i) = -1;
        end
    end
end

%Set value of omega, can change for tests
omega = 0.01:0.01:1.99;
x = zeros(N^2, 1);

least_time = 10000000;
best_omega = 0;
time_taken = zeros(length(omega),1);
for i = i:length(omega)
    tic;
    x = SORMethod(A,b,x, omega(i));
    time = toc;
    time_taken(i) = time;
    if time < least_time
        best_omega = omega(i);
        least_time = time;
    end
end

plot(omega, time_taken)
xlabel("Values of Omega")
ylabel("Time Taken to Complete")

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
