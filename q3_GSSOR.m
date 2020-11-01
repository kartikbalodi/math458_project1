%------------------------------------------
%Generating matrix A and b
%------------------------------------------
N = 3; % N

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
        %    0 I B           0 1 4   0 is 3x3 zero matrix
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


%Range of possible omega values by 0.01 step
x_init = zeros(length(A));
omega = 0:0.01:2.0;
iter_nums = zeros(length(omega));



%Stopping Conditions: a tolerance limit and max number of iterations
%before quitting
tolerance = 1e-3;
max_iterations = 10000;
  
%LDU Decomoposition of A
%Isolates lower diagonal (below the main diagonal)
L = -tril(A, -1);
%Isolates main diagonal
D = diag(diag(A));
%Isolates upper diagonal (above the main diagonal)
U = -triu(A, 1);

%Creates initial solution array of zeros
x = zeros(N^2);

%Stored in reference since it is used repeatedly
%A in the equation Ax_i = b where x_i represents the currrent
%estimation of x


num_iterations = 0;
%Repeats until stopping criterion of iterations is met
for j = 1:length(omega)
    num_iterations = 0;
    for i = 1:max_iterations
        ref = (D - omega(j) * L);
        num_iterations = num_iterations + 1;
        %Solving the iteration step
        x = mldivide(ref,((1-omega(j))*D+omega(j)*U)*x_init)+omega(j)*mldivide(ref,b);

        %Ends the loop if stopping criterion of tolerance is met
        %Tolerance is the max. acceptable distance between x and x_init
        %(previous estimate of x)
        if norm(x - x_init) < tolerance
            break;
        end
        %x_init stores previous value of x
        x_init = x;
    end
    iter_nums(j) = num_iterations;
end

    
plot(omega, iter_nums)
xlabel("Values of Omega")
ylabel("Number of Iterations Required to Solve")
