N = 200;

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
omega = 1.930;
atol = 1e-3; % absolute tolerance
programTimeLimit = 300; % runtime limit in seconds
x_init = zeros(N^2, 1);
r = norm(b - A * x_init);
num_iter = 0;

figure(1);
ax = axes();
hold(ax)
xlabel(ax, 'iterations')
ylabel (ax, 'residual')
%title(ax,'residual against time/s and iteration/count')
ax_top = axes();
hold(ax_top)
ax_top.Position = ax.Position;
ax_top.YAxis.Visible = 'off';
ax_top.XAxisLocation = 'top';
xlabel(ax_top, 'time')


%Stopping Conditions: a tolerance limit and max number of iterations
%before quitting
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
tic
for i = 1:max_iterations
    plot(ax,num_iter,r,'kx');
    time = toc;
    plot(ax_top,time,r,'kx');
    if time > programTimeLimit
        disp('Program time limit '+programTimeLimit+' reached');
        break
    end
    drawnow
	
	%Solving the iteration step
	x = mldivide(ref,((1-omega)*D+omega*U)*x_init)+omega*mldivide(ref,b);
        
	%Ends the loop if stopping criterion of tolerance is met
	%Tolerance is the max. acceptable distance between x and x_init
	%(previous estimate of x)
	if r < atol
        break;
    end
	%x_init stores previous value of x
	x_init = x;
    r = norm(b - A*x_init);
    num_iter = num_iter + 1;
end
plot(ax,num_iter,r,'kx');
time = toc;
plot(ax_top,time,r,'kx');
drawnow

disp('N = ');
disp(N);
disp('time = ');
disp(time);
disp('iterations = ');
disp(num_iter);


