% Assignment:
% Generates matrix A and b for any N and solves each iteratively
    % write uij as a vector and then Au=b is a system of equations
    % We find A and b for any N
% Jacobi, qn(1) and print plots(Kartik)
% Gauss-Seidel with over-relaxation parameter (Saket)
% Conjugate Gradient method (Heather)

%--------------------------------------------------------------------------
%Start of qn(1) code
%--------------------------------------------------------------------------

N = 20; % N

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

%--------------------------------------------------------------------------
% Start of Jacobi Method
%--------------------------------------------------------------------------

% these are initializations which are used for every program
atol = 1e-3; % absolute tolerance
programTimeLimit = 100; % runtime limit in seconds
x0 = zeros(N^2, 1); % initial guess
r = norm(b - A *x0, 1); % initial residual (b-A*xk)
xk = x0; % kth iteration
iterCount = 0; % number of iterations completed


figure(1);
ax = axes();
hold(ax)
xlabel(ax, 'iterations')
ylabel (ax, 'residual')
%title(ax,'residual against time/s and iteration/count')
% ax_top = axes();
% hold(ax_top)
% ax_top.Position = ax.Position;
% ax_top.YAxis.Visible = 'off';
% ax_top.XAxisLocation = 'top';
% xlabel(ax_top, 'time')

tic
% repeat iterations until r <= atol or program time limit reached
while r > atol 
    plot(ax,iterCount,r,'kx');
    time = toc;
%     plot(ax_top,time,r,'kx');
    title(time)
    if time > programTimeLimit
        break
    end
    drawnow
    
    % find the x in next iteration by given formula
    iterNext = zeros(N^2,1);
    for i = 1:N^2
        sum = 0;
        for j = 1:N^2
            if j ~= i
                sum = sum + A(i,j)* xk(j, 1);
            end
        end
        iterNext(i, 1) = (b(i,1)-sum)/A(i,i);
    end
    xk = iterNext;
    r = norm(b - A *xk, 1); % update residual
    iterCount = iterCount + 1; % increment iteration count
end 
% plot one more time for the last pt i.e. first r > atol
plot(ax,iterCount,r,'kx');
time = toc;
% plot(ax_top,time,r,'kx');
drawnow

%hold(ax_top,'off')
%hold(ax,'off')



