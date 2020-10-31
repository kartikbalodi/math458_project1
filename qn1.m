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
        % A= I B I where B = 1 4 1 , I is I3x3 identity matrix and
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

atol = 1e-3; % absolute tolerance
programTimeLimit = 50; % runtime limit in seconds
x0 = zeros(N^2, 1); % initial guess
