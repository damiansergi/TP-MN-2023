% El sistema a resolver es Ax=b dado por las siguientes matrices
A = [3, -1, 0, 0, 0;
     0, 5, 0, -1, 2;
     1, -1, 9, 0, 0;
     -2, 3, -1, 12, 1;
     0, 0, -1, 1, 15];
     
b = [-1; 0; 1; 1; 0];

% Las funciones utilizadas para resolver el sistema son 
% Eliminacion Gaussiana
function x = gaussianElimination(A, b)
    n = length(b);

    % Check determinant
    if det(A) == 0
        error("The determinant of matrix A is zero. The system may not have a unique solution.");
    end

    % Forward elimination
    for k = 1:n-1
        for i = k+1:n
            factor = A(i,k) / A(k,k);
            A(i,:) = A(i,:) - factor * A(k,:);
            b(i) = b(i) - factor * b(k);
        end
    end

    % Back substitution
    x = zeros(n, 1);
    x(n) = b(n) / A(n,n);
    for i = n-1:-1:1
        x(i) = (b(i) - A(i,i+1:n) * x(i+1:n)) / A(i,i);
    end
  
endfunction
  
% Con descomposicion LU
function x = LU(A, b)
    % LU factorization
    n = length(b);
    L = eye(n);
    U = A;

    for k = 1:n-1
        for i = k+1:n
            factor = U(i, k) / U(k, k);
            L(i, k) = factor;
            U(i, k:n) = U(i, k:n) - factor * U(k, k:n);
        end
    end

    % Forward substitution
    y = zeros(n, 1);

    for i = 1:n
        y(i) = b(i) - L(i, 1:i-1) * y(1:i-1);
    end

    % Back substitution
    x = zeros(n, 1);

    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
  
endfunction
  
% Con descomposicion QR
function x = QR(A, b)
    % QR factorization
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n);

    for j = 1:n
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
    
    % Back substitution
    y = Q' * b;
    
    n = size(R, 2);
    x = zeros(n, 1);

    for i = n:-1:1
        x(i) = (y(i) - R(i, i+1:n) * x(i+1:n)) / R(i, i);
    end
endfunction

% Con el metodo iterativo de Jacobi 
function x = Jacobi(A, b, maxIt, error)
    n = length(b);
    x = zeros(n, 1);
    xPrev = x;

    for iteration = 1:maxIt
        for i = 1:n
            sum = 0;
            for j = 1:n
                if i ~= j
                    sum = sum + A(i, j) * xPrev(j);
                end
            end
            x(i) = (b(i) - sum) / A(i, i);
        end

        if norm(x - xPrev) < error
            break;
        end

        xPrev = x;
    end
    
endfunction


%______________________________________________________________________________
% El sistema a resolver
disp("El sistema a resolver:")
[n, m] = size(A);
for i = 1:n
    equation = "";
    for j = 1:m-1
        equation = strcat(equation, num2str(A(i, j)), "x", num2str(j), " + ");
    end
    equation = strcat(equation, num2str(A(i, m)), "x", num2str(m), " = ", num2str(b(i)));
    disp(equation);
end

% 1a) Resolver por eliminacion Gaussiana
x = gaussianElimination(A, b);

disp("\nCon Eliminacion Gaussiana:")
disp("x = ");
disp(x);

% 1b) Resolver con factorizacion LU
x = LU(A, b);

disp("Con factorizacion LU:")
disp("x = ");
disp(x);

% 1c) Resolver con factorizacion QR
x = QR(A, b);

disp("Con factorizacion QR:")
disp("x = ");
disp(x);

% 1d) Resolver con el metodo de Jacobi con Error<10^-6
x = Jacobi(A, b, 1000, 1e-06);

disp("Con el metodo de Jacobi:")
disp("x = ");
disp(x);