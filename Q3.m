% Gauss-Seidel Method for linear equations

% check if method is convergant

A = [20 12 -3 1;
    2 -15 5 -5;
    1 3 15 -8;
    1 1 -2 -10;];
b = [7; 12; 2; 8];
x = ones(length(b),1);
TOL = 10^-10;
N = 2^10;

disp("Is the orignal matrix convergent?")
is_convergant(A)
disp("Solution for X from gauss seidel method")
gauss_seidel(A, b, x, TOL, N)
disp("Solution for X from jacobi method")
jacobi(A, b, x, TOL, N)

[A b] = setAandBforCase2A(10);
x = ones(length(b),1);
disp("Is the matrix from 2A convergent?")
is_convergant(A)
disp("Solution for X from gauss seidel method")
gauss_seidel(A, b, x, TOL, N)
disp("Solution for X from jacobi method")
jacobi(A, b, x, TOL, N)

[A b] = setAandBforCase2B(10);
% does not work
x = ones(length(b),1);
disp("Is the matrix from 2B convergent?")
is_convergant(A)
disp("Solution for X from gauss seidel method")
gauss_seidel(A, b, x, TOL, N)
disp("Solution for X from jacobi method")
jacobi(A, b, x, TOL, N)

function [x] = jacobi(A, b, x0, tol, iter)
    norms = [];
    x = x0;
    
    D = diag(diag(A));
    invD = inv(D);
    L = tril(A, -1);
    U = triu(A, 1);
    LU = L+U;
    
    for k = 1:iter
        x0 = x;
        x = -invD*LU*x0 + invD*b;
        
        norms(k) = norm(((A*x)-b));
        if norm(x-x0) < tol
            break
        end
        x0 = x;
    end
    
    xAxis = 1:k;
    figure
    plot(xAxis, norms);
    xlabel('Iteration'); ylabel('Norm of solution Error');
    title("Jacobi Method");
end



function [x] = gauss_seidel(A, b, x0, tol, iter)
    norms = [];
    x = x0;
    for k = 1:iter
        for i = 1:length(A)
            s1 = dot(A(i, 1:i-1), x(1:i-1));
            s2 = dot(A(i, i+1:end), x(i+1:end));
            x(i) = (b(i) - s1 - s2) / A(i, i);
        end
        norms(k) = norm(((A*x0)-b));
        if norm(x-x0) < tol
            break
        end
        x0 = x;
    end
    xAxis = 1:k;
    figure
    plot(xAxis, norms);
    xlabel('Iteration'); ylabel('Norm of solution Error');
    title("Gauss Seidel Method");
    
end

function [is_conv] = is_convergant(A)
    
    is_conv = false;
    % is symmetric positive-definite
    if issymmetric(A)
        is_conv = all(eig(A) > eps);
        return
    end
    
    
    % is striclty or irredivulby diagonally dominant
    %row
    row_dominant = all((2*abs(diag(A))) > sum(abs(A),2), 'all');
    % column
    col_dominant = all((2*abs(diag(A))) > sum(abs(A),1), 'all');
    
    if row_dominant || col_dominant
        is_conv = true;
    end
    
end

function [A, b] = setAandBforCase2A (n)
A = zeros(n);

% fill in a
for i = 1:n
    for j = 1:n
        if i == j
            A(i, j) = 8;
        elseif (i == j-1) || (i == j+1)
            A(i, j) = 2;
        end
    end
end

b = [4:0.7:(0.7*n+3.3)]';

end

function [A, b] = setAandBforCase2B (n)
A = zeros(n);
b = zeros(n, 1);

% fill in a
for i = 1:n
    for j = 1:n
        A(i, j) = 3/(5*(i+j+1));
    end
end

% fill in b
for i = 1:n
    if mod(i, 2) == 1
        b(i) = 0;
    else 
        b(i) = 5/(4*i);
    end
end

end

