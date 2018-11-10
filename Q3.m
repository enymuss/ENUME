% Gauss-Seidel Method for linear equations

% check if method is convergant

A = [16 3; 7 -11]
b = [11; 13]
x = [1; 1]
TOL = 2^-13;
N = 2^10;

gauss_seidel(A, b, x, TOL, N)

B = [1 2; 2 1]
is_convergant(B)

function [x] = gauss_seidel(A, b, x0, tol, iter)
    x = x0;
    for k = 1:iter
        for i = 1:length(A)
            s1 = dot(A(i, 1:i-1), x(1:i-1));
            s2 = dot(A(i, i+1:end), x(i+1:end));
            x(i) = (b(i) - s1 - s2) / A(i, i);
        end
        if norm(x-x0) < tol
            break
        end
        x0 = x;
    end
end

function [is_conv] = is_convergant(A)
    is_conv = false;
    % is symmetric positive-definite
    if issymmetric(A)
        is_conv = all(eig(A) > eps)
        return
    end
    
    
    % is striclty or irredivulby diagonally dominant
    %row
    row_dominant = all((2*abs(diag(A))) > sum(abs(A),2), 'all')
    % column
    col_dominant = all((2*abs(diag(A))) > sum(abs(A),1), 'all')
    
    if row_dominant || col_dominant
        is_conv = true
    end
    
end


