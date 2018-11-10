%Gaussian elimination with partial pivoting
clc;

% Plot graphs

[iterations, norms, matlabComparison] = iterateforGraph('A', 6);
figure;
plot(iterations, norms, iterations, matlabComparison, '--');
xlabel('n'); ylabel('Norm of Residuum');
title("Case A");
legend("Progam", "Matlab")


[iterations, norms, matlabComparison] = iterateforGraph('B', 6);
figure;
plot(iterations, norms, iterations, matlabComparison, '--');
xlabel('n'); ylabel('Norm of Residuum');
title("Case B");
legend("Progam", "Matlab");

function [iterations, norms, matlabComparison] = iterateforGraph(matrixCase, numOfIterations)
iterations = zeros(numOfIterations);
norms = zeros(numOfIterations);
matlabComparison = zeros(numOfIterations);
A = [];
b = [];
for i = 1:numOfIterations
    n = 10*(2^(i-1));
    iterations(i) = n;
    
    if (matrixCase == 'A')
        [A, b] = setAandBforCaseA(n);
    else
        [A, b] = setAandBforCaseB(n);
    end
    
    x = gaussDivide(A, b);
    residuum = (A*x)-b;
    norms(i) = norm(residuum);
    
    if (n == 10)
        %special case for n = 10, make residual correction
        disp("Case " + matrixCase)
        norms(i) = norm(makeResidualCorrection(A, b));
    end
    
    matlabX = mldivide(A, b);
    residuum2 = (A*matlabX)-b;
    matlabComparison(i) = norm(residuum2);
end
end

function [residuum] = makeResidualCorrection(A, b)
disp("Case for n = 10")
x = gaussDivide(A, b);
residuum = (A*x)-b;
delX = gaussDivide(A, residuum);
newX = x - delX;
newR = (A*newX)-b;
nResiduum =  norm(residuum);
nNewR = norm(newR);

disp("Solution")
x'
disp("Solution Error")
residuum'
disp("New Solution")
newX'
disp("New Solution Error")
newR'

while (nResiduum - nNewR > eps)
    disp("Residumm is decreasing")
    residuum = newR;
    delX = gaussDivide(A, residuum);
    newX = newX - delX;
    newR = (A*newX)-b;
    nResiduum =  norm(residuum)
    nNewR = norm(newR)
    
    disp("Solution")
    newX'
    disp("Solution Error")
    newR'
end

disp("Residumm is not decreasing")

if nResiduum < nNewR;
    residuum = residuum;
else
    residuum = newR;
end

end




%matlabX = mldivide(A, b);
%solutionErrorandNorm(A, b, matlabX)'


%{
delX = mldivide(A,r)
newX = x - delX
newR = (A*newX)-b
norm(r)
norm(newR)
%}


function [x] = gaussDivide(A, b)
    n = length(b);
    
    %Step 1 - Put into Triangluar Form
    triangle = upperTriangleWithPartialPivoting([A b], n);

    %Step 2- Solve by back substitution
    x = backSubstitution(triangle, n);
end

function [newRow] = zeroFristElement (row1, row2)
    rowMultiplier = row2(1,1) / row1(1,1);
    newRow = row2' - rowMultiplier * row1';
    newRow = newRow';
end


function [newSum] = sumPreviousElements (colA, colX)
    newVec = colA .* colX;
    newSum = sum(newVec, 'all');
end

function [A, b] = setAandBforCaseA (n)
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

function [A, b] = setAandBforCaseB (n)
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

function [upper] = upperTriangleWithPartialPivoting (triangle, n)
upper = triangle;
for i = 1:n-1
    % find row of largest element in column
    [mavValue maxValueIndex] = max(abs(upper(i:end, i)));
    maxValueIndex = maxValueIndex + (i-1);
    
    % switch largest element with the top element
    tempRow = upper(maxValueIndex , :);
    upper(maxValueIndex , :) = upper(i, :);
    upper(i, :) = tempRow;
    for j = i:n-1
        row1 = upper(i, i:end);
        row2 = upper(j+1, i:end);
        upper(j+1, i:end) = zeroFristElement(row1, row2);
    end
end
end

function [x] = backSubstitution (upper, n)
x = zeros(n,1);
for i = n:-1:1
    numerator = upper(i, end) - sumPreviousElements (upper(i, i:end-1)', x(i:end));
    x(i) = numerator/upper(i, i);
end
end


