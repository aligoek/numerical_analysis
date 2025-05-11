#include <stdio.h>
#include <math.h>
#define MAX 20
#define SIZE 15

typedef struct term {
    int degree;
    double coefficient;
} TERM;

double evaluatePolynomial(TERM polynomial[MAX], double x, int degree);
int getPolynomial(TERM polynomial[MAX]);
void bisectionMethod();
void falsePositionMethod();
void newtonRaphsonMethod();
void numericalDerivative();
void trapezoidalMethod();
void simpsonsMethod();
void matrixInverse();
void gaussSeidelMethod();
void gaussEliminationMethod();
void gregoryNewtonInterpolation();
double productTerm(double x, int n, int x0, int h);
int factorial(int n);

int main() {
    int continueProgram = 1;

    while (continueProgram == 1) {
        int method;
        printf("\nWhich method would you like to try?\n");
        printf("1. Bisection\n");
        printf("2. False Position\n");
        printf("3. Newton Raphson\n");
        printf("4. Matrix Inverse\n");
        printf("5. Gauss Elimination\n");
        printf("6. Gauss Seidel\n");
        printf("7. Numerical Derivative\n");
        printf("8. Simpson's Method\n");
        printf("9. Trapezoidal Method\n");
        printf("10. Gregory Newton\n");
        printf("0. Exit\n");
        printf("Your choice: ");
        scanf("%d", &method);

        switch (method) {
            case 1:
                bisectionMethod();
                break;
            case 2:
                falsePositionMethod();
                break;
            case 3:
                newtonRaphsonMethod();
                break;
            case 4:
                matrixInverse();
                break;
            case 5:
                gaussEliminationMethod();
                break;
            case 6:
                gaussSeidelMethod();
                break;
            case 7:
                numericalDerivative();
                break;
            case 8:
                simpsonsMethod();
                break;
            case 9:
                trapezoidalMethod();
                break;
            case 10:
                gregoryNewtonInterpolation();
                break;
            case 0:
                printf("Exiting program.\n");
                continueProgram = 0;
                break;
            default:
                printf("Invalid selection.\n");
        }
    }
    return 0;
}

int getPolynomial(TERM polynomial[MAX]) {
    int i, degree;
    printf("Enter the degree of the polynomial: ");
    scanf("%d", &degree);

    for (i = 0; i <= degree; i++) {
        printf("Enter the degree of the %d. term: ", i + 1);
        scanf("%d", &polynomial[i].degree);
        printf("Enter the coefficient of the %d. term: ", i + 1);
        scanf("%lf", &polynomial[i].coefficient);
    }

    return degree;
}

double evaluatePolynomial(TERM polynomial[MAX], double x, int degree) {
    double total = 0.0;
    double factor;
    int i, j;
    for (i = 0; i <= degree; i++) {
        factor = 1.0;
        for (j = 0; j < polynomial[i].degree; j++) {
            factor = factor * x;
        }
        total = total + (polynomial[i].coefficient * factor);
    }
    return total;
}

void bisectionMethod() {
    int iteration;
    double errorTolerance, a, b, c;
    TERM polynomial[MAX];
    int degree = getPolynomial(polynomial);
    iteration = 0;

    printf("Enter the error tolerance: ");
    scanf("%lf", &errorTolerance);

    do {
        printf("Enter the start value of the interval: ");
        scanf("%lf", &a);

        printf("Enter the end value of the interval: ");
        scanf("%lf", &b);
    } while (evaluatePolynomial(polynomial, a, degree) * evaluatePolynomial(polynomial, b, degree) > 0);

    do {
        double Fa = evaluatePolynomial(polynomial, a, degree);
        double Fb = evaluatePolynomial(polynomial, b, degree);

        c = (a + b) / 2;
        double Fc = evaluatePolynomial(polynomial, c, degree);

        if (Fa * Fc < 0) {
            b = c;
        } else if (Fa * Fc > 0) {
            a = c;
        }
        iteration++;

    } while ((b - a) / pow(2, iteration) > errorTolerance);

    printf("\nApproximated root: %lf\n", c);
}

void falsePositionMethod() {
    int iteration;
    double errorTolerance, a, b, c;
    TERM polynomial[MAX];
    int degree = getPolynomial(polynomial);
    iteration = 0;

    printf("Enter the error tolerance: ");
    scanf("%lf", &errorTolerance);

    do {
        printf("Enter the start value of the interval: ");
        scanf("%lf", &a);

        printf("Enter the end value of the interval: ");
        scanf("%lf", &b);
    } while (evaluatePolynomial(polynomial, a, degree) * evaluatePolynomial(polynomial, b, degree) > 0);

    do {
        double Fa = evaluatePolynomial(polynomial, a, degree);
        double Fb = evaluatePolynomial(polynomial, b, degree);

        c = (a * Fb - b * Fa) / (Fb - Fa);
        double Fc = evaluatePolynomial(polynomial, c, degree);

        if (Fa * Fc < 0) {
            b = c;
        } else if (Fa * Fc > 0) {
            a = c;
        }
        iteration++;

    } while ((b - a) / pow(2, iteration) > errorTolerance);

    printf("\nApproximated root: %lf\n", c);
}

void newtonRaphsonMethod() {
    int iteration;
    double errorTolerance, a;
    TERM polynomial[MAX];
    TERM derivative[MAX];
    double oldValue, newValue;
    int degree = getPolynomial(polynomial);
    printf("\n---Enter the derivative of the polynomial:--\n");
    int degree2 = getPolynomial(derivative);
    int exitFlag = 0;

    printf("\nEnter the error tolerance: ");
    scanf("%lf", &errorTolerance);

    printf("Enter the initial guess: ");
    scanf("%lf", &a);

    oldValue = a;

    do {
        double Fpolynomial = evaluatePolynomial(polynomial, oldValue, degree);
        double Fderivative = evaluatePolynomial(derivative, oldValue, degree2);

        newValue = oldValue - (Fpolynomial / Fderivative);

        if (fabs(newValue - oldValue) > errorTolerance) {
            oldValue = newValue;
        } else {
            exitFlag = 1;
        }

    } while (exitFlag == 0);

    printf("\nApproximated root: %lf\n", newValue);
}

void numericalDerivative() {
    int choice;
    double x, h, derivativeValue;
    TERM polynomial[MAX];
    int degree = getPolynomial(polynomial);

    printf("\nWhich method would you like to use to calculate the derivative?(1. Backward Difference 2. Forward Difference 3. Central Difference):");
    scanf("%d", &choice);
    printf("\nEnter the point at which you want to calculate the derivative:");
    scanf("%lf", &x);
    printf("Enter the value of h:");
    scanf("%lf", &h);

    switch (choice) {
        case 1:
            derivativeValue = (evaluatePolynomial(polynomial, x, degree) - evaluatePolynomial(polynomial, x - h, degree)) / h;
            break;
        case 2:
            derivativeValue = (evaluatePolynomial(polynomial, x + h, degree) - evaluatePolynomial(polynomial, x, degree)) / h;
            break;
        case 3:
            derivativeValue = (evaluatePolynomial(polynomial, x + h, degree) - evaluatePolynomial(polynomial, x - h, degree)) / (2.0 * h);
            break;
        default:
            printf("Invalid choice.\n");
            return;
    }

    printf("\nNumerically calculated derivative value:%lf\n", derivativeValue);
}

void trapezoidalMethod() {
    double integral, xi, a, b, h, ySum = 0;
    int i, n;
    TERM polynomial[MAX];
    int degree = getPolynomial(polynomial);

    printf("\nEnter the lower limit of the integral:");
    scanf("%lf", &a);

    printf("\nEnter the upper limit of the integral:");
    scanf("%lf", &b);

    printf("\nEnter the value of n:");
    scanf("%d", &n);

    h = (b - a) / n;
    xi = a;

    ySum = ySum + (evaluatePolynomial(polynomial, a, degree) + evaluatePolynomial(polynomial, b, degree)) / 2.0;
    for (i = 1; i <= n - 1; i++) {
        xi = xi + h;
        ySum = ySum + evaluatePolynomial(polynomial, xi, degree);
    }

    integral = h * ySum;

    printf("\nIntegral value calculated by Trapezoidal method:%lf\n", integral);
}

void simpsonsMethod() {
    double integral, xi, a, b, h, ySum = 0;
    int i, n;
    TERM polynomial[MAX];
    int degree = getPolynomial(polynomial);

    printf("\nEnter the lower limit of the integral:");
    scanf("%lf", &a);

    printf("\nEnter the upper limit of the integral:");
    scanf("%lf", &b);
    do {
        printf("\nEnter the value of n (must be even):");
        scanf("%d", &n);
    } while (n % 2 != 0);

    h = (b - a) / n;
    xi = a;

    ySum = ySum + (evaluatePolynomial(polynomial, a, degree) + evaluatePolynomial(polynomial, b, degree));
    for (i = 1; i <= n - 1; i++) {
        xi = xi + h;
        if (i % 2 == 1) {
            ySum = ySum + 4.0 * evaluatePolynomial(polynomial, xi, degree);
        } else {
            ySum = ySum + 2.0 * evaluatePolynomial(polynomial, xi, degree);
        }
    }

    integral = (h / 3.0) * ySum;

    printf("\nIntegral value calculated by Simpson's method:%lf\n", integral);
}

void gaussSeidelMethod() {
    double equations[SIZE][SIZE], values[SIZE], tempValues[SIZE];
    double errorTolerance, sum;
    int n, i, j, k;
    int iteration = 0;
    int errorCheck = 1;

    printf("Number of equations:");
    scanf("%d", &n);
    printf("Error tolerance:");
    scanf("%lf", &errorTolerance);

    for (i = 1; i <= n; i++) {
        printf("Enter coefficients for equation %d:", i);
        for (j = 1; j <= n; j++) {
            scanf("%lf", &equations[i][j]);
        }
        printf("\n");
    }

    printf("Enter the result column:");
    for (i = 1; i <= n; i++) {
        scanf("%lf", &equations[i][n + 1]);
    }

    printf("Enter initial values for the variables (x1, x2, ...):");
    for (i = 1; i <= n; i++) {
        scanf("%lf", &values[i]);
    }

    for (i = 1; i <= n; i++) {
        for (k = i + 1; k <= n; k++) {
            if (fabs(equations[i][i]) < fabs(equations[k][i])) {
                for (j = 1; j <= n + 1; j++) {
                    double temp = equations[i][j];
                    equations[i][j] = equations[k][j];
                    equations[k][j] = temp;
                }
            }
        }
    }

    while (errorCheck == 1) {
        errorCheck = 0;
        iteration++;

        for (i = 1; i <= n; i++) {
            sum = 0;
            tempValues[i] = values[i];

            for (j = 1; j <= n; j++) {
                if (j != i)
                    sum += equations[i][j] * values[j];
            }

            values[i] = (equations[i][n + 1] - sum) / equations[i][i];
        }

        for (i = 1; i <= n; i++) {
            if (fabs(tempValues[i] - values[i]) >= errorTolerance)
                errorCheck = 1;
        }
    }
    printf("\nConverged values after %d iterations:\n\n", iteration);

    for (i = 1; i <= n; i++) {
        printf("x%d = %lf\n", i, values[i]);
    }
}

void matrixInverse() {
    double a[SIZE][SIZE * 2];
    int i, j, k, x, n;
    int currentRow = 0;
    int dividingRow = 1;
    int column = 0;

    printf("Enter the size of the square matrix (N):");
    scanf("%d", &n);

    for (i = 1; i <= n; i++) {
        for (j = 1; j <= 2 * n; j++) {
            if (j == (i + n)) {
                a[i][j] = 1;
            } else {
                a[i][j] = 0;
            }
        }
    }

    for (i = 1; i <= n; i++) {
        printf("Enter row %d:", i);
        for (j = 1; j <= n; j++) {
            scanf("%lf", &a[i][j]);
        }
    }

    for (i = 1; i <= n * n; i++) {

        if ((i + n) % n == 1) {
            column++;
            currentRow++;
            double divisor = a[dividingRow][column];
            for (k = 1; k <= 2 * n; k++) {
                a[dividingRow][k] = a[dividingRow][k] / divisor;
            }
            dividingRow++;

        } else {
            for (k = 1; k <= n; k++) {
                if (k != dividingRow - 1) {
                    double factor = a[k][column];
                    for (x = 1; x <= 2 * n; x++) {
                        a[k][x] = a[k][x] - a[currentRow][x] * factor;
                    }
                }
            }
            i = i + (n - 2);
        }
    }
    printf("\nInverse of the matrix:\n");
    for (i = 1; i <= n; i++) {
        for (j = 1 + n; j <= 2 * n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
}

void gaussEliminationMethod() {
    double a[SIZE][SIZE + 1], x[SIZE];
    double temp;
    int i, j, k, n;

    printf("Number of equations:");
    scanf("%d", &n);

    for (i = 0; i < n; i++) {
        printf("Coefficients for row %d:", i + 1);
        for (j = 0; j < n; j++) {
            scanf("%lf", &a[i][j]);
        }
    }

    printf("Result column:");
    for (i = 0; i < n; i++) {
        scanf("%lf", &a[i][n]);
    }

    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (fabs(a[i][i]) < fabs(a[j][i])) {
                for (k = 0; k < n + 1; k++) {
                    temp = a[i][k];
                    a[i][k] = a[j][k];
                    a[j][k] = temp;
                }
            }
        }
        for (j = i + 1; j < n; j++) {
            double ratio = a[j][i] / a[i][i];
            for (k = 0; k < n + 1; k++) {
                a[j][k] = a[j][k] - ratio * a[i][k];
            }
        }
    }

    for (i = n - 1; i >= 0; i--) {
        x[i] = a[i][n];
        for (j = i + 1; j < n; j++) {
            x[i] = x[i] - a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }

    printf("\nSolutions:\n");
    for (i = 0; i < n; i++) {
        printf("x%d = %lf ", i + 1, x[i]);
    }
    printf("\n");
}

void gregoryNewtonInterpolation() {
    int n, i, j, h;
    double k, evaluatedY, targetX;
    double x[SIZE], y[SIZE][SIZE];

    printf("Enter the number of x and y values:");
    scanf("%d", &n);
    printf("Enter the x value for which you want to see the result:");
    scanf("%lf", &targetX);

    for (i = 0; i < n; i++) {
        printf("%d. x - > f(x) = ", i + 1);
        scanf("%lf %lf", &x[i], &y[i][0]);
    }

    h = x[1] - x[0];

    for (i = 1; i < n; i++) {
        for (j = 0; j < n - i; j++) {
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }
    }

    evaluatedY = y[0][0];

    for (i = 1; i < n; i++) {
        evaluatedY = evaluatedY + (productTerm(targetX, i, x[0], h) * y[0][i]) / (pow(h, i) * factorial(i));
    }

    printf("\nf(%d) = %lf\n", (int)targetX, evaluatedY);
}

int factorial(int n) {
    int i;
    int f = 1;

    for (i = 2; i <= n; i++) {
        f = f * i;
    }
    return f;
}

double productTerm(double x, int n, int x0, int h) {
    int i;
    double temp = x - x0;

    for (i = 1; i < n; i++) {
        x0 = x0 + h;
        temp = temp * (x - x0);
    }

    return temp;
}
