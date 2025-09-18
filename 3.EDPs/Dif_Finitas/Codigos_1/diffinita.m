% Ejemplo de una funcion para resolver la ecuacion diferencial parcial en 1D
% Uxx = -Pi*Pi*cos(Pi*x)
% xi <= U <= xf
% U(xi) = vi y U(xf) = vf

function [A, b, x] = diffinita(n)
    xi = -1;          % Inicio de dominio
    xf = 2;           % Fin de dominio
    vi = -1;          % Valor en la frontera xi
    vf = 1;           % Valor en la frontera xf
    N = n - 2;        % Nodos interiores
    h = (xf - xi) / (n - 1);  % Incremento en la malla
    A = zeros(N, N);     % Matriz A
    b = zeros(N, 1);     % Vector b

    R = 1 / (h^2);
    P = -2 / (h^2);
    Q = 1 / (h^2);

    % Primer renglón de la matriz A y vector b
    A(1, 1) = P;
    A(1, 2) = Q;
    b(1) = LadoDerecho(xi) - vi * R;

    % Renglones intermedios de la matriz A y vector b
    for i = 2:N - 1
        A(i, i - 1) = R;
        A(i, i) = P;
        A(i, i + 1) = Q;
        b(i) = LadoDerecho(xi + h * (i - 1));
    end

    % Renglón final de la matriz A y vector b
    A(N, N - 1) = R;
    A(N, N) = P;
    b(N) = LadoDerecho(xi + h * N) - vf * Q;

    % Resuelve el sistema lineal Ax = b
    x = inv(A) * b;

    % Prepara la graficación
    xx = zeros(n, 1);
    zz = zeros(n, 1);
    for i = 1:n
        xx(i) = xi + h * (i - 1);
        zz(i) = SolucionAnalitica(xx(i));
    end

    yy = zeros(n, 1);
    yy(1) = vi;         % Condición inicial
    for i = 1:N
        yy(i + 1) = x(i); % Solución numérica en los nodos interiores
    end
    yy(n) = vf;        % Condición final

    % Graficar la solución de la Ecuación Diferencial Parcial en 1D
    figure; % Crea una nueva ventana para la gráfica
    plot(xx, yy, 'b-', 'DisplayName', 'Solución numérica'); % Solución numérica
    hold on;
    plot(xx, zz, 'r--', 'DisplayName', 'Solución analítica'); % Solución analítica
    xlabel('x'); % Etiqueta del eje x
    ylabel('U(x)'); % Etiqueta del eje y
    title('Solución de la Ecuación Diferencial Parcial en 1D'); % Título de la gráfica
    legend show; % Muestra la leyenda
    grid on; % Activa la cuadrícula
endfunction

% Lado derecho de la ecuación
function y = LadoDerecho(x)
    y = -pi^2 * cos(pi * x);
endfunction

% Solución analítica a la ecuación
function y = SolucionAnalitica(x)
    y = cos(pi * x);
endfunction

% Ejecutar la función
[A, b, x] = diffinita(30);  % Aquí se pasa el valor de n correctamente






