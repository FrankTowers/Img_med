
close all
fprintf(1,'\nStarting ARTdemo:\n\n');


function RMSE = calcularRMSE(x_exact, X_rec)
% CALCULARRMSE Calcula el RMSE entre la imagen exacta y la reconstrucción.
%
%   RMSE = calcularRMSE(x_exact, X_rec)
%
% Entrada:
%   x_exact   - Imagen exacta, un vector con los valores de la imagen original.
%   X_rec     - Imagen reconstruida, un vector con los valores de la imagen reconstruida.
%
% Salida:
%   RMSE      - El valor del RMSE entre la imagen exacta y la imagen reconstruida.

    % Asegurarse de que las entradas sean vectores columna
    x_exact = x_exact(:);
    X_rec = X_rec(:);

    % Calcular el RMSE
    RMSE = sqrt(mean((x_exact - X_rec).^2));
end

function [errEM] = calcErrEM(theta, p, eta,k)
    N=50;
    [A, b_ex, x_ex] = paralleltomo(N,theta,p)

    delta= eta * norm(b_ex);
     % Agregar ruido a la derecha del sistema
    randn('state',0)
    e = randn(size(b_ex));
    e = delta * e / norm(e);
    b = b_ex + e;

    % Número de iteraciones
    %k = 10;

    X_EM = em(A,b,k)
    errEM = sqrt(mean((X_EM - x_ex).^2));
end




function [errKacz,errSart, errSymk, errRand] = calcularErroresReconstruccion(theta, p, eta)
% CALCULARERRORESRECONSTRUCCION Calcula el error de reconstrucción (RMSE)
% para los métodos Kaczmarz, EM, Symmetric Kaczmarz y Randomized Kaczmarz.
%
%   [errKacz, errEm, errSymk, errRand] = calcularErroresReconstruccion(theta, p, eta)
%
% Entrada:
%   theta  - Vector de ángulos a utilizar en la tomografía.
%   p      - Número de rayos paralelos por ángulo.
%   eta    - Nivel de ruido relativo que se agregará a los datos.
%
% Salida:
%   errKacz - Error de reconstrucción (RMSE) usando Kaczmarz.
%   errEm   - Error de reconstrucción (RMSE) usando EM.
%   errSymk - Error de reconstrucción (RMSE) usando Symmetric Kaczmarz.
%   errRand - Error de reconstrucción (RMSE) usando Randomized Kaczmarz.

    % Parámetros del problema de prueba
    N = 50; % Discretización de la imagen

    % Crear el problema de tomografía
    [A, b_ex, x_ex] = paralleltomo(N, theta, p);

    % Nivel de ruido
    delta = eta * norm(b_ex);

    % Agregar ruido a la derecha del sistema
    randn('state',0)
    e = randn(size(b_ex));
    e = delta * e / norm(e);
    b = b_ex + e;

    % Número de iteraciones
    k = 10;

    % Realizar las iteraciones de Kaczmarz
    Xkacz = kaczmarz(A, b, k);

    % Realizar las iteraciones de EM
    %XEm = ex_max(A, b, k);

    % Realizar las iteraciones de Symmetric Kaczmarz
    Xsymk = symkaczmarz(A, b, k);

    % Realizar las iteraciones de Randomized Kaczmarz
    Xrand = randkaczmarz(A, b, k);

    % Realizar Sart
    Xsart = sart(A, b, k);

    % Calcular el error (RMSE) para cada método
    errKacz = sqrt(mean((Xkacz - x_ex).^2));
    errSart = sqrt(mean((Xsart - x_ex).^2));
    errSymk = sqrt(mean((Xsymk - x_ex).^2));
    errRand = sqrt(mean((Xrand - x_ex).^2));
end

##% Parámetros para el barrido
##theta = 0:5:179;    % Ángulos a utilizar
##eta = 0.05;         % Nivel de ruido
##
##% Número de puntos para el barrido en p (de 10 a 500)
##p_values = linspace(10, 500, 20);
##
##% Inicializar los vectores para almacenar los errores
##errKacz_vals = zeros(size(p_values));
##errSart_vals = zeros(size(p_values));
##errSymk_vals = zeros(size(p_values));
##errRand_vals = zeros(size(p_values));
##
##% Realizar el barrido en p
##for i = 1:length(p_values)
##    p = round(p_values(i));  % Redondear el valor de p
##    [errKacz,errSart, errSymk, errRand] = calcularErroresReconstruccion(theta, p, eta);
##
##    % Guardar los errores en los vectores correspondientes
##    errKacz_vals(i) = errKacz;
##    errSart_vals(i) = errSart;
##    errSymk_vals(i) = errSymk;
##    errRand_vals(i) = errRand;
##end
##
##% Graficar la evolución de los errores
##figure;
##
##% Graficar el error de Kaczmarz
##plot(p_values, errKacz_vals, '-o', 'DisplayName', 'Kaczmarz');
##hold on;
##
##% Graficar el error de sart
##plot(p_values, errSart_vals, '-s', 'DisplayName', 'Sart');
##
##% Graficar el error de Symmetric Kaczmarz
##plot(p_values, errSymk_vals, '-^', 'DisplayName', 'Symmetric Kaczmarz');
##
##% Graficar el error de Randomized Kaczmarz
##plot(p_values, errRand_vals, '-d', 'DisplayName', 'Randomized Kaczmarz');
##
##% Añadir etiquetas y leyenda
##xlabel('Número de Rayos Paralelos (p)');
##ylabel('Error RMSE');
##title('Evolución de los Errores con el Número de Rayos Paralelos');
##legend show;
##
##% Mejorar la presentación de la gráfica
##grid on;
##hold off;

##%Parametros barrido en eta...
##theta = 0:5:179;
##p = 75;
##
##eta_values = linspace(0.05,1,20)
##
##
####% Inicializar los vectores para almacenar los errores
##errKacz_vals = zeros(size(p_values));
##errSart_vals = zeros(size(p_values));
##errSymk_vals = zeros(size(p_values));
##errRand_vals = zeros(size(p_values));
##
##
##
##% Realizar el barrido en eta
##for i = 1:length(eta_values)  % Asegúrate de que el tamaño del bucle es el correcto
##    eta_actual = eta_values(i);  % Cambiar el nombre de la variable
##    [errKacz, errSart, errSymk, errRand] = calcularErroresReconstruccion(theta, p, eta_actual);
##
##    % Guardar los errores en los vectores correspondientes
##    errKacz_vals(i) = errKacz;
##    errSart_vals(i) = errSart;
##    errSymk_vals(i) = errSymk;
##    errRand_vals(i) = errRand;
##end
##
##% Graficar la evolución de los errores
##figure;
##
##% Graficar el error de Kaczmarz
##plot(eta_values, errKacz_vals, '-o', 'DisplayName', 'Kaczmarz');
##hold on;
##
##% Graficar el error de sart
##plot(eta_values, errSart_vals, '-s', 'DisplayName', 'Sart');
##
##% Graficar el error de Symmetric Kaczmarz
##plot(eta_values, errSymk_vals, '-^', 'DisplayName', 'Symmetric Kaczmarz');
##
##% Graficar el error de Randomized Kaczmarz
##plot(eta_values, errRand_vals, '-d', 'DisplayName', 'Randomized Kaczmarz');
##
##% Añadir etiquetas y leyenda
##xlabel('eta');
##ylabel('Error RMSE');
##title('Evolución de los Errores con eta');
##legend show;
##
##% Mejorar la presentación de la gráfica
##grid on;
##hold off;

##
##%Parametros barrido en eta...
##
##p = 75;
##eta = 0.05;
##
##theta_values = linspace(1,20,20)
##
##% Número de iteraciones
##num_iter = 20;
##
##% Rango de valores de theta
##theta_min = 1;
##theta_max = 179;
##
####% Inicializar los vectores para almacenar los errores
##errKacz_vals = zeros(size(p_values));
##errSart_vals = zeros(size(p_values));
##errSymk_vals = zeros(size(p_values));
##errRand_vals = zeros(size(p_values));
##
##
##
##% Realizar el barrido en eta
##for i = 1:20  % Asegúrate de que el tamaño del bucle es el correcto
##    theta = linspace(theta_min,theta_max,i)
##    [errKacz, errSart, errSymk, errRand] = calcularErroresReconstruccion(theta, p, eta);
##
##    % Guardar los errores en los vectores correspondientes
##    errKacz_vals(i) = errKacz;
##    errSart_vals(i) = errSart;
##    errSymk_vals(i) = errSymk;
##    errRand_vals(i) = errRand;
##end
##
##% Graficar la evolución de los errores
##figure;
##
##% Graficar el error de Kaczmarz
##plot(theta_values, errKacz_vals, '-o', 'DisplayName', 'Kaczmarz');
##hold on;
##
##% Graficar el error de sart
##plot(theta_values, errSart_vals, '-s', 'DisplayName', 'Sart');
##
##% Graficar el error de Symmetric Kaczmarz
##plot(theta_values, errSymk_vals, '-^', 'DisplayName', 'Symmetric Kaczmarz');
##
##% Graficar el error de Randomized Kaczmarz
##plot(theta_values, errRand_vals, '-d', 'DisplayName', 'Randomized Kaczmarz');
##
##% Añadir etiquetas y leyenda
##xlabel('Número de angulos');
##ylabel('Error RMSE');
##title('Evolución de los Errores con el número de angulos');
##legend show;
##
##% Mejorar la presentación de la gráfica
##grid on;

#### Iteraciones con EM:
##
%Parametros barrido en eta...
iteraciones = 20
p = 75;
eta = 0.05;
theta = 0:5:179;
lala = 1:iteraciones;


##% Inicializar los vectores para almacenar los errores
errEM_vals = zeros(1,iteraciones);



% Realizar el barrido en eta
for i = 1:iteraciones  % Asegúrate de que el tamaño del bucle es el correcto
    [errEM] = calcErrEM(theta, p, eta,i);

    % Guardar los errores en los vectores correspondientes
    errEM_vals(i) = errEM;
end

% Graficar la evolución de los errores
figure;

% Graficar el error de EM
plot(lala, errEM_vals, '-o', 'DisplayName', 'Expectation maximization');
hold on;


% Añadir etiquetas y leyenda
xlabel('Número de iteraciones');
ylabel('Error RMSE');
title('Evolución de los Errores con el número de iteraciones');
legend show;

% Mejorar la presentación de la gráfica
grid on;
pause;




