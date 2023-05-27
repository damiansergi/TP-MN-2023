%La función en la cual se aplica el método de punto fijo es f(t)=t-(12*(1-e^(-0.04t))+5*e^(-0.04t)-10.2)

%Para encontrar un intervalo en donde se encuentre la solución basta con
%evaluar f(t) = 12*(1-e^(-0.04*(t)))+5*e^(-0.04*(t))-10.2; (donde 10.2 = 0.85*Ce)
%y encontrar dos puntos x1, x2 tal que f(x1)*f(x2)<0. Tales podrían ser por ejemplo
%x1=0 y x2=50
%Es válido que la raíz sea única en este punto ya que f'(t)>0 para todo t>0
%Por lo que la función es monótona creciente y la raíz es única


format long;
function y=g(t) %Funcion para el metodo de punto fijo
  y = (t)-(12*(1-e^(-0.04*(t)))+5*e^(-0.04*(t))-10.2);
endfunction

function y=f(t) %Funcion para el metodo de Newton-Raphson
  y = 12*(1-e^(-0.04*(t)))+5*e^(-0.04*(t))-10.2;
endfunction

function y=df(t)  %Funcion para el metodo de Newton-Raphson
  y = (0.04)*12*e^(-0.04*(t))+(-0.04)*5*e^(-0.04*(t));
endfunction

%Parametros para el método de punto fijo
x0FP = 25;  %Tomo el punto medio del intervalo previamente definido
cotaErrorFP = 10^(-6);

%Funcion para el método de punto fijo
xAntFP = x0FP;
xSigFP = g(xAntFP);
contFP = 1;
while abs(xSigFP-xAntFP)>cotaErrorFP
    xAntFP=xSigFP;
    xSigFP=g(xAntFP);
    contFP=contFP+1;
end

%Soluciones para el método de punto fijo
disp("Las soluciones con el método del punto fijo son:")
solFP=xSigFP
nIterFP=contFP
errorFP=abs(xSigFP-xAntFP)

%Parametros para el método de Newton Raphson
x0NR = 25; %Tomo el punto medio del intervalo previamente definido
cotaErrorNR = 10^(-6);

%Funcion para el método de Newton Raphson
xAntNR = 0;
xSigNR = x0NR;
contNR = 0;
while abs(xSigNR-xAntNR)>cotaErrorNR
    xAntNR=xSigNR;
    xSigNR=xAntNR - ( f(xAntNR) / df(xAntNR) );
    contNR=contNR+1;
end

%Soluciones para el método de Newton Raphson
disp("Las soluciones con el método de Newton Raphson son:")
solNR=xSigNR
nIterNR=contNR
errorNR=abs(xSigNR-xAntNR)

