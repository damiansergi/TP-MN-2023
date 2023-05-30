x = [0.5 0.8 1.5 2.5 4];
k = [1.1 2.4 5.3 7.6 8.9];

%{
k(x) = kmax x**2/(x**2+cs) -> 
f(x) = 1/k(x) =1/kmax*(x**2+cs)/x**2

x -> u, x > 0 
x**2 = 1/u

f(u) = 1/kmax+cs/kmax*u
%}
f = 1./k;
u = 1./(x.*x);


%{
 A'A*x=A'b
 x = [cs/kmax 1/kmax]
%}
A = [u' ones(columns(u), 1)];
y = f';

params = linsolve(A'*A, A'*y);
disp("Los parametros estimados son: ")
kmax = 1/params(2)
cs = params(1)*kmax
disp("Ahora reemplazo en la funcion con los parametros reemplazados con c = 2:")
c = 2;
printf("k(2) = %f \n", kmax*c^2/(c^2+cs))
