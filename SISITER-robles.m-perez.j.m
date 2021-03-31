
%Programa final: Jacobi, SOR y Gradiente conjugado.

%Entradas
n = 4;
m = 3;
a = 10;
b = 2;
epsilon = 0.01;
x = zeros(n*m,1);

%Vector B Mx=B
B = ones(n*m,1);

%omega para SOR
w = 1;

JACOBI_ProyectoNumerico(n,m,a,b,epsilon,B);
SOR_Numerico(n,m,a,b,epsilon,B,w);
Gradiente_Numerico(n,m,a,b,epsilon,B);
Gradiente_Num(n,m,a,b,B);


function [] = JACOBI_ProyectoNumerico(n,m,a,b,epsilon,B)

%Guardamos el valor del xk
xk = zeros(n*m,1);
%Guardamos el valor de xk+1
xk_1 = zeros(n*m,1);
%Residuo
r = zeros(n*m,1);

para = false;


%Jacobi
k = 1;
while (para == false )
    for i = 1: n*m

        %El caso en que i es menor a n
       if(i <= n)
           if (i == 1)
               xk_1(i) = 1/a * ( -(b*xk(i+1) - xk(i+n)) + B(i)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               xk_1(i) = 1/a * ( -(b*xk(i-1) + b*xk(i+1) - xk(i+n)) + B(i));
           else 
               xk_1(i) = 1/a * ( -(b*xk(i-1) - xk(i+n)) + B(i)); 
           end

       %El caso en que i es mayor a nm-n, es decir el último caso
       elseif(i > m*n-n)
           if (mod(i,n) == 1)
               xk_1(i) = 1/a * ( -(b*xk(i+1) - xk(i-n)) + B(i)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               xk_1(i) = 1/a * ( -(b*xk(i-1) + b*xk(i+1) - xk(i-n)) + B(i));
           else
               xk_1(i) = 1/a * ( -(b*xk(i-1) - xk(i-n)) + B(i)); 

           end

       %El caso en que i es mayor a nm-n
       else 
           if (mod(i,n) == 1)
               xk_1(i) = 1/a * ( -( -xk(i-n) + b*xk(i+1) - xk(i+n)) + B(i)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               xk_1(i) = 1/a * ( - (-xk(i-n) + b*xk(i-1) + b*xk(i+1) - xk(i+n)) + B(i)); 
           else
               xk_1(i) = 1/a * ( -( -xk(i-n) + b*xk(i-1) - xk(i+n)) + B(i)); 

           end
       end
    end 
    k = k+1;
    
    xk = xk_1;
    
    %Verifica si para o sigue calculando 
    %Se debe calcular el residuo. Si es resiuo es menor a epsilon, debe
    %parar la iteración.
    
    for i = 1: n*m

        %El caso en que i es menor a n
       if(i <= n)
           if (i == 1)
               r(i) = B(i) - (a*xk(i) + b*xk(i+1) - xk(i+n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (b*xk(i-1) + a*xk(i) + b*xk(i+1) - xk(i+n));
           else 
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) - xk(i+n)); 
           end

       %El caso en que i es mayor a nm-n, es decir el último caso
       elseif(i > m*n-n)
           if (mod(i,n) == 1)
               r(i) =  B(i) - (a*xk(i) + b*xk(i+1) - xk(i-n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) + b*xk(i+1) - xk(i-n));
           else
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) - xk(i-n)); 

           end

       %El caso en que i es mayor a nm-n
       else 
           if (mod(i,n) == 1)
               r(i) =  B(i) - (-xk(i-n) +  a*xk(i) + b*xk(i+1) - xk(i+n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (-xk(i-n) + b*xk(i-1) +  a*xk(i) + b*xk(i+1) - xk(i+n)); 
           else
               r(i) =  B(i) - (-xk(i-n) + b*xk(i-1) +  a*xk(i) - xk(i+n)); 

           end
       end
    end 
    
    if (norm(r) <= epsilon )
        para = true;
    end
    
    
end

   fprintf("\nResultado de JACOBI:");
   fprintf("\nNúmero de iteraciones necesarias: %d", k);
   fprintf("\nNorma del residuo final: %d", norm(r));
   fprintf("\nVector resultante xk:");
   for j = 1:n*m
       fprintf("\n %d",xk(j));
   end


end

%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%----------------------------- AQUI EMPIEZA SOR
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------


function [] = SOR_Numerico(n,m,a,b,epsilon,B,w)

%Guardamos el valor del xk
xk = zeros(n*m,1);
%Guardamos el valor de xk+1
xk_1 = zeros(n*m,1);
%Residuo
r = zeros(n*m,1);

para = false;


%Jacobi
k = 1;
while (para == false)
    
    for i = 1: n*m

        %El caso en que i es menor a n
       if(i <= n)
           if (i == 1)
               xk_1(i) = 1/a * ( w*B(i) + a*xk(i) - w*(a*xk(i) + b*xk(i+1) - xk(i+n)) ); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               xk_1(i) = 1/a * ( w*B(i) - w*(b*xk_1(i-1)) + a*xk(i) - w*(a*xk(i) + b*xk(i+1) - xk(i+n)) );
           else 
               xk_1(i) = 1/a * ( w*B(i) - w *(b*xk_1(i-1)) + a*xk(i) - w* ( a*xk(i)- xk(i+n) ) ); 
           end

       %El caso en que i es mayor a nm-n, es decir el último caso
       elseif(i > m*n-n)
           if (mod(i,n) == 1)
               xk_1(i) = 1/a * ( w*B(i) - w*(-xk_1(i-n)) + a*xk(i) - w*(a*xk(i) + b*xk(i+1) ) );
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               xk_1(i) = 1/a * ( w*B(i) - w*(-xk_1(i-n) + b*xk_1(i-1)) + a*xk(i) - w*(a*xk(i) + b*xk(i+1) ) );
           else
               xk_1(i) = 1/a * ( w*B(i) - w *(-xk_1(i-n) + b*xk_1(i-1)) + a*xk(i) - w* (a*xk(i)) ); 

           end

       %El caso en que i es mayor a nm-n
       else 
           if (mod(i,n) == 1)
               xk_1(i) = 1/a * ( w*B(i) - w*(-xk_1(i-n)) + a*xk(i) - w*(a*xk(i) + b*xk(i+1) - xk(i+n)) );
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               xk_1(i) = 1/a * ( w*B(i) - w*(-xk_1(i-n) + b*xk_1(i-1)) + a*xk(i) - w*(a*xk(i) + b*xk(i+1) -xk(i+n) ) ); 
           else
               xk_1(i) = 1/a * ( w*B(i) - w*(-xk_1(i-n) + b*xk_1(i-1)) + a*xk(i) - w*(a*xk(i) -xk(i+n) ) ); 

           end
       end
    end 
    k = k+1;
    
    xk = xk_1;
    
    %Verifica si para o sigue calculando 
    %Se debe calcular el residuo. Si es resiuo es menor a epsilon, debe
    %parar la iteración.
    
     for i = 1: n*m

        %El caso en que i es menor a n
       if(i <= n)
           if (i == 1)
               r(i) = B(i) - (a*xk(i) + b*xk(i+1) - xk(i+n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (b*xk(i-1) + a*xk(i) + b*xk(i+1) - xk(i+n));
           else 
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) - xk(i+n)); 
           end

       %El caso en que i es mayor a nm-n, es decir el último caso
       elseif(i > m*n-n)
           if (mod(i,n) == 1)
               r(i) =  B(i) - (a*xk(i) + b*xk(i+1) - xk(i-n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) + b*xk(i+1) - xk(i-n));
           else
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) - xk(i-n)); 

           end

       %El caso en que i es mayor a nm-n
       else 
           if (mod(i,n) == 1)
               r(i) = B(i) - (-xk(i-n) +  a*xk(i) + b*xk(i+1) - xk(i+n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (-xk(i-n) + b*xk(i-1) +  a*xk(i) + b*xk(i+1) - xk(i+n)); 
           else
               r(i) =  B(i) - (-xk(i-n) + b*xk(i-1) +  a*xk(i) - xk(i+n)); 

           end
       end
    end 
    
    if (norm(r) <= epsilon )
        para = true;
    end
    
end

   fprintf("\n\nResultado de SOR:");
   fprintf("\nNúmero de iteraciones necesarias: %d", k);
   fprintf("\nNorma del residuo final: %d", norm(r));
   fprintf("\nVector resultante xk:");
   for j = 1:n*m
       fprintf("\n %d",xk(j));
   end

end


%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%----------------------------- AQUI EMPIEZA GRADIENTE CONJUGADO
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------

function [] = Gradiente_Numerico(n,m,a,b,epsilon,B)


%Guardamos el valor del xk
xk = zeros(n*m,1);

%Residuo
rk = B;

%Los p
pk = B;

%Jacobi
k = 0;

para = false;

while (para == false)
   multPorPk = multiplicacionPorA(pk,n,m,a,b); %Bien
   alfak = (norm(rk)^2)/dot(multPorPk,pk); %Bien
   xk_1 = xk + alfak*pk; %Bien
   rk_1 = rk - alfak*multPorPk; %Bien
   betak_1 = -(norm(rk_1)^2)/(norm(rk)^2);
   pk_1 = rk_1+(betak_1*pk);
   
   rk = rk_1;
   pk = pk_1;
   xk = xk_1;
   
   k = k+1;
   
   
    %Verifica si para o sigue calculando 
    %Se debe calcular el residuo. Si es resiuo es menor a epsilon, debe
    %parar la iteración.
    
     for i = 1: n*m

        %El caso en que i es menor a n
       if(i <= n)
           if (i == 1)
               r(i) = B(i) - (a*xk(i) + b*xk(i+1) - xk(i+n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (b*xk(i-1) + a*xk(i) + b*xk(i+1) - xk(i+n));
           else 
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) - xk(i+n)); 
           end

       %El caso en que i es mayor a nm-n, es decir el último caso
       elseif(i > m*n-n)
           if (mod(i,n) == 1)
               r(i) =  B(i) - (a*xk(i) + b*xk(i+1) - xk(i-n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) + b*xk(i+1) - xk(i-n));
           else
               r(i) = B(i) - (b*xk(i-1) +  a*xk(i) - xk(i-n)); 

           end

       %El caso en que i es mayor a nm-n
       else 
           if (mod(i,n) == 1)
               r(i) = B(i) - (-xk(i-n) +  a*xk(i) + b*xk(i+1) - xk(i+n)); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               r(i) = B(i) - (-xk(i-n) + b*xk(i-1) +  a*xk(i) + b*xk(i+1) - xk(i+n)); 
           else
               r(i) =  B(i) - (-xk(i-n) + b*xk(i-1) +  a*xk(i) - xk(i+n)); 

           end
       end
    end 
    
    if (norm(r) <= epsilon )
        para = true;
    end
   
   
end

   fprintf("\n\nResultado de GRADIENTE CONJUGADO CON EPSILON:");
   fprintf("\nNúmero de iteraciones necesarias: %d", k);
   fprintf("\nNorma del residuo final: %d", norm(rk));
   fprintf("\nVector resultante xk:");
   for j = 1:n*m
       fprintf("\n %d",xk(j));
   end
   


function [vector] = multiplicacionPorA(x,n,m,a,b)
    vector = zeros(n*m,1);
    
    for i = 1: n*m
        %El caso en que i es menor a n
       if(i <= n)
           if (i == 1)
               vector(i) = a*x(i) + b*x(i+1) - x(i+n); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               vector(i) =  b*x(i-1) + a*x(i) + b*x(i+1) - x(i+n);
           else 
               vector(i) = b*x(i-1) +  a*x(i) - x(i+n); 
           end

       %El caso en que i es mayor a nm-n, es decir el último caso
       elseif(i > m*n-n)
           if (mod(i,n) == 1)
               vector(i) =  a*x(i) + b*x(i+1) - x(i-n); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               vector(i) = b*x(i-1) +  a*x(i) + b*x(i+1) - x(i-n);
           else
               vector(i) = b*x(i-1) +  a*x(i) - x(i-n); 

           end

       %El caso en que i es mayor a nm-n
       else 
           if (mod(i,n) == 1)
               vector(i) =  -x(i-n) +  a*x(i) + b*x(i+1) - x(i+n); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               vector(i) = -x(i-n) + b*x(i-1) +  a*x(i) + b*x(i+1) - x(i+n); 
           else
               vector(i) =  -x(i-n) + b*x(i-1) +  a*x(i) - x(i+n); 

           end
       end
       
    end 
      
end 
end

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------    GRADIENTE 2: SIN EL EPSILON
%-------------------------------------------------------------------------------

function [] = Gradiente_Num(n,m,a,b,B)


%Guardamos el valor del xk
xk = zeros(n*m,1);

%Residuo
rk = B;

%Los p
pk = B;

%Jacobi
k = 0;

para = false;

while (k <= n*m)
   multPorPk = multiplicacionPorA(pk,n,m,a,b); %Bien
   alfak = (norm(rk)^2)/dot(multPorPk,pk); %Bien
   xk_1 = xk + alfak*pk; %Bien
   rk_1 = rk - alfak*multPorPk; %Bien
   betak_1 = -(norm(rk_1)^2)/(norm(rk)^2);
   pk_1 = rk_1+(betak_1*pk);
   
   rk = rk_1;
   pk = pk_1;
   xk = xk_1;
   
   k = k+1;
   
   
end

   fprintf("\n\nResultado de GRADIENTE CONJUGADO SIN EL EPSILON:");
   fprintf("\nNúmero de iteraciones realizadas: %d", k);
   fprintf("\nNorma del residuo final: %d", norm(rk));
   fprintf("\nVector resultante xk:");
   for j = 1:n*m
       fprintf("\n %d",xk(j));
   end
   


function [vector] = multiplicacionPorA(x,n,m,a,b)
    vector = zeros(n*m,1);
    
    for i = 1: n*m
        %El caso en que i es menor a n
       if(i <= n)
           if (i == 1)
               vector(i) = a*x(i) + b*x(i+1) - x(i+n); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               vector(i) =  b*x(i-1) + a*x(i) + b*x(i+1) - x(i+n);
           else 
               vector(i) = b*x(i-1) +  a*x(i) - x(i+n); 
           end

       %El caso en que i es mayor a nm-n, es decir el último caso
       elseif(i > m*n-n)
           if (mod(i,n) == 1)
               vector(i) =  a*x(i) + b*x(i+1) - x(i-n); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               vector(i) = b*x(i-1) +  a*x(i) + b*x(i+1) - x(i-n);
           else
               vector(i) = b*x(i-1) +  a*x(i) - x(i-n); 

           end

       %El caso en que i es menor a nm-n
       else 
           if (mod(i,n) == 1)
               vector(i) =  -x(i-n) +  a*x(i) + b*x(i+1) - x(i+n); 
           elseif (mod(i,n) ~= 1 && mod(i,n) ~= 0)
               vector(i) = -x(i-n) + b*x(i-1) +  a*x(i) + b*x(i+1) - x(i+n); 
           else
               vector(i) =  -x(i-n) + b*x(i-1) +  a*x(i) - x(i+n); 

           end
       end
       
    end 
      
end 
end


