clear all
close all
clc
%Resistencias y Voltaje de sat
Rc=     64e3;
Rf=     1e6;
Ri=     1e3;
Rix=    10e3;
Vsat=   6.4;
%Calculos para k,alpha,s,h,Isat
k=      (Rix*Vsat)/(Rc)
alpha=  (Ri*Vsat)/(Rf)
s=      k/alpha
Isat=   Vsat/Rc
h=      2/(1 + (Ri/Rf))
n=      h/k
%Condiciones para escalon de en medio
x=  alpha-(2*h):0.0001:(2*h)-alpha;
Fx= (((x>=(-alpha-h))&(x<=(alpha-h))).*((((k/alpha)*(x+(h)))-h)))+(((x>=(h-alpha))&(x<=(h+alpha))).*((h+((k/alpha)*(x-(h))))))+(((x>=h+alpha)&(x<=(2*h)-alpha)).*(k + n*k))+(((x>=alpha)&(x<=h-alpha)).*(k))+(((x>=alpha-(2*h))&(x<-h-alpha)).*(-k-n*k))+(((x>=alpha-h)&(x<=-alpha)).*(-k))+(((x>-alpha)&(x<alpha)).*(x./alpha));
figure
plot(x,Fx)
%Ecuaciones de estado ejemplo 1.
%Parametros
a=0.7;
b=a;
c=a;
d=a;
%Numero de soluciones
n=3000000;
%Separacion entre puntos
h= 0.001;
%Estados inicilaes
X(1,1)= 0.1;
Y(1,1)= 0;
Z(1,1)= 0;
%Soluciones2rollos
% % for f=1:n
% %     X(1,f+1)= X(1,f) + h*(Y(1,f));
% %     Y(1,f+1)= Y(1,f) + h*(Z(1,f));
% %     if X(1,f) >= alpha
% %         Z(1,f+1)= Z(1,f) + h*((-a*X(1,f))-(b*Y(1,f))-(c*Z(1,f))+(d*k));
% %     elseif X(1,f) <= -alpha
% %         Z(1,f+1)= Z(1,f) + h*((-a*X(1,f))-(b*Y(1,f))-(c*Z(1,f))-(d*k));
% %     else
% %         Z(1,f+1)= Z(1,f) + h*((X(1,f)*(-a+(d*(k/alpha))))-(b*Y(1,f))-(c*Z(1,f)));
% %     end
% % end
% % figure
% % plot(X,Y)
% % title('Gráfica X,Y')
% % xlabel('Xn')
% % ylabel('Yn')
% % grid
%Soluciones4rollos
E2= 2*(k^2)/(alpha) - 2*k;
E1= E2*(-1);
for f=1:n
    X(1,f+1)= X(1,f) + h*(Y(1,f));
    Y(1,f+1)= Y(1,f) + h*(Z(1,f));
    if X(1,f) >= (2*k + alpha)
        Z(1,f+1)= Z(1,f) + h*(-a*(X(1,f))-b*(Y(1,f))-c*(Z(1,f))+d*(3*k));
    elseif ((2*k - alpha) < X(1,f))&(X(1,f) < 2*k + alpha)
        Z(1,f+1)= Z(1,f) + h*(-a*(X(1,f))-b*(Y(1,f))-c*(Z(1,f))+d*(X(1,f)*(k/alpha) + E1));
    elseif (alpha < X(1,f))&(X(1,f) < 2*k - alpha)
        Z(1,f+1)= Z(1,f) + h*(-a*(X(1,f))-b*(Y(1,f))-c*(Z(1,f))+d*(k));
    elseif (-alpha < X(1,f))&(X(1,f) < alpha)
        Z(1,f+1)= Z(1,f) + h*(-a*(X(1,f))-b*(Y(1,f))-c*(Z(1,f))+d*(X(1,f)*(k/alpha)));
    elseif ((alpha - 2*k) < X(1,f))&(X(1,f) <= -alpha)
        Z(1,f+1)= Z(1,f) + h*(-a*(X(1,f))-b*(Y(1,f))-c*(Z(1,f))+d*(-k));
    elseif ((-alpha - 2*k )< X(1,f))&(X(1,f) < (alpha - 2*k))
        Z(1,f+1)= Z(1,f) + h*(-a*(X(1,f))-b*(Y(1,f))-c*(Z(1,f))+d*(X(1,f)*(k/alpha) + E2));
    else
        Z(1,f+1)= Z(1,f) + h*(-a*(X(1,f))-b*(Y(1,f))-c*(Z(1,f))+d*(-3*k));
    end
end
figure
plot(X,Y)
title('Gráfica X,Y')
xlabel('Xn')
ylabel('Yn')
grid
