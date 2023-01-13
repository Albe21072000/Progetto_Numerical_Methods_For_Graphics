clear ; %cancello i risultati delle esecuzioni precednti
a=[1/3,5,3]; %vettori dei termini dei coefficienti lineari
b=[3,5,1/3];
a1=a;
b1=b;
n=length(a)+1; %calcolo il numero di punti necessario
if(n~=length(b1)+1)
    error("I vettori contenenti i termini a e b dei coefficienti lineari hanno lunghezze diverse")
end
t=0.5; %parametro nel quale si vuole valutare la curva
tt=1-t;
L=a.*tt+b.*t; %ricavo i valori dei coefficienti lineari in t
w=ones([1 n]);
w(1)=a(1);
w(2)=b(1);
for i=2:n-1   %calcolo del valore dei pesi associati ai punti
    w(i+1)=b(i)*w(i);
    for j=i:-1:1
        if j==1
            w(j)=a(i)*((i-j+1)/i)*w(j);
        else
            w(j)=a(i)*((i-j+1)/i)*w(j)+b(i)*((j-1)/i)*w(j-1);
        end
    end
end
disp(w)
xy=ginput(n);  %prendo in input i punti del poligono di controllo
%xy=[0,0;-1,1;2,3;1,0]; %punti di test per verificare la correttezza della
%mia implementazione (uguale all'esempio riportato nel paper)
xy1=xy;
hold on; box on; axis equal; axis([0 1 0 1])
plot(xy(:,1),xy(:,2),'x')
plot(xy(:,1),xy(:,2))
a=a./L;
b=b./L;
for j=1:n-1
    pause(2);
    xy(1:n-j,:)=a(j)*tt*xy(1:n-j,:)+t*b(j)*xy((1:n-j)+1,:);
    plot(xy(1:n-j,1),xy(1:n-j,2))
end
plot(xy(1,1),xy(1,2),'x');
disp(xy(1,:))
k=100; %divido l'intervallo parametrico [0,1] in 100 parti dove calcolare
% il valore della curva con lo stesso algoritmo visto in precedenza
t1=t;
B=zeros([k+1,2]);
for t=1:k+1
    B(t,:)=decasteljaurationalcurvecalc((t-1)/k,xy1,a1,b1);
end
plot(B(:,1),B(:,2)) %disegno la curva completa

function b0=decasteljaurationalcurvecalc(t,xy,a,b)
%function che applica l'algortimo visto in precedenza senza graficare i
%risultati parziali delle iterazioni.
n=length(xy(:,1));
tt=1-t;
L=a.*tt+b.*t;
a=a./L;
b=b./L;
for j=1:n-1
    xy(1:n-j,:)=a(j)*tt*xy(1:n-j,:)+t*b(j)*xy((1:n-j)+1,:);
end
b0=xy(1,:);
return;
end