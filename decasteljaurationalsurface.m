clear
%scelgo i punti del net di controllo
xyz=[0 0 1.5; 0 1 1.5; 0 2 1.5; 0 3 1.5;
    1 0 1.5;  1 1 0.5; 1 2 2; 1 3 1.5;
    2 0 1.5; 2 1 1.5; 2 2 1.5; 2 3 1.5];
an=[1,1,1;1/2,10,10;3,2,1/3]; %vettori dei coefficienti lineari per i pesi dei punti delle righe
bn=[1,1,1;10,10,1/2;1/3,2,3];
am=[1/2,2]; %vettori dei coefficienti lineari per i pesi delle righe
bm=[2,1/2];
n=length(am(1,:))+1; %numero di righe
m=length(an(1,:))+1; %numero di colonne
n2=length(an(:,1))+1;
numpunti=m*n; %ricavo il numero di punti necessario
if(n2~=n+1)
    error("La lunghezza dei vettori dei coefficienti lineari non corrisponde!")
end
u=0.5; %parametri in cui valutare la superficie
v=0.5;
%Ricavo i pesi per colonna e i pesi per punto nella riga e li salvo nella
%matrice W
wm=ones([1 n]);
wm(1)=am(1);
wm(2)=bm(1);
disp("Peso righe del net dei punti di controllo: ")
W=zeros(n,m);
for i=2:n-1
    wm(i+1)=bm(i)*wm(i);
    for j=i:-1:1
        if j==1
            wm(j)=am(i)*((i-j+1)/i)*wm(j);
        else
            wm(j)=am(i)*((i-j+1)/i)*wm(j)+bm(i)*((j-1)/i)*wm(j-1);
        end
    end
end
disp(wm');
disp("Peso dei punti del net dei punti di controllo per riga: ")
for c=1:n2-1
    w=ones([1 m]);
    w(1)=an(c,1);
    w(2)=bn(c,1);
    for i=2:m-1
        w(i+1)=bn(c,i)*w(i);
        for j=i:-1:1
            if j==1
                w(j)=an(c,i)*((i-j+1)/i)*w(j);
            else
                w(j)=an(c,i)*((i-j+1)/i)*w(j)+bn(c,i)*((j-1)/i)*w(j-1);
            end
        end
    end
    disp(w/(sum(w)))
    W(c,:)=w.*wm(c)/(sum(w)); %standardizzo i pesi in maniera da poter confrontare anche i pesi
    %che si trovano su righe diverse
end
if(length(xyz)~=numpunti)
    error("Numero di Punti non corretto!");
end
ris=zeros(n,3);
%mostro sul grafico i punti del net di controllo
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'or');
hold on;
xlabel("x");
ylabel("y")
zlabel("z")
pause
%mostro sul grafico le colonne del net di controllo
for i=1:m
    plot3(xyz(i:m:numpunti,1),xyz(i:m:numpunti,2),xyz(i:m:numpunti,3),'y-')
end
for i=1:n
    %mostro le righe del net di controllo
    plot3(xyz(((m)*(i-1))+1:(m)+((m)*(i-1)),1),xyz(((m)*(i-1))+1:(m)+((m)*(i-1)),2),xyz(((m)*(i-1))+1:(m)+((m)*(i-1)),3),'y-')
    %applico l'algoritmo su ogni riga del net
    ris(i,:)=decrazplot(v, xyz(((m)*(i-1))+1:(m)+((m)*(i-1)),:),an(i,:),bn(i,:));
end
plot3(ris(:,1),ris(:,2),ris(:,3),'-.*b'); %mostro sul grafico i punti ottenuti
%applico l'algoritmo una volta sui punti che rappresentano le righe
fin=decrazplot(u,ris,am,bm);
i=0;
j=0;
k=100;
B=zeros([k+1,k+1,3]);
%dividido l'intervallo parametrico [0,1] in 100 parti per i parametri u e v andando ad
%applicare l'algoritmo 10000 volte per ogni coppia di valori u e v
for v=1:k+1
    for u=1:k+1
        B(v,u,:)=rationaldecasteljausurface(xyz,am,bm,an,bn,(u-1)/k,(v-1)/k);
    end
end
%mostro sul grafico il punto sulla superficie ottenuto dall'applicazione dell'algoritmo
plot3(fin(1),fin(2),fin(3),'r*');
pause
%mostro il resto della superficie
surf(B(:,:,1),B(:,:,2),B(:,:,3));
disp("Risultato: ")
disp(fin);


function b0=decrazplot(t,xyz,a,b)
%Function che applica l'algoritmo di tipo De Casteljau per valutare le curve razionali
%di Bezièr in t mostrando nel grafico i risultati parziali delle iterazioni
n=length(xyz(:,1));
tt=1-t;
L=a.*tt+b.*t;
a=a./L;
b=b./L;
for j=1:n-1
    pause
    xyz(1:n-j,:)=a(j)*tt*xyz(1:n-j,:)+t*b(j)*xyz((1:n-j)+1,:);
    plot3(xyz(1:n-j,1),xyz(1:n-j,2),xyz(1:n-j,3)) %stampo i risultati parziali dell'iterazione
end
plot3(xyz(1,1),xyz(1,2),xyz(1,3), 'xr')
b0=xyz(1,:);
return;
end

function fin=rationaldecasteljausurface(xyz,am,bm,an,bn, u, v)
%Function che applica l'algoritmo di tipo De Casteljau per valutare le superfici razionali
%di Bezièr in u e v senza mostrare i risultati parziali.
ml=length(am(1,:))+1;
n1=length(an(1,:))+1;
n2=length(an(:,1))+1;
numpunti=n1*ml;
if(n2~=ml+1)
    error("La lunghezza dei vettori dei coefficienti lineari non corrisponde!")
end
if(length(xyz)~=numpunti)
    error("Numero di Punti non corretto!");
end
ris=zeros(ml,3);
for i=1:ml
    ris(i,:)=decraz(v, xyz(((n1)*(i-1))+1:(n1)+((n1)*(i-1)),:),an(i,:),bn(i,:));
end
fin=decraz(u,ris,am,bm);
end