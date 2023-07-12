% Modelizacion del quench en un iman impregnado
% Caso adiabatico: se desprecia la conveccion
% Se supone el quench en el conductor mas interior
% Se transita dandole una temperatura por encima de la critica a ese elemento
% Imanes en doble pancake con un solo conductor
% Espaciadores
% Modela el conductor con un nucleo de cobre, rodeado por un anillo de superconductor, y un anillo mas exterior de cobre
% Heaters

% clear all

%% 	PARAMETROS DEL MODELO

clear all
% cd('E:\Documentos_E\MCBX\SQUID')

%%%%%Parametros geometricos

%Horizontal
Wc=0.882*0.860*1e-3;       % anchura del conductor desnudo en m %(anchura media 0.86). 0.882= filling factor 
Wah=0.135e-3;	           % espesor del aislante de un conductor en horizontal
Wph=0.002e-3;              % espesor del pegamento de un conductor en horizontal 40
Wnbtih=Wc-0.05e-3;         % anchura de la zona interior del conductor con filamentos de NbTi
Wcuexth=(Wc/2)-(Wnbtih/2); % anchura de la parte exterior de cobre sin NbTi en horizontal
Wcuinh=0.1e-3;             % anchura de la parte interior de cobre sin NbTi en horizontal

%Vertical
Hc=4.37e-3;		           % altura del conductor desnudo en m 
Wav=0.335e-3;	           % espesor del aislante de un conductor en vertical  (real 0.135)0.435e-3
Wpv=0.002e-3;              % espesor del pegamento de un conductor en vertical
Wnbtiv=Hc-0.05e-3;         % altura de la zona interior del conductor con filamentos de NbTi
Wcuextv=(Hc/2)-(Wnbtiv/2); % altura de la parte exterior de cobre sin NbTi en vertical
Wcuinv=0.1e-3;             % altura de la parte interior de cobre sin NbTi en vertical




Nlay=2;		% numero de capas
Nturn=10;    % numero de espiras por capa del dipolo exterior
Lint=4.6;   % longitud espira en capa interior en m 2.6
Lext=4.6;   % longitud espira en capa exterior en m 2.6



%%%%%Parametros del circuito electrico

I(1)=1457;                  % corriente inicial en A del dipolo interior. 
                            %nom 1625. 
                            %ultimate 1742 A.
                            
B0=1.4142*2.23*I(1)/1625;	% 2.23 T. 1.4142* para adicion de campo de OD.  
                            %campo magnetico inicial en conductores en T (luego matriz?)
                            %B0=2.24*I(1)/1630; 
                            %B0=2.46*I(1)/1790; 
                            %B0=2.73*I(1)/1630; 
                            %*1.4142= raiz de 2, ambos dipolos en corriente
                            
B=B0;                       % campo magnetico segï¿½n se produce el quench

% Inductancia variable segun baja la corriente debido a saturaciï¿½n en hierro
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Inductancia Variable\ILV1.csv";
data=readmatrix(nombreArchivo);
IL=data(:, 1);
Lv=data(:, 2);

Rdamp=0.15;     % resistencia exterior de proteccion. Resistencia del circuito Rdamp=0.007 at 30 ms;	
tdamp=0.027;	% tiempo de retardo de conexion de resistencia exterior de proteccion;
HeaterStatus=1; % variable de estado de los heaters (0 están listos para disparar, 1 están inhibidos)
Vdet=0.1;       % voltaje de deteccion de quench
tdet=0.013;     % inicializacion de variable para tiempo de deteccion de quench
tval=0.01;      % tiempo para validacion de quench
tdel=0.014;     % tiempo de calentamiento de los heaters (desde validation hasta quench en el cable)

%%%%%Parametros fisicos de los materiales

Rsc=1-(1.75/(1.75+1));                   % relacion volumetrica de superconductor
Rcuinttot=Wcuinh*Wcuinv/(Wc*Hc*(1-Rsc)); % relacion entre el area interior de cobre y todo el cobre del cable
Rins=0.4;                                %?% porcentaje de aislamiento en la seccion recta de la bobina
Rfillet=0.01;                            % porcentaje de redondeo que tiene en cuenta que la seccion metalica no es rectangular Hc*Wc
To=1.9;                                  % temperatura del baño en K
% Jco=2584e6;                            % J critica en NbTi 5T y 4.2K
Jco=3000e6;                              %ok
Lo=2.45e-8;                              % constante de Lorentz
%tempLo=[0,5,8,10,20,30,40,50,60,70,80,90,100,200,300];
%Lov=[2.45,2.45,2.2,2.4,2,1.57,1.38,1.35,1.35,1.4,1.5,1.6,1.69,2.2,2.32].*1e-8;
Jc=Jco*(3.57-0.377*To+(0.012*To-0.25)*B);% Jc NbTi segun Spigo
% Jc=(-900*(B-4) +5600)*1E6;             % cable supercon
Ic=Jc*Wc*Hc*Rsc;                         % intensidad critica en parte metalica
Tc=9.4*(1-B/14.5)^0.59;                  % temperatura critica
Tcs=To+(Tc-To)*(1-I(1)/Ic);              % temperatura current-sharing
RRR=200;                                 % triple R. Measured 268. Initial considered 120.		                            % triple R. Measured 268. Initial considered 120.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% capacidad calorifica del cobre
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Capacidad Calorifica Cu\CC_Cu.csv";
data=readmatrix(nombreArchivo);
tempcuv=data(:, 1);
cpcuv=data(:, 2);

% capacidad calorifica del NbTi (CERN)
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Capacidad Calorifica NbTi\CC_NbTi_CERN.csv";
data=readmatrix(nombreArchivo);
tempnbtiv=data(:, 1);
cpnbtiv=data(:, 2);

% capacidad calorifica del G10
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Capacidad Calorifica G10\CC_G10.csv";
data=readmatrix(nombreArchivo);
tempg10=data(:, 1);
cpg10=data(:, 2);

%capacidad calorï¿½fica del adhesivo epoxy ISR-D1
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Capacidad Calorifica EPOXY\CC_EPOXY_Cryocomp.csv";
data=readmatrix(nombreArchivo);
tempcpEpo=data(:, 1);
cpEpo=data(:, 2);

%capacidad calorï¿½fica del barniz
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Capacidad Calorifica Barniz\CC_Barniz_GE7031_iwasa.csv";
data=readmatrix(nombreArchivo);
tempvar=data(:, 1);
cpvar=data(:, 2);



% NbTi thermal conductivity (GSI)
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Conductividad Termica NbTi\CT_NbTi_GSI.csv";
data=readmatrix(nombreArchivo);
tempnbtik=data(:, 1);
Knbtik=data(:, 2);

% conductividad termica del barniz GE
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Conductividad Termica Barniz\CT_Barniz_GE.csv";
data=readmatrix(nombreArchivo);
tempGE=data(:, 1);
KaGE=data(:, 2);

% conductividad termica de la fibra de vidrio
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Conductividad Termica Vidrio\CT_Glass.csv";
data=readmatrix(nombreArchivo);
tempGlass=data(:, 1);
KaGlass=data(:, 2);

% conductividad termica del EPOXY
nombreArchivo="X:\Carpetas personales\Luis_Gonzalez\Quench\Conductividad Termica EPOXY\CT_EPOXY.csv";
data=readmatrix(nombreArchivo);
tempKaepo=data(:, 1);
Kaepo=data(:, 2);

% 	Parï¿½metros de la malla de diferencias finitas
h=1e-3;             % paso de integracion (s)
m=10;               % 60 numero de cortes transversales. m=80, Tmax=180 K vs 182 K before, negligible. m=20, Tmax=200 K, qualitatively like Roxie
toler=.01;          % tolerancia para convergencia de temperaturas en K
VecesEstancado=0;   % Veces que la soluciï¿½n del sistema no ha convergido
VecesNoConv=0;

%% VARIABLES

Ne=Nturn*Nlay;				% numero de nodos en cada corte ( o de espiras)
T=zeros(Ne*m,1);			% tension inicial
Tnew=To*ones(Ne*m,1);       % tension en el paso siguiente

NumIt=2e4; %Numero mï¿½x. de iteraciones esperado para evitar crecimiento de variables en el loop
t=zeros(1,NumIt);       %Tiempo
intII=zeros(1,NumIt);
intE=zeros(1,NumIt);
I=I+zeros(1,NumIt);     %Corriente
pI=zeros(1,NumIt);
V=zeros(1,NumIt);       % Voltage
V1=zeros(1,NumIt);      % Voltage
V2=zeros(1,NumIt);      % Voltage
V3=zeros(1,NumIt);      % Voltage
R=zeros(1,NumIt);       % Resistencia
VlongCalculado=0;       % Velocidad longitudinal todavï¿½a no calculada
VtransHCalculado=0;     % Velocidad transversal horizontal todavï¿½a no calculada
VtransVCalculado=0;     % Velocidad transversal vertical todavï¿½a no calculada

% Elecciï¿½n del nodo inicio del quench

quenchPoint="Pole"

if quenchPoint=="Pole"
    ElemQuench=1;
elseif quenchPoint=="Midplane"
    ElemQuench=Nlay;
else quenchPoint=="Middle"
    ElemQuench=Ne*idivide(int32(m),int32(2))+idivide(int32(Nturn),int32(2))+Nturn*idivide(int32(Nlay),int32(2)); % elegimos un elemento del centro de la bobina;
end

Tnew(ElemQuench)=Tcs;       % ponemos su temperatura a la temperatura de current sharing
%Tnew(ElemQuench+Ne)=Tcs;   % elemento siguiente por si  todavia no se llega al minimo tamaño de quench
%Tnew(ElemQuench+2*Ne)=Tcs; % elemento siguiente por si  todavia no se llega al minimo tamaño de quench
Etran=Tnew >= Tcs;          % condicion de elementos transitados

% Proceso iterativo de resoluciï¿½n del problema
%Tmax=zeros(m,1);	% matriz de temperaturas maximas
Tmax=zeros(NumIt,1);
ii=0;               % contador de iteraciones

%Matriz campo magnético en cada punto
Bmat=zeros(Ne*m,1);
Baux=ones(Nturn,1)*ones(1,Nlay);
%Baux=[linspace(0.9,1,29),linspace(1,0.9,29)]'*[linspace(1,0.01,3*idivide(int32(Nlay),int32(4))),linspace(0.01,0.4,Nlay-3*idivide(int32(Nlay),int32(4)))];
Baux=reshape(Baux,Nturn*Nlay,1);
Baux=Baux*ones(1,m);
Baux=reshape(Baux,Nturn*Nlay*m,1);
Bmat=B*Baux;



Lturn=zeros(Ne*m,1);
Lturn=ones(Nturn,1)*linspace(Lint,Lext,Nlay);
Lturn=reshape(Lturn,Nturn*Nlay,1);
Lturn=Lturn*ones(1,m);
Lturn=reshape(Lturn,Nturn*Nlay*m,1);

%% CONSTRUCCION DE MATRICES

% Definicion de la caja de la matriz de coeficientes (corte transversal).
% Se generan dos matrices A_v, A_h y Aux con los coeficientes de la conducciï¿½n

% DIAGONAL HORIZONTAL -------------------
vaux=2*ones(Ne,1);			% vector auxiliar para diagonal
vaux(1:Nturn)=vaux(1:Nturn)-ones(Nturn,1);
vaux(Ne-Nturn+1:Ne)=vaux(Ne-Nturn+1:Ne)-ones(Nturn,1);
Amh=sparse(1:Ne,1:Ne,vaux,Ne,Ne,Ne*13); %Pone vaux en los elementos de la diagonal de una matriz dispersa (todo lo demï¿½s=0)

figure;
imagesc(Amh);
title('Diagonal Horizontal');

% elementos de conduccion horizontal entre conductores,
% y se añaden los elementos simetricos en la matriz
Amh=Amh+sparse(1:Ne-Nturn,Nturn+1:Ne,-ones(Ne-Nturn,1),size(Amh,1),size(Amh,2));
Amh=Amh+sparse(Nturn+1:Ne,1:Ne-Nturn,-ones(Ne-Nturn,1),size(Amh,1),size(Amh,2));

figure;
imagesc(Amh);
title('Diagonal Horizontal con elementos de conducción horizontal');


% ensamblamos las cajas para tener la matriz completa horizontal A_h
[Nfil,Ncol,valor]=find(Amh);
A_h=sparse(Nfil,Ncol,valor,Ne*m,Ne*m,Ne*m*13);
% repite la matriz del corte transversal horizontal m veces
for jj=1:m-1
    A_h=A_h+sparse(Nfil+jj*Ne,Ncol+jj*Ne,valor,size(A_h,1),size(A_h,2));
end

figure;
imagesc(A_h);
title('Cajas ensambladas horizontal');



% DIAGONAL VERTICAL -------------------

vaux=2*ones(Ne,1);			% vector auxiliar para diagonal
vaux(1:Nturn:Ne)=vaux(1:Nturn:Ne)-ones(Nlay,1);
vaux(Nturn:Nturn:Ne)=vaux(Nturn:Nturn:Ne)-ones(Nlay,1);
Amv=sparse(1:Ne,1:Ne,vaux,Ne,Ne,Ne*13); %Pone vaux en los elementos de la diagonal de una matriz dispersa (todo lo demï¿½s=0)

figure;
imagesc(Amv);
title('Diagonal Vertical');

% elementos de conduccion vertical entre conductores y
% tambien los simetricos en la matriz
Maux=sparse(1:Ne-1,2:Ne,-ones(Ne-1,1),Ne,Ne,Ne*13);
Maux=Maux+sparse(Nturn:Nturn:Ne-Nturn,Nturn+1:Nturn:Ne-Nturn+1,ones(Nlay-1,1),size(Maux,1),size(Maux,2));
Maux=Maux+sparse(2:Ne,1:Ne-1,-ones(Ne-1,1),size(Maux,1),size(Maux,2));
Maux=Maux+sparse(Nturn+1:Nturn:Ne-Nturn+1,Nturn:Nturn:Ne-Nturn,ones(Nlay-1,1),size(Maux,1),size(Maux,2));
Amv=Amv+Maux;

figure;
imagesc(Amv);
title('Diagonal Verical con elementos de conducción vertical');



% ensamblamos las cajas para tener la matriz completa vertical A_v
[Nfil,Ncol,valor]=find(Amv);
A_v=sparse(Nfil,Ncol,valor,Ne*m,Ne*m,Ne*m*13);
% repite la matriz del corte transversal vertical m veces
for jj=1:m-1
    A_v=A_v+sparse(Nfil+jj*Ne,Ncol+jj*Ne,valor,size(A_v,1),size(A_v,2));
end

figure;
imagesc(A_v);
title('Cajas ensambladas vertical');


%A_aux=A_v;
%A_v=A_h;
%A_h=A_aux;


% hacemos ahora las conexiones en serie (conductividad longitudinal)
Aux=sparse(Nfil,Ncol,0,size(A_h,1),size(A_h,2),Ne*m*5);
% entre cortes transversales
for jj=1:m-1
    Aux=Aux+sparse(1+(jj-1)*Ne:jj*Ne,1+jj*Ne:(jj+1)*Ne,-ones(Ne,1),size(A_h,1),size(A_h,2));
    Aux=Aux+sparse(1+jj*Ne:(jj+1)*Ne,1+(jj-1)*Ne:jj*Ne,-ones(Ne,1),size(A_h,1),size(A_h,2));
end

figure;
imagesc(Aux);
title('Conexiones en serie Conductividad Longitudinal');

% bobina windowframe
%Maux=sparse(2:Nturn,1+(m-1)*Ne:(m-1)*Ne+Nturn-1,-ones(Nturn-1,1),size(A_h,1),size(A_h,2)); % conectamos la primera capa con el ultimo corte
%Maux=Maux+sparse(1+(m-1)*Ne:(m-1)*Ne+Nturn-1,2:Nturn,-ones(Nturn-1,1),size(A_h,1),size(A_h,2)); % simetricos
%for jj=2:Nlay
%    if rem(jj,2)==0 %si el numero de capa es par
%        Maux=Maux+sparse(jj*Nturn,(m-1)*Ne+(jj-1)*Nturn,-ones(1,1),size(A_h,1),size(A_h,2)); % salto de capa
%        Maux=Maux+sparse((m-1)*Ne+(jj-1)*Nturn,jj*Nturn,-ones(1,1),size(A_h,1),size(A_h,2)); % salto de capa simetrico
%        Maux=Maux+sparse((jj-1)*Nturn+1:(jj)*Nturn-1,2+(m-1)*Ne+(jj-1)*Nturn:(m-1)*Ne+(jj)*Nturn,-ones(Nturn-1,1),size(A_h,1),size(A_h,2));
%        Maux=Maux+sparse(2+(m-1)*Ne+(jj-1)*Nturn:(m-1)*Ne+(jj)*Nturn,(jj-1)*Nturn+1:(jj)*Nturn-1,-ones(Nturn-1,1),size(A_h,1),size(A_h,2)); %simetrico
%    else % si es impar
%        Maux=Maux+sparse((jj-1)*Nturn+2:(jj)*Nturn,1+(m-1)*Ne+(jj-1)*Nturn:(m-1)*Ne+(jj)*Nturn-1,-ones(Nturn-1,1),size(A_h,1),size(A_h,2));
%        Maux=Maux+sparse(1+(m-1)*Ne+(jj-1)*Nturn:(m-1)*Ne+(jj)*Nturn-1,(jj-1)*Nturn+2:(jj)*Nturn,-ones(Nturn-1,1),size(A_h,1),size(A_h,2)); %simetrico
%        Maux=Maux+sparse((jj-1)*Nturn+1,(m-1)*Ne+(jj-2)*Nturn+1,-ones(1,1),size(A_h,1),size(A_h,2)); %salto de capa
%        Maux=Maux+sparse((m-1)*Ne+(jj-2)*Nturn+1,(jj-1)*Nturn+1,-ones(1,1),size(A_h,1),size(A_h,2)); %salto de capa simetrico
%    end
%end

% hacemos ahora las conexiones electricas para doble galleta
Nrib=1;
	% entre conductores de una misma capa de abajo
	Maux=     sparse((m-1)*Ne+1:m*Ne-Nrib*Nturn-1,2:Nrib*Nturn,-ones(Ne/2-1,1),size(Aux,1),size(Aux,2));
	Maux=Maux+sparse(2:Nrib*Nturn,(m-1)*Ne+1:m*Ne-Nrib*Nturn-1,-ones(Ne/2-1,1),size(Aux,1),size(Aux,2));

	% le quitamos los finales de capa de abajo
	Maux=Maux+sparse((m-1)*Ne+Nturn:Nturn:(m-0.5)*Ne-1,Nturn+1:Nturn:Nrib*Nturn,ones(Nrib-1,1),size(Aux,1),size(Aux,2));
	Maux=Maux+sparse(Nturn+1:Nturn:Nrib*Nturn,(m-1)*Ne+Nturn:Nturn:(m-0.5)*Ne-1,ones(Nrib-1,1),size(Aux,1),size(Aux,2));

	% entre conductores de una misma capa de arriba
	Maux=Maux+sparse(Nrib*Nturn+1:Ne-1,(m-1)*Ne+Nrib*Nturn+2:m*Ne,-ones(Ne/2-1,1),size(Aux,1),size(Aux,2));
	Maux=Maux+sparse((m-1)*Ne+Nrib*Nturn+2:m*Ne,Nrib*Nturn+1:Ne-1,-ones(Ne/2-1,1),size(Aux,1),size(Aux,2));

	% le quitamos los finales de capa de arriba
	Maux=Maux+sparse((Nrib+1)*Nturn:Nturn:Ne-1,(m-1)*Ne+(Nrib+1)*Nturn+1:Nturn:m*Ne-1,ones(Nrib-1,1),size(Aux,1),size(Aux,2));
	Maux=Maux+sparse((m-1)*Ne+(Nrib+1)*Nturn+1:Nturn:m*Ne-1,(Nrib+1)*Nturn:Nturn:Ne-1,ones(Nrib-1,1),size(Aux,1),size(Aux,2));

	% cambio de capa
	Maux=Maux+sparse((m-1)*Ne+Nrib*Nturn+1:Nturn:m*Ne,1:Nturn:Nrib*Nturn,-ones(Nrib,1),size(Aux,1),size(Aux,2));
	Maux=Maux+sparse(1:Nturn:Nrib*Nturn,(m-1)*Ne+Nrib*Nturn+1:Nturn:m*Ne,-ones(Nrib,1),size(Aux,1),size(Aux,2));

	% soldaduras interiores
	%Maux=Maux+sparse((Nrib+1)*Nturn:Nturn:Ne-1,(m-1)*Ne+Nturn*2:Nturn:(m-1)*Ne+Nrib*Nturn,-ones(Nrib-1,1),size(A,1),size(A,2));
	%Maux=Maux+sparse((m-1)*Ne+Nturn*2:Nturn:(m-1)*Ne+Nrib*Nturn,(Nrib+1)*Nturn:Nturn:Ne-1,-ones(Nrib-1,1),size(A,1),size(A,2));

%

Aux=Aux+Maux; % formamos la matriz de coeficientes de conducciï¿½n longitudinal

figure;
imagesc(Aux);
title('Conexiones en serie Coef Conductividad Longitudinal');

% Waitbar initialization
waitbar_handle = waitbar(0,'Computing...');
set(waitbar_handle,'Units','Normalized',...
    'Position', [0.3869    0.4692    0.2250    0.1] )
get(waitbar_handle,'Position')


%% BUCLE DE SOLUCION
tic
while((I(ii+1)/I(1))>1e-2)
    T=Tnew;
    ii=ii+1;
    t(ii+1)=t(ii)+h;
    
    B=I(ii)*B0/I(1); % variaciï¿½n aproximada del campo magnetico al caer la corriente
    Bmat=B*Baux;    % matriz con campo magnetico para solenoide largo
    L=interp1(IL,Lv,I(ii)); % variaciï¿½n de la autoinducciï¿½n al caer la corriente
    Jc=Jco*(3.57-0.377*To+(0.012*To-0.25)*Bmat); 	% recalcular Jc
    % Jc=(-526.85*Bmat+5314.4)*1E6;
    Ic=Jc*Wc*Hc*Rsc;            % recalcular Ic
    Tc=9.4*(1-Bmat./14.5).^0.59;		% recalcular la temperatura critica
    Tcs=To+(Tc-To).*(1-I(ii)./Ic);	% recalcular la temperatura current-sharing
    
    % Propiedades termicas de los materiales
    % resistividad del cobre en funcion de RRR,T,B
    roo=15.53e-9/RRR;
    roi=1.171e-17*T.^4.49./(1+4.498e-7*T.^3.35.*exp(-(50./T).^6.428));
    roio=0.4531*roo.*roi./(roo+roi);
    rocuo=roo+roi+roio;
    
    s=15.53e-9*Bmat./rocuo;
    r=10.^(-2.662+0.3168.*log10(s)+0.6229.*log10(s).^2- ...
        0.1839*log10(s).^3+0.01827*log10(s).^4);
    rocu=(1+r).*rocuo.*(Bmat > 0.1)+rocuo.*(Bmat <= 0.1);
    
    
    % CONDUCTIVIDADES TERMICAS para cada elemento a la temperatura de esta iteraciï¿½n
    
    % conductividad del cobre en funcion de RRR,T,B
    %Lo=interp1(tempLo,Lov,T);
    Kcu=Lo.*T./rocu;
    %Kcu=10.^((2.2154+(-0.88068).*T.^0.5+0.29505.*T+(-0.04831).*T.^1.5+0.003207.*T.^2)./(1+(-0.47461).*T.^0.5+0.13871.*T+(-0.02043).*T.^1.5+0.001281.*T.^2));
    
    % conductividad del NbTi en funcion de T
    Knbti=interp1(tempnbtik,Knbtik,T);
    
    % conductividad longitudinal del cable en funcion de RRR,T,B (sï¿½lo parte metï¿½lica)
    Kc=(1-Rsc)*Kcu+Rsc*Knbti;
    
    % conductividad barniz GE
    Kpva=interp1(tempGE,KaGE,T);
    
    % conductividad de la fibra de vidrio
    Kglass=interp1(tempGlass,KaGlass,T);
    
    % conductividad del epoxy
    Kepo=interp1(tempKaepo,Kaepo,T);
    
    % CAPACIDAD CALORIFICA para cada elemento a la temperatura de esta iteracion
    
    % capacidad calorifica del cobre a la temperatura de esta iteracion
    Ccu=interp1(tempcuv,cpcuv,T);
    
    % capacidad calorï¿½fica del NbTi
    Cnbti=interp1(tempnbtiv,cpnbtiv,T);
    
    % capacidad calorï¿½fica del epoxy
    Cepo=interp1(tempcpEpo,cpEpo,T);
   
    % capacidad calorï¿½fica del barniz GE 7031 (Iwasa + extrapolacion)
    Cpva=interp1(tempvar,cpvar,T);
    
    % capacidad calorï¿½fica del G10
    Cg10=interp1(tempg10,cpg10,T);
    
%%%%% capacidad calorï¿½fica del cable
%     Cc=((1-Rsc)*Ccu+Rsc*Cnbti)*(1-Rins)+Cepo*(Rins);%Wrong aisl. de epoxy
%     Cc=((1-Rsc)*Ccu+Rsc*Cnbti)*(1-Rins)+(0.5*Cg10+0.5*Cepo)*(Rins); %Aislamiento de G10 con epoxy aumentado (0.5*Cg10+0.5*Cepo)
    Cc=((1-Rsc)*Ccu+Rsc*Cnbti)*(1-Rins)+Cg10*(Rins); %Aislamiento de G10 
    %Cc=((1-Rsc)*Ccu+Rsc*Cnbti);
    
    
    
    
    % PARAMETROS DEL CIRCUITO ELECTRICO EQUIVALENTE
    % Resistencia y admitancia de conduccion longitudinal (sï¿½lo parte metï¿½lica)
    Rcl=Lturn./(m*Kc*Wc*Hc*(1-Rfillet));
    Ycl=1./Rcl;
    
    % Resistencia y admitancia de conduccion transversal en vertical
    Repoxy_v=Wpv./(Kepo*Wc*0.8.*Lturn/m);
%     Rpva_v=Wav./(Kpva*Wc*0.8.*Lturn/m); % Wrong: This cable has not pva glue
%     Rpva_v=Wav./((0.5*Kglass+0.5*Kepo)*Wc*0.8.*Lturn/m); % ~KG10 con epoxy aumentado:(0.5*Kglass+0.5*Kepo)
    Rpva_v=Wav./(Kglass*Wc*0.8.*Lturn/m); % Kglass es KG10 con fit de NIST
    Rcuext_v=Wcuextv./(Kcu*Wc*0.8.*Lturn/m);
    Rnbti_v=(Wnbtiv/2-Wcuinv/2)./(Knbti*Wnbtih.*Lturn/m);
    Rcuin_v=(Wcuinv/2)./(Kcu*Wcuinh.*Lturn/m);
    Rctv=2*(Repoxy_v+Rpva_v+Rcuext_v+Rnbti_v*Rcuinttot+Rcuin_v*Rcuinttot);
    %Rctv=(2*Wav./Kglass+Hc./Kc)./(0.25*Wc*Lturn/m); % Calculo antiguo de resistencia transversal vertical
    Yctv=1./Rctv;
    
    % Resistencia y admitancia de conduccion transversal en horizontal
    Repoxy_h=Wph./(Kepo*Hc*0.8.*Lturn/m);
%     Rpva_h=Wah./(Kpva*Hc*0.8.*Lturn/m); % Wrong: This cable has not pva glue
%     Rpva_h=Wah./((0.5*Kglass+0.5*Kepo)*Hc*0.8.*Lturn/m); % ~KG10 con epoxy aumentado:(0.5*Kglass+0.5*Kepo)
    Rpva_h=Wah./(Kglass*Hc*0.8.*Lturn/m); % Kglass es KG10 con fit de NIST
    Rcuext_h=Wcuexth./(Kcu*Hc*0.8.*Lturn/m);
    Rnbti_h=(Wnbtih/2-Wcuinh/2)./(Knbti*Wnbtiv.*Lturn/m);
    Rcuin_h=(Wcuinh/2)./(Kcu*Wcuinv.*Lturn/m);
    Rcth=2*(Repoxy_h+Rpva_h+Rcuext_h+Rnbti_h*Rcuinttot+Rcuin_h*Rcuinttot);
    %Rcth=(2*Wah./Kglass+Wc./Kc)./(0.25*Hc*Lturn/m); % Calculo antiguo de resistencia transversal horizontal
    Ycth=1./Rcth;
    
    Re=rocu.*(Lturn/m)/(Wc*Hc*(1-Rsc)*(1-Rfillet));		% resistencia de los elementos transitados en ohm
    C=Cc*(Wc+2*Wah)*(Hc+2*Wav)*cos(30*pi/180).*Lturn/m; 	% capacidad equivalente en F
    % C=Cc*Wc*Hc*cos(30*pi/180).*Lturn/m; 	% capacidad equivalente en F
    
    % Se meten los datos de las admitancias en la matriz AA segï¿½n los elementos de A_h A_v y Aux.
    % se pone a 0 la matriz AA en cada iteraciï¿½n
    AA=sparse(1,1,0,Ne*m,Ne*m,Ne*m*13);
    
    % Se construyen y se suman
    % multiplicamos por el valor de la admitancia transversal horizontal
    [Nfil,Ncol,valor]=find(A_h);
    valor=valor*0.5.*(Ycth(Nfil)+Ycth(Ncol));
    AA=AA+sparse(Nfil,Ncol,valor,size(AA,1),size(AA,2));
    
    % multiplicamos por el valor de la admitancia transversal vertical
    [Nfil,Ncol,valor]=find(A_v);
    valor=valor*0.5.*(Yctv(Nfil)+Yctv(Ncol));
    AA=AA+sparse(Nfil,Ncol,valor,size(AA,1),size(AA,2));
    
    % multiplicamos por el valor de la admitancia longitudinal
    [Nfil,Ncol,valor]=find(Aux);
    valor=valor*0.5.*(Ycl(Nfil)+Ycl(Ncol));
    AA=AA+sparse(Nfil,Ncol,valor,size(AA,1),size(AA,2));
    
    % sumamos el termino de admitancias de conduccion longitudinales en la diagonal
    AA=AA+sparse(1:Ne*m,1:Ne*m,2*Ycl,size(AA,1),size(AA,2));
    % Quitamos una admitancia en la diagonal del primer elemento (no conectado a nada)
    AA(1,1)=AA(1,1)-Ycl(1);
    % Quitamos una admitancia en la diagonal del ï¿½ltimo elemento (solo para el caso de bobinado windowframe)
    AA(Ne*m-Nturn+1,Ne*m-Nturn+1)=AA(Ne*m-Nturn+1,Ne*m-Nturn+1)-Ycl(Ne*m-Nturn+1);
    
    % sumamos el termino de capacidad calorifica
    AA=AA+sparse(1:Ne*m,1:Ne*m,C./h,size(AA,1),size(AA,2));
    
    % Calor generado en la zona transitada
    q=I(ii)^2*Re;
    % solo los elementos que no son superconductores
    var=sparse(C./h.*T+q.*Etran);
    [Tnew,bandera,NULL,iter] = bicgstab(AA,var,1e-8,200,[],[],T);
    h=1/(40*iter); % modificamos el tamaï¿½o de paso segï¿½n la convergencia anterior
    switch bandera
        case 1
            sprintf('Metodo no converge en el nï¿½mero mï¿½ximo de iteraciones %d',iter)
            VecesNoConv=VecesNoConv+1;
        case 3
            sprintf('Metodo estancado en iteracion %d',iter)
            VecesEstancado=VecesEstancado+1;
        case 4
            sprintf('Matriz mal condicionada en iteraciï¿½n %d',iter)
        otherwise
    end
    
    % Quench heaters en t=0.1
    %     if (t(ii+1)>(0.1-(h/2)))&&(t(ii+1)<(0.1+(h/2)))
    %         Tnew=Tnew+Tcs;
    %     end
    
        % Quench heaters in t=tdet+tval+tdel
    if (tdet~=0 && HeaterStatus==0 && t(ii)>tdet+tval+tdel)
    %  vauxh=zeros(Ne,1);                                       % sin qh
%        vauxh=[zeros(5,1);ones(Nturn-5,1);zeros(Ne-Nturn,1)];  % qh en bloques 4, 5
%      vauxh=[zeros(42,1);ones(Nturn-42,1);zeros(Ne-Nturn,1)];  % qh en bloque 4
%      vauxh=[zeros(14,1);ones(Nturn-42,1);zeros(28,1);zeros(Ne-Nturn,1)];  % qh en bloque 5
%        vauxh=vauxh*ones(1,m);                                   % qh con bloques simetricos
%      vauxh=vauxh*[ones(1,floor(m/2)),zeros(1,ceil(m/2))];  % qh sin bloques simetricos
%        vauxh=reshape(vauxh,Ne*m,1);
%        Tnew=Tnew+Tcs.*~Etran.*vauxh;
       Tnew=Tnew+Tcs;
       Re=Re.*2;                     % transitan las dos bobinas del dipolo por los heaters
       HeaterStatus=1;
    end
    
    % Elementos transitados y calculo de resistencia total
    Etran=Tnew > Tcs;
    R(ii)=sum(Re.*Etran);
    
%     Calculo velocidades. Vï¿½lido sï¿½lo para quenches iniciados en el primer corte transversal (m=1)
%     La veloc. transversal vertical solo es correcta si es mï¿½s rï¿½pida que la longitudinal hasta la siguiente vuelta
%     Mtemp=reshape(1:Nturn*Nlay,Nturn,Nlay);
%     [fila,colum]=find(Mtemp==ElemQuench);
%     Calculo de la velocidad longitudinal media de quench a un tercio de las divisiones mï¿½ximas
%     if ((Etran(ElemQuench+(idivide(int32(m),int32(2)))*Ne)==1)&&(VlongCalculado==0))
%         Vlong=(double(idivide(int32(m),int32(2))).*Lturn/m)/t(ii+1);
%         VlongCalculado=1;
%     end
%     Calculo de la velocidad horizontal media de quench
%     if ((Etran(Mtemp(fila,Nlay))==1)&&(VtransHCalculado==0))
%         VtransH=(Nlay-colum)*(Wc+2*Wah+2*Wph)/t(ii+1);
%         VtransHCalculado=1;
%     end
%     Calculo de la velocidad vertical media de quench
%     if ((Etran(Mtemp(Nturn,colum))==1)&&(VtransVCalculado==0))
%         VtransV=(Nturn-fila)*(Hc+2*Wav+2*Wpv)/t(ii+1);
%         VtransVCalculado=1;
%     end
%     
%     Calculo del tiempo en que toda la bobina ha quencheado
%     if min(Etran)==0
%         TiempoQuenchCompleto=t(ii+1);
%     end
    
    % ecuacion del circuito electrico L+R. Variaciï¿½n de la intensidad respecto del tiempo
    % if t(ii+1)>tdamp
    if V3(ii)>=Vdet
        if tdet==0
            tdet=t(ii);
        end
         if t(ii)>tdet+tval
             pI(ii)=-(R(ii)+0.)*I(ii)/L;
         end
          if t(ii)>tdamp+tdet+tval
            pI(ii)=-(R(ii)+Rdamp)*I(ii)/L;
        end
    else
        pI(ii)=-(R(ii)+0.)*I(ii)/L;
    end
    
    % actualizacion
    I(ii+1)=I(ii)+pI(ii)*h;
    V(ii+1)=L*pI(ii); % Tension de todas las bobinas juntas debido a su inductancia
    %V1(ii+1)=I(ii+1)*sum(Re(1:Ne*m).*Etran(1:Ne*m))+L/4*pI; % Tension de la bobina que quenchea (inductiva + resistiva). Caso de quadrupolo.
    V3(ii+1)=I(ii+1)*sum(Re(1:Ne*m).*Etran(1:Ne*m));
    
    % integral de MIITS
    intII(ii+1)=intII(ii)+I(ii)*I(ii)*h;
    
    % energia disipada
    intE(ii+1)=intE(ii)+I(ii)*I(ii)*R(ii)*h;
    
    Tmax(ii)=max(Tnew);
    
    % Updating waitbar
    waitbar((I(1)-I(ii))/I(1),waitbar_handle,['Iter: ' num2str(ii) char(10)...
        'Total number of transited elem: ' num2str(nnz(Etran)) char(10)...
        'IterConv: ' num2str(iter) char(10)...
        't = ' num2str(t(ii+1)) 's    I = ' num2str(I(ii)) 'A' char(10)...
        'Temp Max = ' num2str(Tmax(ii)) '    Energy dissipated = ' num2str(intE(ii+1))]);
    
    %     sprintf('Iteraciï¿½n nï¿½mero: %i, Nï¿½mero total de elementos transitados: %i, IterConv: %i',ii,nnz(Etran),iter)
    %     sprintf('Tiempo: %.4f, Intensidad: %.3f, Max. Temp.: %.2f, Energï¿½a disipada: %.2f',t(ii+1),I(ii),Tmax(ii),intE(ii+1))
    %     find(Etran)
    %     pause
end
toc



%% RESULTADOS
MIITS=intII(ii+1);
sprintf('MIITs: %.2f',MIITS)
sprintf('Dissipated total energy: %.2f',intE(ii+1))
%sprintf('Velocidad longitudinal: %.2f; Velocidad transversal horizontal: %.3f; Velocidad transversal vertical: %.3f',Vlong,VtransH,VtransV)

% Plot settings
%**************************************************************************
norm_size=0.04;
plot_title='Quench Simulation Results';

% Figure Generation
%**************************************************************************

% Fer
% figure
% plot(t(2:ii+1),I(1:ii),'+y');
% xlabel('Time (s)')
% ylabel('Current (A)')
% title('Current decay in quench')
% figure
% plot(t(2:ii+1),Tmax(1:ii),'r');
% xlabel('Time (s)')
% ylabel('Temperature (K)')
% title('Maximum temperature in quench')
% figure
% plot(t(2:ii+1),R(1:ii),'+');
% xlabel('Time (s)')
% ylabel('Resistance (Ohm)')
% title('Coil resistance in quench')
% figure
% plot(t(2:ii+1),V3(1:ii),'+');
% xlabel('Time (s)')
% ylabel('Voltage (V)')
% title('Quenching coil voltage')
% figure
% plot(t(2:ii+1),intII(1:ii),'+');
% figure
% plot(t(2:ii+1),intII(1:ii)/(Wc*Hc)^2,'+');


%Chus
h_fig=figure('Name',plot_title,...
    'Units','Normalized',...
    'OuterPosition',[0.025,0.04,0.95,0.95],...
    'NumberTitle','off');

% Axis and labels
h_ax(1)=axes('Parent',h_fig,...
    'Position',[0.05 0.55 0.4 0.4]);
plot(h_ax(1),t(2:ii+1),I(1:ii),'+y');


title('Parent',h_ax(1),'Current decay in quench',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);
ylabel(h_ax(1),'Current (A)',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);


h_ax(2)=axes('Parent',h_fig,...
    'Position',[0.55 0.55 0.4 0.4]);
plot(h_ax(2),t(2:ii+1),Tmax(1:ii),'r');

title('Parent',h_ax(2),'Maximum temperature in quench',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);
ylabel(h_ax(2),'Temperature (K)',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);


h_ax(3)=axes('Parent',h_fig,...
    'Position',[0.05 0.05 0.4 0.4]);
plot(h_ax(3),t(2:ii+1),R(1:ii),'+');

title('Parent',h_ax(3),'Coil resistance in quench',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);
ylabel(h_ax(3),'Resistance (Ohm)',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);


h_ax(4)=axes('Parent',h_fig,...
    'Position',[0.55 0.05 0.4 0.4]);
plot(h_ax(4),t(2:ii+1),V3(1:ii),'+');

title('Parent',h_ax(4),'Quenching coil voltage',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);
ylabel(h_ax(4),'Voltage (V)',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);


xlabel(h_ax(1),'Time (s)',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);
xlabel(h_ax(2),'Time (s)',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);
xlabel(h_ax(3),'Time (s)',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);
xlabel(h_ax(4),'Time (s)',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);

set( h_ax,   'Xgrid','on',...
    'Ygrid','on',...
    'Yscale','linear',...
    'FontWeight','bold',...
    'FontUnits','normalized',...
    'FontSize',norm_size);




save 'Quench_MCBXFA_ext'