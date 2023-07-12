
Se han guardado las curvas de Capacidad calorifica y conductividad térmica en archivos .csv fuera del código

Se llama a estos datos con un read matrix y los guarda como dos arrays

Parámetros físicos de los materiales:
	Hay que añadir fits de Jc

%%VARIABLES
Se dan tres opciones para ElemQuench:

Con un if

ElemQuench=1     %Polo
ElemQuench=Nlay  %Midplane
ElemQuench= Ne*idivide(int32(m),int32(2))+idivide(int32(Nturn),int32(2))+Nturn*idivide(int32(Nlay),int32(2)); % Centro de la bobina

He quitado la expresión

A_aux=A_v;
A_v=A_h;
A_h=A_aux;

Creo que es simplemente cambiar las matrices horizontales por las verticales