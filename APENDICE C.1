Cálculos de los parámetros moleculares m, σ, ε con el método de
Pedersen

Nombre de rutina: PAREOS
Argumentos de entrada:
gamV : Vector de gammas de cada pseudocomponente
Argumentos de salida:
Parame : (m, σ, ε de NC pseudocomponentes)
Código:

SUBROUTINE PAREOS(gamV,parame)
DOUBLE PRECISION parame(nc,3),MM(3) !mm es masa molar de los psuedos
DOUBLE PRECISION gamV(5),gamma,m,sigma,epsi,epsk,epsks,pm1,pm2,cn,cn2 !gamV es el
vector de las gammas, gamma del psuedo 5: parametro de ajuste densidad, pm: peso
molecular,cn: carbon number
DOUBLE PRECISION zi,zk,zks,zimi,zis,a,b !fraccion calculada con pedersen y el
producto de zi*mi de segmento,suma de zi,parametros de a y b
DOUBLE PRECISION sigmas,pms !suma de: sigma,peso molecular, gamma del pseudo5
Integer i,ic,ini,fin,pseudo(10) !contador,Inicio y fin del número de carbonos
(Carbon Number) de cada pseudocomponente, vector pseudo que indica el inicio y
fin de cada corte
ic=0
zis=0.d0
zimi=0.d0
zks=0.d0
epsks=0.d0
sigmas=0.d0
pms=0.d0
pseudo= (/12,36,37,65,66,84,85,120,121,200/)
!vector que indica el número de carbono inicial y final de cada pseudo, inicia en
C12 hasta C200
! a y b, son valores obtenidos después de minimizar el error con solver de Excel
en este caso hasta minimizar la diferencia de la masa molar C7+ y la composición
C7+.
a=-111.755285d0
b=-28.019565d0
IF (pseudo_zi.eq.0) then !Bandera que indica si se van a calcular composiciones
de los pseudocomponentes
Do i=1,9,2 !5 es el número total de pseudos, 2 el tamaño de paso e.g., toma los
valores del vector pseudo 12,37,66,etc.
ini=pseudo(i)
fin=pseudo(i+1)
ic=ic+1
gamma=gamV(ic)
pms=0.d0
zis=0.d0
zimi=0.d0
zks=0.d0
sigmas=0.d0
epsks=0.d0
Do cn=ini,fin
pm1= 14*cn-4
zi=dexp((cn-a)/b)
m=(1-gamma)*(0.0257*pm1+0.8444)+gamma*(0.0101*pm1+1.7296)
sigma= (1-gamma)*(4.047-4.8013*dlog(pm1)/pm1)+gamma*(4.6169-93.98/pm1)
epsi=(1-gamma)*dexp(5.5769-9.523/pm1)+gamma*(508-234100/pm1**(1.5))
Do cn2=ini,fin
pm2= 14*cn2-4
epsk=(1-gamma)*dexp(5.5769-9.523/pm2)+gamma*(508-234100/pm2**(1.5))
zk=dexp((cn2-a)/b)
zks=zks+(zi*zk)
epsks=epsks+(zi*zk*dsqrt(epsi*epsk))
End do
pms=pms+pm1*zi
zimi=zimi+zi*m
sigmas= sigmas+sigma*zi
zis=zis+zi
End do
m=zimi/zis
sigma=sigmas/zis
pm1=pms/zis
pm(8+ic)=pm1
epsi=epsks/zks
parame(8+ic,1)=m
parame(8+ic,2)=sigma
parame(8+ic,3)=epsi
if (pseudo_zi.EQ.0) ci(NC-5+ic)=zis !composición de cada pseudo ic hasta ic=5
End do
