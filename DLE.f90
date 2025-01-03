C.3. Experimento PVT: Expansión de liberación diferencial
Subroutine DLE(intervaloP,gammaV,parame,x,y,zi)
INTEGER(KIND=1) i,j,k,inic,ite
Real(DP):: IntervaloP(11),parame(nc,3),x(nc),y(nc),zi(nc),gammaV(5),rhoph(nph)
Real(DP):: N,Nvap(11),Nliq(11),Vgas(11),Vliq(11),Vsto,RGA(11),Nlib !Vsto volumen
stock tank oil
ite=20
inic=0
k=size(IntervaloP)
T=102.9d0+273.15d0!en Kelvin, T de yacimiento
N=1000 !N número base de moles
nliq(1)=N
nvap(1)=0
! P= 1.01325d5
Call PAREOS(gammaV,parame)
x=0
Write(4,*) 'RESULTADOS DLE'
P=IntervaloP(1)*1d5 !Se convierte en Pa
Call FLASH(zi,x,y,inic,ite,rhoph,gammaV)
Do i=2,k
P=intervaloP(i)*1d5
Call FLASH(zi,x,y,inic,ite,rhoph,gammaV)
If (BIN.lt.0) bin=0
!write(5,*)BIN
Nvap(i)=Nliq(i-1)*bin
Nliq(i)=Nliq(i-1)-Nvap(i)
zi=x
End Do
T=15.56d0+273.15d0
Do i=1,k-1
Nlib=0
do j = i+1, k
Nlib=Nlib+Nvap(j) !Nlib moles liberados
end do
Vgas(i)=(Nlib*8.314*T/101325)*35.3147 !unidad en m3, si se quiere en pies cúbicos
multiplicar por 35.3147
End Do
Call FLASH(zi,x,y,inic,ite,rhoph,gammaV)
Vsto=Nliq(k)/1000/rhoph(1)*6.2898 !6.2898 factor de conv. de m3 a barriles
estándar
write(5,*) 'P(Bar) NVapor Nliq Vgas RGA'
Do i=1,k
RGA=Vgas(i)/Vsto
write (5,'(5(3x,f10.4))') intervalop(i),Nvap(i),Nliq(i),Vgas(i),RGA(i)
end do
!write (5,'(7(3x,f10.4))') P/1.d5,gammaV,rhoph(1)
End Subroutine
