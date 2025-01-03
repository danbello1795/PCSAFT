C.2. Experimento PVT: Expansión a composición constante
Subroutine CCE(intervaloP,gammaV,parame,x,y,zi)
INTEGER(KIND=1) i,j,inic,ite
Real(DP)::IntervaloP(18),parame(nc,3),x(nc),y(nc),zi(nc),gammaV(5),rhoph(nph),phi
(nph,nc)
real(dp):: Vtot(18),vsat,VR
ite=20
inic=0
x=0
y=0
T=102.9+273.15d0 !en Kelvin
j=size(IntervaloP)
DO i=1,j
P=IntervaloP(i)*1d5 !se convierte a Pa, para que funcione correctamente la
rutina del flash
Call FLASH(zi,x,y,inic,ite,rhoph,gammaV)
IF (BIN.LT.0) BIN=0 !BIN ES LA V/F
VTOT(i)=BIN/RHOPH(2)+(1-BIN)/rhoph(1) !SE SUMAN LOS VOL. DE LIQ Y VAPOR, SE
CONSIDERA 1 KMOL COMO BASE DE CALCULO
IF (I.EQ.9) VSAT= VTOT(i) !7 es la posicion de la presión de burbuja a la T de
Yacimiento.
End Do
DO I=1,j
VR=Vtot(i)/vsat
write (5,'(2(3x,f8.4))') intervaloP(i),VR
END DO
End Subroutine
