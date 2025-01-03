Subroutine PT(gammaV,parame,zi,P0)
INTEGER(KIND=4) i,j,k,ite,Imax,inic,V_F,l
Real(DP)::
gammaV(5),error,zi(nc),P0,x(nc),y(nc),Pnueva,Tnueva,parame(nc,3),rhoph(nph)
Real(DP):: DPT,DT,DP,T0,phi(nph,nc),pmtot(nph),tiempo
Imax=40
DT= 10.D0 !DELTA TEMP. 10 KELVIN
DP=10.D5 !INCREMENTO DE PRESION EN PASCALES IGUAL A 10 BAR
call PAREOS(gammaV,parame) !obtiene parámetros moleculares de pseudos con gammaV
!se calcularon las composiciones de los pseudos y se normalizaron Ci
Do l=1,1
error=0.d0
inic=0
P=P0 !en Bar, se asigna presion inicial
if (l==1) then
IE=1
x=ci
Y=0
v_f=0
else
IE=2 !SE SELECCIONA PUNTOS DE ROCIO
y=ci
x=0
P=180.D0
dt=10.d0
imax=25
v_f=1
end if
write(5,*) T,P
write(9,'(A7,10(f6.4,2x))') 'Gamma= ',gammaV,kij(4,nc-4),kij(4,nc-3),kij(4,nc-
2),kij(4,nc-1),kij(4,nc)
CALL PRERB(x,y,inic,ite,pnueva,gammaV,rhoph,phi,pmtot) !última línea ocultada el
write(9,'(I2,2X,4(3x,f12.4))') V_F,T,P/1.D5,DPT,Error
T0=T
T=T+DT!SE CALCULA SEGUNDO PUNTO
!P=75d5 !Presión inicial supuesta
CALL PRERB(x,y,1,ite,pnueva,gammaV,rhoph,phi,pmtot)
DPT= (P/1.d5-P0)/(T-T0) !P0 ESTABA EN BAR DESDE QUE SE INGRESÓ COMO DATO INICIAL,
POR LO QUE SE CONVIERTE P EN BAR
write(9,'(I2,2X,4(3x,f12.4))') V_F,T,P/1.D5,DPT,Error
!Write(3,'(2(3x,f12.8))') DPT !Primer DeltaPT de los dos primeros puntos
write(3,'(A122)')
' T(K) P(Bar) Comp. x(i) y(i) phi_liq(
i) phi_vap(i) RhoLiq RhoVap'
Write(7,*)
' i V/F Temperatura(K) Presión(Bar) DT/DP Error '
!write(9,'(A7,3(f6.4,2x))') 'Gamma= ',gammaV
Do i=1,Imax
Error=0.d0
T0=T
P0=P
! IF (DPT.LT.0.5) THEN
V_F=0
T=T+DT!SE CALCULA PRIMER PUNTO
IF (DABS(376.05D0-T).LT.8.D0) T=376.05D0 !PARA CALCULAR LA P. DE BURBUJA A
LA T DE YACIMIENTO
CALL PRERB(x,y,1,ite,pnueva,gammaV,rhoph,phi,pmtot) !1 es para que no se
llame a la rutina de inicialización de P en este caso
! ELSE
! V_F=2
! P=(P+DP)/1.d5 !se convierte en Bar
! CALL TERB(x,y,1,ite,tnueva,gammaV,rhoph,phi,pmtot) !SE CALCULA SEGUNDO
PUNTO
!END IF
DPT= (P-P0)/1.D5/(T-T0) !Se convierten a bar las presiones para tener similitud
numérica con la temperatura.
do k = 1,nc
write(3,'(2(2x,f12.4),2x,A7,8(4x,f12.6))')
T,P/1.d5,comp(k),x(k),y(k),phi(1,k),phi(2,k),rhoph,pmtot
end do
write(3,*) ' ' !Salto de linea que divide cada punto de burbuja
Do j=1, nc
Error= error+dabs(y(j)-x(j))
End do
Write(7,'(2X,I3,2x,I2,4X,4(3x,f12.4))') i,V_F,T,P/1.D5,DPT,Error
write(9,'(I2,2X,4(3x,f12.4))') V_F,T,P/1.D5,DPT,Error
IF (ERROR.LT.0.1) THEN
write (*,*) 'LA SUMA DE LA DIFERENCIAS DE COMPOSICIONES L-V ES MENOR A 0.1'
END IF
END DO
write(7,*) 'PUNTOS DE ROCIO '!salto de linea para separar puntos de burbuja con
rocio
End Do
End Subroutine
