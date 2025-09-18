! Determina los autovalores y autovectores de una matriz tridiagonal simetrica
! d(:) contiene los elementos de la diagonal como entrada y como salida
!       los autovalores
! e(:) contiene los elementos no diagonales como entrada y es 
!      destruida
! z(:,:) como entrada es una matriz diagonal o la matriz de transformacion
!        si se utilizo tridiagon.f90 para reducir una matriz a tridiagonal.
!        Como salida contiene los autovectores normalizados de d(k) en la
!        columna k-esima.
! vectors: si es .false. no se calcula los autovectores y el numero de
!          operaciones se reduce de ~ 3*n^3 a ~ 30*n^2
! version modificada desde Numerical Recipes (tqli).
SUBROUTINE diagotri(d,e,N,z,vectors)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::N 
  REAL(8),INTENT(INOUT)::d(N),e(N),z(N,N)
  LOGICAL,INTENT(IN)::vectors
  integer::i,iter,k,l,m
  real(8)::b,c,dd,f,g,p,r,s,pythag,ff(N)

  e(:)=eoshift(e(:),1) !renumera elementos de e (los desplaza un sitio)
  do l=1,N
     iter=0
     ITERATE: do
        do m=l,N-1 !busca un elem no-diagonal pequeño para dividir la matriz
           dd=abs(d(m))+abs(d(m+1))
           if(abs(e(m))+dd.eq.dd)EXIT
        enddo
        if(m.eq.l) EXIT ITERATE
        if(iter.eq.30)then
           write(6,*)'demasiadas iteraciones (30)'; stop
        end if
        iter=iter+1
        g=(d(l+1)-d(l))/(2.d0*e(l))
        r=pythag(g,1.d0)
        g=d(m)-d(l)+e(l)/(g+sign(r,g))
        s=1.d0; c=1.d0; p=0.d0
        do i=m-1,l,-1  !restaura la forma tridiagonal 
           f=s*e(i); b=c*e(i)
           r=pythag(f,g); e(i+1)=r
           if(r.eq.0.d0) then   !se recupera de un underflow
              d(i+1)=d(i+1)-p
              e(m)=0.d0
              CYCLE ITERATE
           endif
           s=f/r; c=g/r; g=d(i+1)-p
           r=(d(i)-g)*s+2.d0*c*b
           p=s*r; d(i+1)=g+p; g=c*r-b
           if(vectors)then  !calcula autovectores 
              ff(:)=z(:,i+1)
              z(:,i+1)=s*z(:,i)+c*ff(:)
              z(:,i)=c*z(:,i)-s*ff(:)
           end if
        end do
        d(l)=d(l)-p; e(l)=g; e(m)=0.d0
     enddo ITERATE
  enddo

END SUBROUTINE diagotri
!*****************************************************
REAL(8) FUNCTION pythag(a,b)
  !calcula sqrt(a**2+b**2) sin problemas de underflow y overflow
  REAL(8)::a,b
  REAL(8)::absa,absb
  absa=abs(a)
  absb=abs(b)
  if(absa.gt.absb)then
     pythag=absa*sqrt(1.+(absb/absa)**2)
  else
     if(absb.eq.0.d0)then
        pythag=0.d0
     else
        pythag=absb*sqrt(1.d0+(absa/absb)**2)
     end if
  end if
END FUNCTION pythag
