************************************************************************
* This file is part of DynamicSiffness, a Fortran library that         * 
* implements the Dynamic Stiffness Method                              *
* Copyright (C) 2021  Tanguy BEVANCON,                                 * 
* Quartz Laboratory - Supmeca                                          *
* 3 rue Ferand Hainaut                                                 *
* 93407 SAINT-OUEN - FRANCE                                            *      
* tanguy.bevancon@edu.supmeca.fr                                       *
*                                                                      *
* This program is free software: you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* This program is distributed in the hope that it will be useful,      *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.*
************************************************************************
      
************************************************************************
*     This function computes the components of the dynamic stiffness   *
*     matrix for a circular beam in a local basis.                     *
*                                                                      *
*     This matrix relates traction and flexion displacement vector D   *
*     and external force vector F at the tips of the AB beam for a     *
*     given frequency w according to KW . D = F                        *
*     where T = (VA,TA,VB,TB)^T and F = (FA,MA,FB,MB)^T                *
*     VA is the vertical displacement along Y-axis of tip A            *
*     TA is the rotation about Z-axis of tip A                         *      
*     VB is the vertical displacement along Y-axis of tip B            *
*     TB is the rotation about Z-axis of tip B                         *
*     FA is the vertical external force along Y-axis applied on tip A  *      
*     FB is the vertical external force along Y-axis applied on tip B  *      
*     MA is the external moment about Z-axis applied on tip A          *
*     MB is the external moment about Z-axis applied on tip B          *     
*                                                                      *
*     Input Args :                                                     * 
*          W : circular frequency                                      *
*          S : section area                                            *
*          IZ : quadratic moment of the section about Z-axis           *
*          L : length of the beam                                      *      
*          RHO : mass density                                          *
*          E : complex Young's modulus including structural damping    *
*		   Nu : Poisson's ratio                                        *
*          kY : Timoshenko's Section Reduction                         *
*          R : Radius of the element                                   *
*                                                                      *
*     Output Args :                                                    *
*            KW : Computed dynamic stiffness matrix                    *
*                                                                      *      
*     Return value :                                                   *
*            unused logical error flag                                 *
************************************************************************

      
      FUNCTION XYCIRCULARBEAM(W,S,IZ,L,RHO,E,NU,R,KY,KW)
      IMPLICIT NONE
      LOGICAL XYCIRCULARBEAM, DYNAMICDIFFERENTIALMATRIX, RET
      
*     Circular frequency                                               *
      DOUBLE PRECISION W
      
*     Geometrical properties of the beam                               *      
      DOUBLE PRECISION S,IZ,L,R
      
*     Timoshenko property
      DOUBLE PRECISION KY
      
*     Material properties of the beam                                  *
      DOUBLE PRECISION RHO,NU,G
      COMPLEX*16 E

*     Intermediate matrices                                            *
      COMPLEX*16 EXPDW(6,6), EIGENVALUES(6,6)
      
*     Transpose Matrices                                               *
      COMPLEX*16 TW(3,3),TW2(3,3),TW2_INV(3,3),TW3(3,3),TW4(3,3)
      
*     paramaters of zgetri for Matrix 3x3                              *
      COMPLEX*16 WORK3(3)
      INTEGER IPIV(3)
      INTEGER INFO

*     paramaters of zgetri for Matrix 6x6                              *
      COMPLEX*16 WORK6(6)
      INTEGER IPIV6(6)
      
*     XY-Bending Dynamic Stiffness Matrix                              *      
      COMPLEX*16 KW(6,6)

      INTEGER          N,I,J
      PARAMETER        ( N = 6 )
      INTEGER          LDA, LDVL, LDVR
      PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )

*     .. Local Scalars ..                                              *
      INTEGER          LWORK

*     .. Local Arrays ..                                               *
*     RWORK dimension should be at least 2*N                           *
      DOUBLE PRECISION RWORK( 2*N )
      COMPLEX*16       A( LDA, N ), VL( LDVL, N ), VR( LDVR, N ),
     $                 VR_INV( LDVR, N ), WVAL( N ), WORK( LWMAX ),
     $                 DW(N,N)
     
*     .. External Subroutines ..                                       *
      EXTERNAL         ZGEEV
      
*     Build the dynamicDifferentialMatrix                              *
      RET=DYNAMICDIFFERENTIALMATRIX(W,S,IZ,L,RHO,E,NU,R,KY,DW)
      
*     Subroutine to get the eigenvalues and the eigenvectors           *
      LWORK = 780
      
*     Solve eigenproblem.

*     Function that gets the eigenvalues and eigenvectors matrixes     *
      CALL ZGEEV( 'N', 'Vectors', N, DW, LDA, WVAL, VL, LDVL,
     $            VR, LDVR, WORK, LWORK, RWORK, INFO )

*     Check for convergence.

      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
      
*     Initialization of the eigenvalues Matrix                         *
      DO I=1,6
          DO J=1,6
              EIGENVALUES(I,J)=0;    
          ENDDO
      ENDDO
*     Build the exponential matrix of eigenvalues                      *
      DO I=1,6
            EIGENVALUES(I,I)=EXP(CMPLX(WVAL(I)))
      ENDDO
      
      
*     Compute the inverse matrix of eigenvectors                       *
      VR_INV=VR
      CALL ZGETRF(6,6,VR_INV,6,IPIV6,INFO)
      CALL ZGETRI(6,VR_INV,6,IPIV6,WORK,6,INFO)
      
*     Compute the exponential matrix                                   *
      EXPDW=MATMUL(MATMUL(VR,EIGENVALUES),VR_INV)
      
*     Formulas to get the dynamic stiffness matrix                     *      
      TW=EXPDW(1:3,1:3)
      TW2=EXPDW(1:3,4:6)
      TW3=EXPDW(4:6,1:3)
      TW4=EXPDW(4:6,4:6)
      
      TW2_INV=TW2
      CALL ZGETRF(3,3,TW2_INV,3,IPIV,INFO)
      CALL ZGETRI(3,TW2_INV,3,IPIV,WORK3,3,INFO)
      
      
      KW(1:3,1:3)=MATMUL(TW2_INV,TW)
      KW(1:3,4:6)=-TW2_INV
      KW(4:6,1:3)=TW3-MATMUL(MATMUL(TW4,TW2_INV),TW)
      KW(4:6,4:6)=MATMUL(TW4,TW2_INV)
      
      
      XYCIRCULARBEAM=.TRUE.
      
      END
      
* Function to get the coordinates of the center from 2 points		   *
      SUBROUTINE GETCENTER(X1,Y1,X2,Y2,L,R,XCENTER,YCENTER)
      IMPLICIT NONE
	  DOUBLE PRECISION X1,Y1,X2,Y2,R,XMIDDLE,YMIDDLE,XCENTER,YCENTER
	  DOUBLE PRECISION XVECTOR,YVECTOR,L,L2
	 
	  XMIDDLE = (X1+X2)/2
	  YMIDDLE = (Y1+Y2)/2
	  L2 = SQRT(R**2-(L/2)**2)
* the center is founded with a vector product		     			   *
	  YVECTOR = -(X2-X1)/L
	  XVECTOR = (Y2-Y1)/L
	  XCENTER = XMIDDLE + L2*XVECTOR
	  YCENTER = YMIDDLE + L2*YVECTOR
	 
      END SUBROUTINE GETCENTER

* Function to get the angle of the arc                      		   *        
      SUBROUTINE GETTHETA(RADIUS,LENGTH,THETA)
      IMPLICIT NONE
	  DOUBLE PRECISION LENGTH,LENGTH2,RADIUS,THETA  
      
      LENGTH2=LENGTH/2
*	  Theta is the opening angle between the 2 points			       *
      THETA = 2*ATAN(LENGTH/(2*SQRT((RADIUS**2-LENGTH2**2))))
      
      END SUBROUTINE GETTHETA
  
* Function to get the angle between the center of a circle and a point *
* of his radius														   *
      SUBROUTINE GETANGLE(X,Y,XC,YC,RADIUS,ALPHA)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,XC,YC,RADIUS,ALPHA
	  DOUBLE PRECISION XVECT,YVECT,PDTSCAL,PDTVECT,PI
	
	  PI = 4.D0*DATAN(1.D0)
	  XVECT = (X-XC)/RADIUS
	  YVECT = (Y-YC)/RADIUS
* the scalar product give the cosinus								   *
	  PDTSCAL = XVECT
* the vector product give the sinus									   *
	  PDTVECT = -YVECT
      IF (PDTSCAL.GT.0) THEN
		ALPHA = ATAN(PDTVECT/PDTSCAL)
	  ELSEIF (PDTSCAL.LT.0) THEN
		ALPHA = PI + ATAN(PDTVECT/PDTSCAL)
	  ELSE
		ALPHA = ASIN(PDTVECT)
	  ENDIF
	
      END SUBROUTINE GETANGLE
