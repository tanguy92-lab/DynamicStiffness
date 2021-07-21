************************************************************************
* This file is part of DynamicSiffness, a Fortran library that         * 
* implements the Dynamic Stiffness Method                              *
* Copyright (C) 2017  Jean-Baptiste CASIMIR,                           *
* Quartz Laboratory - Supmeca                                          *
* 3 rue Ferand Hainaut                                                 *
* 93407 SAINT-OUEN - FRANCE                                            *      
* jean-baptiste.casimir@supmeca.fr                                     *
*                                                                      *
* This program is distributed in the hope that it will be useful,      *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* This program is free software: you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.*
************************************************************************
      
************************************************************************
*     This function computes the dynamic stiffness matrix for a        *
*     straight planar beam in a global basis. (xy) is the              *
*     plane of the beam. The first principal direction of inertia  Y  *
*     belongs to global xy-plane
*                                                                      *
*     This matrix relates displacement vector d and external force     *
*     vector f at the tips of the AB beam for a given circular         *
*     frequency w according to KW . d = f where                        *
*     d = (uA,vA,tA,uB,vB,tB)^T and f = (fxA,fyA,mA,fxB,fyB,mB)^T      *
*     uA, vA are displacements along x-axis and y-axis of tip A resp.  * 
*     tA is the rotation about z-axis of tip A                         *
*     uB, vB are displacements along x-axis and y-axis of tip B resp.  * 
*     tB is the rotation about z-axis of tip B                         *    
*                                                                      *
*     Input Args :                                                     * 
*          W : circular frequency                                      *
*          S : section area                                            *
*          IZ : quadratic moment of the section about Z-axis           *
*          L : length of the beam                                      *      
*          RHO : mass density                                          *
*          E : complex Young's modulus including structural damping    *
*          R : Radius of the element (for circular beam)               *
*          X : abscissa direction of the beam.                         *
*          X : ordinate direction of the beam.                         *
*          TE : bending theory                                         *
*                                                                      *
*     Output Args :                                                    *
*            KW : Computed dynamic stiffness matrix                    *
*                                                                      *      
*     Return value :                                                   *
*            unused logical error flag                                 *
*                                                                      *
*     Update in inputs for the use of Timoshenko's theory              *
************************************************************************

************************************************************************
*     Update for Timoshenko's theory and circular beam theory          *
* from 04/2021 to 07/2021 by Tanguy BEVANCON                           *
* tanguy.bevancon@edu.supmeca.fr                                       *
************************************************************************
      
      FUNCTION GLOBALPLANARBEAM(W,S,IZ,L,RHO,E,NU,KY,R,X,Y,TE,KW)
      IMPLICIT NONE
      LOGICAL GLOBALPLANARBEAM

*     Local Dynamic Stiffness Martrix Function                         *
      LOGICAL PLANARBEAM
      
*     Circular frequency                                               *
      DOUBLE PRECISION W
      
*     Geometrical properties of the beam                               *      
      DOUBLE PRECISION S,IZ,L,R
      
*     Material properties of the beam                                  *
      DOUBLE PRECISION RHO,NU
      COMPLEX*16 E

*     Direction of the beam                                            *
      DOUBLE PRECISION X(2),Y(2),XC,YC
      
*     Angle of Direction
      DOUBLE PRECISION ALPHA1, ALPHA2
      
*     Direction Vectors
      DOUBLE PRECISION E1(3),E2(3),V1(3),V2(3),VERIF
      
*     Bending Theory                                                   *
      INTEGER TE
      
*     Timoshenko's theory   :   04/2021                                *
      DOUBLE PRECISION KY
      
*     Global Dynamic Stiffness Matrix                                  *
      COMPLEX*16 KW(6,6)

*     Local Dynamic Stiffness Matrices                                 *
      COMPLEX*16 LKW(6,6)

*     Rotation matrix and its transpose                                *
      DOUBLE PRECISION P(6,6),PT(6,6)
      
      INTEGER I,J
      LOGICAL RET
      COMPLEX*16 MAT(6,6)
      
      DATA P/36*0/

*     Computation of the local dynamic stiffness matrix                *
*     Update inputs :   04/2021                                        *   
      RET=PLANARBEAM(W,S,IZ,L,RHO,NU,KY,R,E,TE,LKW)
     
*     Computation of the rotation matrix                               *

      IF (TE.NE.4) THEN
      
          P(1,1)=(X(2)-X(1))/L
          P(2,1)=(Y(2)-Y(1))/L
          P(1,2)=-P(2,1)
          P(2,2)=P(1,1)
          P(3,3)=1
          P(4,4)=P(1,1)
          P(5,4)=P(2,1)
          P(4,5)=P(1,2)
          P(5,5)=P(2,2)
          P(6,6)=1 
          
          ELSE
          
*     Comutation of the rotation matrix for circular beam              *
*     06/2021                                                          *

*     Calculate the coordinates of the vectors with the angle          *
          CALL GETCENTER(X(1),Y(1),X(2),Y(2),L,R,XC,YC)
          
          CALL GETANGLE(X(1),Y(1),XC,YC,R,ALPHA1)
          E2(1)=COS(ALPHA1)
          E2(2)=SIN(ALPHA1)
          E2(3)=0
          
          CALL GETANGLE(X(2),Y(2),XC,YC,R,ALPHA2)
          V2(1)=COS(ALPHA2)
          V2(2)=SIN(ALPHA2)
          V2(3)=0

          VERIF=(E2(2)*(X(2)-X(1))-E2(1)*(Y(2)-Y(1)))*R
          
          IF(VERIF.LT.0) THEN
                E2(1)=-E2(1)
                E2(2)=-E2(2)
                V2(1)=-V2(1)
                V2(2)=-V2(2)
          
          ENDIF
          
          E1(1)=E2(2)
          E1(2)=-E2(1)
          V1(1)=V2(2)
          V1(2)=-V2(1)

*    The passage matrix has to be calculated as shown below but for now*
*    we use the identity matrix                                        *          
*          P(1,1)=E1(1)
*          P(2,1)=E1(2)
*          P(1,2)=E2(1)
*          P(2,2)=E2(2)
*          P(3,3)=1
*          P(4,4)=V1(1)
*          P(5,4)=V1(2)
*          P(4,5)=V2(1)
*          P(5,5)=V2(2)
*          P(6,6)=1
          P(1,1)=1
          P(2,2)=1
          P(3,3)=1
          P(4,4)=1
          P(5,5)=1
          P(6,6)=1
          
      
      ENDIF
      
*     End of modification for circular beam                            *
      
      
*     Basis change                                                     *
      KW=MATMUL(TRANSPOSE(P),MATMUL(LKW,P))
      
      
      GLOBALPLANARBEAM=.TRUE.   
            
      END
