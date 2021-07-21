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
*     This function computes the dynamic differential matrix of the    *
*     dynamic stiffness matrix for a circular beam in a local basis.   *
*                                                                      *
*     Input Args :                                                     * 
*          W : circular frequency                                      *
*          S : section area                                            *
*          IZ : quadratic moment of the section about Z-axis           *
*          L : length of the beam                                      *      
*          RHO : mass density                                          *
*          E : complex Young's modulus including structural damping    *
*		   Nu : Poisson's ratio                                        *
*          R : Radius of the element                                   *
*          KN : Timoshenko's Section Reduction                         *
*                                                                      *
*     Output Args :                                                    *
*            DW : Computed dynamic differential matrix                 *
*                                                                      *
*     Return value :                                                   *
*            unused logical error flag                                 *
************************************************************************

      FUNCTION DYNAMICDIFFERENTIALMATRIX(W,S,IZ,L,RHO,E,NU,R,KN,DW)
      IMPLICIT NONE
      LOGICAL DYNAMICDIFFERENTIALMATRIX
      
*     Circular frequency                                               *
      DOUBLE PRECISION W
      
*     Geometrical properties of the beam                               *      
      DOUBLE PRECISION S,IZ,L,R
      
*     Timoshenko property
      DOUBLE PRECISION KN
      
*     Material properties of the beam                                  *
      DOUBLE PRECISION RHO,NU,G
      COMPLEX*16 E

*     Intermediate matrice                                             *
      COMPLEX*16 DW(6,6)
      
      G=E/2/(1+NU)
      
      DW(1,1)=0
      DW(1,2)=L/R
      DW(1,3)=0
      DW(1,4)=IZ*L/(S*R**3)
      DW(1,5)=0
      DW(1,6)=0
      DW(2,1)=-L/R
      DW(2,2)=0
      DW(2,3)=L/R
      DW(2,4)=0
      DW(2,5)=KN*E*IZ*L/(G*S*R**3)
      DW(2,6)=0
      DW(3,1)=0
      DW(3,2)=0
      DW(3,3)=0
      DW(3,4)=0
      DW(3,5)=0
      DW(3,6)=L/R
      DW(4,1)=-RHO*S*R**3*W**2*L/(E*IZ)
      DW(4,2)=0
      DW(4,3)=0
      DW(4,4)=0
      DW(4,5)=0
      DW(4,6)=0
      DW(5,1)=0
      DW(5,2)=-RHO*S*R**3*W**2*L/(E*IZ)
      DW(5,3)=0
      DW(5,4)=-L/R
      DW(5,5)=0
      DW(5,6)=0
      DW(6,1)=0
      DW(6,2)=0
      DW(6,3)=-RHO*R*W**2*L/E
      DW(6,4)=0
      DW(6,5)=-L/R
      DW(6,6)=0

      
      DYNAMICDIFFERENTIALMATRIX=.TRUE.   

      
      END
