subroutine quad2d4n(IGAU,NGAU,GAU,WEI,maxnsd,maxnquad)
implicit none
integer I,J
integer IGAU,NGAU,maxnquad,maxnsd
real(8) ::  GAU(maxnsd,maxnquad), WEI(maxnquad)
real(8) :: GC
	  DO I = 1,maxnquad
		WEI(I) = 0.0d0
		DO J = 1,maxnsd
		  GAU(J,I) = 0.0d0
		ENDDO
	  ENDDO

      IF(IGAU.EQ.1) THEN
		NGAU = 1
		GAU(1,1) = 0.0d0
		GAU(2,1) = 0.0d0

		WEI(1) = 4.0d0
		RETURN  
	  ENDIF
	  
      IF(IGAU.EQ.2) THEN 
		NGAU = 4
		GC      =+.577350269189626d0

        GAU(1,1)=-GC
        GAU(2,1)=-GC
        GAU(1,2)=+GC
        GAU(2,2)=-GC
        GAU(1,3)=-GC
        GAU(2,3)=+GC
        GAU(1,4)=+GC
        GAU(2,4)=+GC
		WEI(1)= 1.0d0
		WEI(2)= 1.0d0
		WEI(3)= 1.0d0
		WEI(4)= 1.0d0
      ENDIF
END
 
