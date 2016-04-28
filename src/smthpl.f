	SUBROUTINE SMTHPL(Y,NP,NA)

c	SMOOTHS OUT ARRAY Y BY A MEAN RUNNING AVERAGE OF NA POINTS.
c	TAKES CARE OF BORDER EFFECTS.
c	NP IS THE TOTAL NUMBER OF POINTS IN THE ARRAY. WORKS FOR UP TO
c	200000 POINTS.
c	RETURNS SMOOTHED ARRAY INTO SAME ARRAY.

	REAL Y(NP),BUF(200000)
	NB=NA/2
	DO I=1,NB
	BUF(I)=AVERAG(Y,I,I-1)
	M=NP-I+1
	BUF(M)=AVERAG(Y,M,I-1)
	ENDDO
	I1=NB+1
	I2=NP-NB
	DO I=I1,I2
	BUF(I)=AVERAG(Y,I,NB)
	ENDDO
	DO I=1,NP
	Y(I)=BUF(I)
	ENDDO
	RETURN
	END

	FUNCTION AVERAG(Y,I,J)
c	TAKES AVERAGE OF (2*J+1) POINTS OF ARRAY Y CENTERED ON Y(I)
	REAL Y(20000)
	AVERAG=0.
	L1=I-J
	L2=I+J
	DO L=L1,L2
	AVERAG=AVERAG+Y(L)
	ENDDO
	AVERAG=AVERAG/(2*J+1.)
	RETURN
	END