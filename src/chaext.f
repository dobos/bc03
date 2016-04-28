	SUBROUTINE CHAEXT(S,E,N)

!	Changes extension (suffix) to filename s
!	SUN/FORTRAN version. G. Bruzual. 08-NOV-1988

	character*(*) s,e
	l=len(s)
	k=len(e)
	ib=largo(s)	
	id=0
	if (n >= 0) then
!		Finds last dot in array s
		do i=ib-1,1,-1
		if (s(i:i).eq.'.') then
			id=i
			if (s(i+1:i+1).eq.'/') id=ib
			goto 1
		endif
		enddo
	endif
1	if (id.eq.0) id=ib+1
	s(id:id)='.'
	if (id+k.gt.l) k=l-id
	do i=1,k
	s(id+i:id+i)=e(i:i)
	enddo
	n=id+k
	do i=id+k+1,l
	s(i:i)=' '
	enddo
	return
	end
