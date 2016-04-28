	SUBROUTINE GAUSSIAN_V_DISP(X,Y,N,SIGMA)

c	Applies a gaussian filter to the sed in (X,Y) to reproduce
c	the effect of the velocity dispersion of stars in a galaxy
c	Enter sigma in Km/s

c	Array declaration
	real x(n),y(n),z(50000),u(10000),g(10000)

c	If sigma = 0, return
	if (sigma.le.0.) return

c	Speed of light in Km/s
	c=3.00E5

c	Number of sigmas to deviate from v=0
	m=6

c	Copy array y to array z
	if (n.gt.50000) then
		write (6,*) 'Too many data points in GAUSSIAN_V_DISP:',n
		stop
	endif
	do i=1,n
	z(i)=y(i)
	enddo

	do i=1,n
	xmax=c*x(i)/(c-m*sigma)
	call locate(x,n,xmax,i1)
	m2=i1+1
	m1=2*i-m2
	if (m1.lt.1) m1=1
	if (m2.gt.n) m2=n
	if (m2-m1+1.gt.10000) then
		write (6,*) 'Too many data points in GAUSSIAN_V_DISP:',m1,m2,m2-m1+1
		stop
	endif
	k=0
	do j=m2,m1,-1
	v=(x(i)/x(j)-1.)*c
	k=k+1
	u(k)=v
	g(k)=z(j)*gauss(v,0.,sigma)
	enddo
	y(i)=trapz1(u,g,k)
	enddo
	return
	end
