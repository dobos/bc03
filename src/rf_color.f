	SUBROUTINE RF_COLOR(io,t,x,y,inw,lun,bolflux,strmass,sf,evflux,snbr,pnbr,bh,sn,wd,rm,xm,bolms,gasmass,galmass,tmlr,tdpr)

!	Array declarations
	parameter (nc=100,lb=2*nc,jid=35,its=400)
	character genfile*96,envfile*96,rfcolorfile*96
!	Negative indices allow for mandatory colors
	integer n1(-4:nc),n2(-4:nc),no(lb),ly(6),kerr(nc)
	real zp(-4:nc),col(-4:nc)
	real umag,bmag,rmag,jmag,kmag,bsun,vsun,ksun,blr1,vlr1,klr1,blr,vlr,klr
	real x(inw),y(inw),fx(lb),balm(3),gwx(jid)
        real m_pn,m_hm,m_bh,m_ns,m_wd
	save no,ly,zp,sl,mb,nv,mc,kb,nb,v0
	include 'stelib.dec'
	include 'filter.dec'
	common /kpercent/ ipcall,iplast,iphead,jc
!	Common to store fluxes used to compute indices
        integer iim(jid)
        real ffb(jid),ffr(jid),ffc(jid),ffl(jid)
        common /fluxes/ iim,ffb,ffr,ffc,ffl
        common /massrem/ m_pn,m_hm,m_bh,m_ns,m_wd
	real*8 mabo(its),mbel(its),mcut
	common /topimf/ mabo,mbel,mcut,ic
	data icall/0/,ly/6*0/,last/0/
!	Define mandatory colors
	data n1/ 15,15,15,12,14, nc*0/
	data n2/125,57,84,13,15, nc*0/
	data kerr/nc*0/

!	Check for zero flux sed
	if (t.eq.0..or.bolflux.le.0.) return
	izero=0
	do i=1,inw
	if (y(i).gt.0) then
		izero=izero+1
	endif
	enddo
	if (izero.eq.0) then
!		write (6,*) 'Exiting rf_color at t =',t,' because of zero flux'
!		return
!		To allow using .?color files with a single number of lines
		do i=1,inw
		y(i) = 1.E-33
		enddo
	endif

!	Compute galaxy mass and mass in gas (except for add_bursts io = 5)
!	if (io.ne.10) then
!		galmass=gal_mass(io,t,sf)
!	endif

!	Compute Worthey indices directly from sed plus other indices
!	Store fluxes used to compute spectral indices in binary file
	do j=1,jid
       	gwx(j)=gw_ix_sed(j,x,y,inw,0)
       	enddo
!	Write results for gwx(j)
	tl=alog10(t)
	write (lun+8,105)  tl,(gwx(j),j= 1,21)
	write (lun+9,106)  tl,(gwx(j),j=22,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
	write (lun+10,110) (tl,j,iim(j),ffb(j),ffr(j),ffc(j),ffl(j),gwx(j),j=1,jid)
105     format (f10.6,21f8.3)
106     format (f10.6,1x,6f8.3,3f10.3,4f9.3,f8.3)
110	format (f10.6,i4,i4,1p4e12.3,0pf12.3)

!106    format (f10.6,8f8.3,14f12.3)
!	write (lun+9,106) tl,(gwx(j),j=14,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
!	write (lun+8,105) tl,(gwx(j),j= 1,13)

!	Check for stelib model
	if (stelib) then
		if (iphead.eq.0) then
!			Report what we are doing
			write (6,*)
!			write (6,'(x,a,$)') 'Computing SEDs and Indices '
			write (6,'(x,a/ )') 'Computing SEDs and Indices '
			iphead=1
		endif
	endif
	if (stelib.or.inw.gt.6600) then
!		Compute indices in Lick system only for stelib or hr csp_galaxev models
		ism=0
		do j=1,jid
       		gwx(j)=gw_ix_sed_lick_system(j,x,y,inw,0,ism)
       		enddo
!		Write results for gwx(j)
		write (lun+11,105) tl,(gwx(j),j= 1,21)
       		write (lun+12,106) tl,(gwx(j),j=22,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,jid)
!		write (lun+11,105) tl,(gwx(j),j= 1,13)
!		write (lun+12,106) tl,(gwx(j),j=14,25),(gwx(j),j=30,31),(gwx(j),j=26,29),(gwx(j),j=32,34)
!		Do not compute colors for _hrs_ models
		if ((inw.eq.6700.or.inw.eq.15011).and.io.gt.0) return
	endif

!	Check if number of points has changed
	if (inw.ne.last) then
!		write (6,*) 'Resetting:',inw,last
!		ireset=1
		last=inw
	endif

	if (icall.eq.0) then
		icall=1

!		Format of file modified 07/01/2003:
!		To avoid recompiling this routine the arrays n1,n2 of
!		up to nc elements each are now read from a file.
!		Get file name from environment variable RF_COLORS_ARRAYS
        	envfile='RF_COLORS_ARRAYS'
        	call getenv(envfile,rfcolorfile)
        	close (1)
        	open (1,file=rfcolorfile,status='old',form='formatted',err=2)
		write (6,*) 'List of filters in file: ',rfcolorfile(1:largo(rfcolorfile))
!		Read one line of header
		read (1,'(a)') genfile
		write (6,*) genfile(1:largo(genfile))
		kb=0
!		Read filter pairs
		do i=1,nc
		read (1,'(2i4)',err=10) n1(i),n2(i)
		kb=kb+1
		fx(kb)=n1(i)
		kb=kb+1
		fx(kb)=n2(i)
		enddo
10		mc=i-1
		close (1)
!		Add mandatory filters
		do i=-4,0
		kb=kb+1
		fx(kb)=n1(i)
		kb=kb+1
		fx(kb)=n2(i)
		enddo
!		Sort filters in numerical order
		call sort(kb,fx)
!		Find independent filters in arrays n1 and n2, and store in no
		nb=1
		no(1)=fx(1)
		do i=2,kb
		if (fx(i).gt.fx(i-1)) then
			nb=nb+1
			no(nb)=nint(fx(i))
		endif
		enddo
!		Write filter ID
                call filterid(fid)
                write (6,'(x,a)') 'Selected filters:'
                do i=1,nb
                write (6,'(i4,'': '',i4,3x,a)') i,no(i),fid(no(i))
                enddo
		write (6,'(i3,a)') mc,' colors selected:'
                write (6,'(32i4)') (n1(i),i=1,mc)
                write (6,'(32i4)') (n2(i),i=1,mc)

!		Log of solar luminosity
		sl=33.+alog10(3.826)

!		Read filter file and compute zero points
		write (6,*)
		write (6,*) 'Computing Zero Points:'
!		Compute V magnitude Vega zero point
		v0=vega_0p_n(15)
		do i=-4,mc
		zp(i)=zerop_n(n1(i),n2(i))

!		Fill arrays with filter numbers
		do j=1,nb
		if (n1(i).eq.no(j)) n1(i)=j
		if (n2(i).eq.no(j)) n2(i)=j
		if (no(j).eq.14) mb=j
		if (no(j).eq.15) nv=j
		enddo
		enddo

!		Find in array x the position of points used to define the
!		continuum at Lyman alpha
		do i=1,inw
		if (x(i).le.1120.) ly(1)=i
		if (x(i).le.1140.) ly(2)=i
		if (x(i).le.1160.) ly(3)=i
		if (x(i).le.1280.) ly(4)=i
		if (x(i).le.1300.) ly(5)=i
		if (x(i).ge.1320.) then
			ly(6)=i
			goto 1
		endif

		enddo
	endif
1	continue

!	Report what we are doing
	if (iphead.eq.0) then
		write (6,*)
		write (6,'(x,a,$)') 'Computing SEDs and Colors. '
		iphead=1
	endif

!	Compute colors for SDSS
!	ireset=1
	call sdss_color(t,x,y,inw,lun,bolflux)

!	Compute flux through each of nb filters
	do i=1,nb
	if (kerr(i).eq.0) then
		fx(i)=filter_n(no(i),x,y,inw,0.,kerr(i))
	else
		fx(i)=0.
	endif
	enddo

!	Compute colors in Vega system
	do i=-4,mc
	if (fx(n1(i)).gt.0..and.fx(n2(i)).gt.0.) then
		col(i)=zp(i)-2.5*alog10(fx(n1(i))/fx(n2(i)))
!		write (6,*) i,n1(i),kerr(n1(i)),n2(i),kerr(n2(i))
	else
		col(i)=-99.99
	endif
	enddo

!	Compute bolometric magnitude
	bolmag=4.75-2.5*alog10(bolflux)
	bolpms=bolflux-bolms
	if (bolms.gt.0.) bolrat=bolpms/bolms

!	Compute V magnitude for a 1 Mo galaxy
!	It is -27.5 magnitudes brighter for a 1E11 Mo galaxy
	vmag=v0-2.5*alog10(fx(nv))

!	Compute U, B, R, K, and J2MASS magnitudes
	bmag=vmag+col( 0)
	umag=bmag+col(-1)
	rmag=vmag-col(-2)
	kmag=vmag-col(-3)
	jmag=vmag-col(-4)

!	Compute mass-to-visual-light ratio in solar units SUPERSEDED (Feb. 27, 2004).
!	Using a G2 V sed and the filters number 14 and 15
!	in the filter file, one derives:
!		fblue(sun) = 0.138Lo
!		fvis (sun) = 0.113Lo
!	This numbers apply only to the filters in this filter library.
!	Express blue and visual flux of model galaxy (also measured in
!	Lo) in units of the blue and visual flux of the sun:
!	fblu=fx(mb)/0.138
!	fvis=fx(nv)/0.113
!	fblu=filter(14,x,y,inw,0.)/0.138
!	fvis=filter(15,x,y,inw,0.)/0.113
!	Total mass in galaxy  = 1 Mo
!	Compute mass-to-visual-light ratio
!	blr=1./fblu
!	vlr=1./fvis
!	Compute stellar-mass-to--light ratios
!	blr=strmass/fblu
!	vlr=strmass/fvis

!	The solar absolute magnitudes for U,B,V,R,I,J,H,K were calibrated against the values
!	of Binney and Merrifield 1998, Galactic Astronomy, Table 2.1 (page 53), assuming
!	Bessell filters, and the offsets used to calibrate the entire set of filters.
!	Some values need checking - particularly those using UV filters (FOCA and Galex).
!	Taken from: http://mips.as.arizona.edu/~cnaw/sun.html

!		Filter	 B&M	 here	 difference
!		U	 5.61	 5.55	 0.06
!		B	 5.48	 5.45	 0.03
!		V	 4.83	 4.80	 0.03
!		R	 4.42	 4.46	 -0.04
!		I	 4.08	 4.11	 -0.03
!		J	 3.64	 3.67	 -0.02
!		H	 3.32	 3.33	 0.01
!		K	 3.28	 3.29	 0.01


!	Compute mass-to-visual-light ratio in solar units. Improved definition (Feb. 27, 2004).
!	Use solar absolute V and B magnitudes and (B-V)sun = 0.65
	vsun=4.80
	bsun=5.45
	ksun=3.29
!	M/L using the total mass in stars (M*), i.e. no remmnants
	blr1=strmass*10.**(0.4*(bmag-bsun))
	vlr1=strmass*10.**(0.4*(vmag-vsun))
	klr1=strmass*10.**(0.4*(kmag-ksun))
!	Modified Jan. 2011 to add mass of remmnants (rm = mBH + mNS + mWD)
	blr=(strmass+rm)*10.**(0.4*(bmag-bsun))
	vlr=(strmass+rm)*10.**(0.4*(vmag-vsun))
	klr=(strmass+rm)*10.**(0.4*(kmag-ksun))

!	Number of Lyman Continuum photons = Cly (log)
!	Flux in Lyman alpha from recombination theory
!	E(Lalpha) = 4.78E-13 * 33.1 * Nuv
!	log E = log(Nuv) -10.8 = cly -10.8
!	Number of Lyman continuum photons
	phly=clyman(x,y,inw)
	if (phly.gt.0.) then
		cly=sl+alog10(phly)
		fa=cly-10.8
	else
		cly=0.
		fa=0.
	endif

!	Number of Helium ionizing photons
	phe=chelium(x,y,inw,phe2)
	if (phe.gt.0.) then
		che=sl+alog10(phe)
	else
		che=0.
	endif
	if (phe2.gt.0.) then
		che2=sl+alog10(phe2)
	else
		che2=0.
	endif

!	Stellar continuum at Lyman alpha
	if (x(1).ge.1320.) then
		scly=0.
	else
		scly=(y(ly(1))+y(ly(2))+y(ly(3))+y(ly(4))+y(ly(5))+y(ly(6)))/6.
	endif
	if (scly.gt.0.) then
		fc=sl+alog10(scly)
	else
		fc=0.
	endif

!	Ly alpha equivalent width assuming that the continuum is the stellar continuum
	if (fa.gt.0..and.fc.gt.0.) then
		ew=10.**(fa-fc)
		ew2=phly/scly/10.**(10.8)
	else
		ew=0.
		ew2=0.
	endif

!	Compute Mg2 index
	ymg2=ymag2(x,y,inw)

!	Compute 912 A break
	b9=b912(x,y,inw)

!	Compute 4000 A break
	b4=b4000(x,y,inw)

!	Compute narrow version of D4000
        b4_n=b4000vn(x,y,inw)

!	Compute SDSS version of D4000
        b4_s=b4000_sdss(x,y,inw)

!	Compute equivalent width of Balmer lines (Hgamma, Hdelta, Hbeta)
	ewbl=ew_balmer(x,y,inw,balm)

!	Compute bolometric magnitude
        bolmag=4.75-2.5*alog10(bolflux)

!	Compute specific flux, snbr, pnbr
	evf=evflux/bolflux
!	write (6,*) t,evflux,evf,xm
	snr=snbr/bolflux
	pnr=pnbr/bolflux

!	SFR/year
!	sf=sfr(t)

!	Compute quantities requested by C. Popescu
!	call popescu(tl,sf,x,y,inw)

!	Compute flux from 1500 to 2800A
!	it=0
!	do i=1,inw
!	if (x(i).ge.1500.0.and.x(i).le.2800.0) then
!		it=it+1
!		xx(it)=x(i)
!		yy(it)=y(i)
!	endif
!	enddo
!	fuv=trapz1(x,y,it)
!	write (9,*) tl, sf,fuv

!	Write results
	write (lun+3 ,102) tl,b4,b4_n,b4_s,b9,cly,che,che2,bolmag,bolflux,snr,bh,sn,pnr,wd,rm
	write (lun+4 ,103) tl,bolmag,bmag,vmag,kmag,strmass,rm,gasmass,galmass,sf,strmass+rm,blr,vlr,klr,blr1,vlr1,klr1
	write (lun+5 ,104) tl,bolmag,evflux,evf,xm,bolrat
	write (lun+1 ,101) tl,bolmag,umag,bmag,vmag,kmag,(col(i),i=1,9)
	write (lun+2 ,101) tl,rmag,jmag,kmag,            (col(i),i=10,20)
	write (lun+14,111) tl,     kmag,                 (col(i),i=21,31),tmlr,tdpr
	write (lun+84,108) tl,vmag,kmag,                 (col(i),i=32,42),cly,blr,vlr,klr
	write (lun+85,108) tl,vmag,kmag,                 (col(i),i=43,53),cly,blr,vlr,klr
	write (lun+86,108) tl,vmag,kmag,                 (col(i),i=54,64),cly,blr,vlr,klr
	write (lun+87,109) tl,vmag,                      (col(i),i=65,80)
	write (lun+88,108) tl,vmag,kmag,                 (col(i),i=81,92),cly,blr,vlr,klr
!	write (lun+17,171) tl,mabo(ic),mbel(ic),mabo(ic)+mbel(ic),mcut
	write (lun+17,171) tl,0.0,0.0,0.0,0.0
101	format (f10.6,14f10.4,1pe13.4)
102	format (f10.6,3f10.4,f10.2,3x,3f10.4,2x,f10.4,1p7e12.4)
103	format (f10.6,4f10.4,1p6e13.4,1x,6e12.4)
104	format (f10.6,f10.4,1p3e12.4,0pf10.4)
108	format (f10.6,18f10.4)
109	format (f10.6,17f10.4)
111	format (f10.6,12f10.4,1p2e13.4)
171	format (f10.6,1p4e12.3)

	return
2	write (6,*) 'Error opening file: ',rfcolorfile(1:largo(rfcolorfile))
	stop
	end
