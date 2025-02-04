!*==sp_AES.f90    
subroutine spAES(iprint,tsst,msst,mssq,tauq,mezz,mezj,gcor,fcor,ecor,szcor,kapcor,mm05cor,nkpcor,ikmcor,izero,itxray,bcor,bcors,     &
             & ncstmax)
!   ********************************************************************
!   *       calculation of     aes  -  spectra                         *
!   ********************************************************************
!   * run ssite always for the  nl+1  angular momentum expansion       *
!   ********************************************************************
!
   use mod_energy , only:etab , nemax , efermi
   use mod_rmesh , only:nrmax , jrws , r2drdi
   use mod_calcmode , only:orbpol
   use mod_angmom , only:nkm , nmemax , nkmmax , nlmax , cgc , nlinq , nkmq , nlq
   use mod_sites , only:nq , iqat , nqmax
   use mod_files , only:datset , system , ldatset , ifilbuildbot , wrbuildbot , taudir
   use mod_types , only:soctl , ntmax , lcxray , ncxray , imt , nlt , nt , txt_t , ctl
   use mod_constants , only:ry_ev , c_au
   use s_amecoul
   use s_cinit
   use s_core
   use s_input_find_section
   use s_mecoul
   use s_rmatstruct
   use s_rradint
   use s_section_find_keyword
   use s_section_get_core_level_info
   use s_section_set_integer
   use s_ssite
   use s_stop_message
   use s_tau_dir_inquire
   use s_tau_dir_read
   use s_wrhead
   use iso_fortran_env
   implicit none
!
! PARAMETER definitions rewritten by SPAG
!
   character(40) , parameter :: ROUTINE = 'aes'
   integer , parameter :: IFILCBWF0 = 60
!
! Dummy argument declarations rewritten by SPAG
!
   integer :: NCSTMAX
   integer :: IPRINT
   complex(REAL64) , dimension(nkmmax,nkmmax,ntmax) :: TSST
   complex(REAL64) , dimension(nkmmax,nkmmax,ntmax) :: MSST
   complex(REAL64) , dimension(nkmmax,nkmmax,nqmax) :: MSSQ
   complex(REAL64) , dimension(nkmmax,nkmmax,nqmax) :: TAUQ
   complex(REAL64) , dimension(nkmmax,nkmmax,ntmax,nmemax) :: MEZZ
   complex(REAL64) , dimension(nkmmax,nkmmax,ntmax,nmemax) :: MEZJ
   real(REAL64) , dimension(nrmax,2,ncstmax) :: GCOR
   real(REAL64) , dimension(nrmax,2,ncstmax) :: FCOR
   real(REAL64) , dimension(ncstmax) :: ECOR
   real(REAL64) , dimension(ncstmax) :: SZCOR
   integer , dimension(ncstmax) :: KAPCOR
   integer , dimension(ncstmax) :: MM05COR
   integer , dimension(ncstmax) :: NKPCOR
   integer , dimension(ncstmax,2) :: IKMCOR
   integer , dimension(ncstmax) :: IZERO
   integer :: ITXRAY
   real(REAL64) , dimension(ntmax) :: BCOR
   real(REAL64) , dimension(ntmax) :: BCORS
!
! Local variable declarations rewritten by SPAG
!
   real(REAL64) , allocatable , dimension(:,:,:) :: amecif , amecig
   logical :: calcint , getirrsol , printame , printmele , readmele
   character(2) :: cl
   real(REAL64) :: de , deme , ee , ei , epsd , imefin , jmc , mj , p3 , q4 , reme , rj , sk , wa , wb , wc , wd , we , wf , wgte
   complex(REAL64) , allocatable , dimension(:,:,:) :: difme , iaes , mmed , mmee , mmetildd , mmetilde
   complex(REAL64) :: ea , eb , ebot , ec , ed , p , rme
   complex(REAL64) , allocatable , dimension(:) :: eme
   complex(REAL64) , save :: eryd
   complex(REAL64) , allocatable , dimension(:,:) :: etabfin , taub , tauc , tautleed
   character(80) :: filnam , spec
   integer :: i , ib , ic , icst , icstme , id , ie , ieb , iec , ied , ieme , ieme30 , ieme40 , ifil , ifilb , ifilc , ifild ,    &
            & ifilfin , ifilme , ifilval , iflag_tau_dir_inquire , ii , ikmb , ikmc , ikmd , il , ilam , ilamp , im , iq , is ,    &
            & it , itx , j , k , kap , l , lam , lam1 , lam2 , lam3 , lamp , lamp1 , lamp2 , lamp3 , mjm05 , ms , n , ncst , ne ,  &
            & nebot , nefin , neint , neme , nemeinp , nememax , nepanel_taudir , nepath_taudir , nkmb , nkmc , nkmd , nl , nlm ,  &
            & ntxrsgrp
   integer , save :: ia_err , iepanel , iepath
   integer , dimension(2) :: igrid_taudir , netab_taudir
   integer , allocatable , dimension(:) :: lammagcpl
   complex(REAL64) , allocatable , dimension(:,:,:,:,:) :: me
   real(REAL64) , allocatable , dimension(:) :: rint
   character(1) , dimension(5) , save :: shell
   complex(REAL64) , dimension(nkmmax,nkmmax,ntmax) :: ssst
   character(4) :: strfile
   character(4) , dimension(12) , save :: strme
   character(3) , dimension(0:4) , save :: subsh , subshp
   complex(REAL64) , dimension(2) :: waes , wtmp
   real(REAL64) , dimension(2) :: xnorm
!
! End of declarations rewritten by SPAG
!
!
!*** start of declarations rewritten by spag
! parameter definitions
! dummy arguments
! local variables
!*** end of declarations rewritten by spag
!
   integer :: spag_nextblock_1
!
   data strme/'me1' , 'me2' , 'me3' , 'me4' , 'me5' , 'me6' , 'me7' , 'me8' , 'me9' , 'me10' , 'me11' , 'me12'/ , ia_err/0/
   data shell/'k' , 'l' , 'm' , 'n' , 'o'/
   data subsh/'1  ' , '2,3' , '4,5' , '6,7' , '8,9'/
   data subshp/'1  ' , '23 ' , '45 ' , '67 ' , '89 '/
   data eryd/(999999D0,999999D0)/
   data iepanel/1/ , iepath/1/
!
   allocate (iaes(ncstmax,2,(2*nemax-1)),rint(nrmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: iaes')
!
   allocate (mmed(nkmmax,nkmmax,nkmmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: mmed')
!
   allocate (mmee(nkmmax,nkmmax,nkmmax),lammagcpl(nkmmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: mmee')
!
   allocate (taub(nkmmax,nkmmax),tauc(nkmmax,nkmmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: tauc')
!
   allocate (mmetildd(nkmmax,nkmmax,nkmmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: mmetildd')
!
   allocate (mmetilde(nkmmax,nkmmax,nkmmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: mmetilde')
!
   allocate (tautleed(nkmmax,nkmmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: tautleed')
!
   allocate (difme(nkmmax,nkmmax,nkmmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: difme')
!
   allocate (amecif(nkmmax,nkmmax,2*nlmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: amecif')
!
   allocate (amecig(nkmmax,nkmmax,2*nlmax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: amecig')
!
   allocate (etabfin((2*nemax-1),2),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: ecortab')
!
!==================================================================start
   write (6,99002)
99002 format (' ',//,10x,62('*'),/,10x,'*',60x,'*',/,10x,'*        ****   *****            ***    ******   ****        *',/,10x,   &
             &'*       *    *  *    *          *   *   *       *    *       *',/,10x,                                              &
             &'*       *       *    *         *     *  *       *            *',/,10x,                                              &
             &'*        ****   *****    ***   *******  *****    ****        *',/,10x,                                              &
             &'*            *  *              *     *  *            *       *',/,10x,                                              &
             &'*       *    *  *              *     *  *       *    *       *',/,10x,                                              &
             &'*        ****   *              *     *  ******   ****        *',/,10x,'*',60x,'*',/,10x,62('*'),//)
!-----------------------------------------------------reading from input
!
   call input_find_section('task',1)
!
   call section_set_integer('it',itxray,1,0)
   it = itxray
!
   call section_set_integer('neme',nemeinp,10,0)
   nememax = nemeinp
!
   allocate (me(nkmmax,nkmmax,nkmmax,nememax,nememax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: me')
!
   allocate (eme(nememax),stat=ia_err)
   if ( ia_err/=0 ) call stop_message(routine,'alloc: me')
!
   call section_find_keyword('printme',printmele)
!
   call section_find_keyword('printame',printame)
!
   call section_get_core_level_info(cl,ncxray(it),lcxray(it))
!
   iq = iqat(1,it)
   im = imt(it)
!
   if ( ldatset/=0 ) then
      filnam = trim(datset)//trim(txt_t(it))//'.'
   else
      filnam = trim(txt_t(it))//'.'
   endif
!
   filnam = trim(filnam)//shell(ncxray(it))
   spec = '  '//shell(ncxray(it))
!
   if ( ncxray(it)/=1 ) then
      spec = trim(spec)//subsh(lcxray(it))
      filnam = trim(filnam)//subshp(lcxray(it))
   endif
!
   filnam = trim(filnam)//'.aes'
   spec = trim(spec)//' - aes spectrum of '//trim(txt_t(it))//' in  '//trim(system)
!
   open (unit=7,file=trim(filnam))
   write (6,'(10x,a,a,/)') 'spec-file :  ( 7) ' , filnam
   write (6,'(a)') spec
   write (6,99003) it , ncxray(it) , lcxray(it)
!----------------------------------------------------------------formats
99003 format (/,4x,' core quantum-numbers  for  it=',i2,':   n=',i2,'  l=',i2,/)
!
! ----------------------------- increase angular momentum expansion by 1
! --------------------- to deal with final states with l_fin = l_ini + 1
   nl = 0
   do iq = 1 , nq
      nlq(iq) = nlq(iq) + 1
      nl = max(nl,nlq(iq))
      if ( nlq(iq)>nlmax ) call stop_message(routine,'nlq(iq)>>nlmax')
      nkmq(iq) = 2*nlq(iq)**2
      nlinq(iq) = 2*nlq(iq)*(2*nlq(iq)-1)
   enddo
   do itx = 1 , nt
      nlt(itx) = nlq(iqat(1,itx))
      ctl(it,nlt(itx)) = ctl(it,nlt(itx)-1)
      soctl(it,nlt(itx)) = soctl(it,nlt(itx)-1)
   enddo
   lam = 0
   do k = 1 , 2*nlmax - 1
      if ( mod(k,2)==0 ) then
         l = k/2
         kap = +l
      else
         l = (k-1)/2
         kap = -l - 1
      endif
      sk = dble(sign(1,kap))
      rj = dble(l) - sk/2D0
      do mjm05 = nint(-rj-0.5D0) , nint(rj-0.5D0)
         mj = dble(mjm05) + 0.5D0
         lam = lam + 1
         jmc = dble(l) + sk/2D0
         if ( abs(mj)<l ) then
            lammagcpl(lam) = nint(2*l*(jmc+0.5D0)+jmc+mj+1)
         else
            lammagcpl(lam) = 0
         endif
      enddo
   enddo
   nlm = nl**2
   nkm = 2*nlm
   nkmb = 2*(nl-1)**2
   nkmc = nkmb
   nkmd = nkm
! ======================================================================
! read tau to get vb-info
! ======================================================================
   ifilval = 99
   open (unit=ifilval,status='scratch',form='unformatted',access='direct',recl=(nkmmax*nkmmax*2*8))
!
   itx = 0
!
   call tau_dir_inquire(taudir,nepanel_taudir,nepath_taudir,igrid_taudir,netab_taudir,iflag_tau_dir_inquire)
!
   ne = netab_taudir(1)
!
   if ( nq>1 ) call stop_message(routine,'nq > 1')
   iq = 1
   do ie = 1 , ne
      eryd = etab(ie,1)
      call tau_dir_read(taudir,eryd,iepanel,iepath,ie,mssq,tauq)
      etab(ie,1) = eryd
      write (ifilval,rec=ie) tauq(:,:,iq)
   enddo
!
   ebot = etab(1,1)
!
!==================================================================
   call wrhead(7,filnam,'sp-aes    ',ne)
!
   ntxrsgrp = 1
   write (7,99012) 'ntxrsgrp  ' , ntxrsgrp
   write (7,99013) 'spectrum  ' , trim(spec)
99013 format (a10,a)
   write (7,99012) 'it        ' , it
   write (7,99012) 'ncxray    ' , ncxray(it)
   write (7,99012) 'lcxray    ' , lcxray(it)
! ======================================================================
! prepare core states
! corecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecoreco
   ncst = 4*lcxray(it) + 2
   if ( ncst>ncstmax ) then
      write (6,99001) lcxray(it) , ncstmax
!
99001 format (//,60('*'),/,10x,' stop in <aes>',/,10x,'lcxray=',i3,' too large for ncstmax=',i3)
      call stop_message(routine,' ')
   endif
!
   call core(iprint,gcor,fcor,ecor,szcor,kapcor,mm05cor,nkpcor,ikmcor,izero,itxray,bcor,bcors,ncstmax)
!
   do ifil = 6 , 7
      write (ifil,99006) ncst , (nkpcor(icst),icst=1,ncst)
99006 format (//,' core states :',//,' ncst:  ',i4,/,' nkpcor:',20I4)
      write (ifil,99007)
99007 format (/,' icst  n   l  kap  mue  ikm     norm   ','     e(ry)       e(ev)    <sigma_z>  i0')
   enddo
   do icst = 1 , ncst
!
      do k = 1 , nkpcor(icst)
         do n = 1 , jrws(im)
            rint(n) = r2drdi(n,im)*(gcor(n,k,icst)**2+fcor(n,k,icst)**2)
         enddo
         call rradint(im,rint,xnorm(k))
      enddo
!
      do ifil = 6 , 7
         write (ifil,99008) icst , ncxray(it) , lcxray(it) , kapcor(icst) , (2*mm05cor(icst)+1) , ikmcor(icst,1) , xnorm(1) ,      &
                          & ecor(icst) , ecor(icst)*ry_ev , szcor(icst) , izero(icst)
99008    format (5I4,'/2',i4,f12.6,f12.4,f12.3,f12.4,i5)
         if ( nkpcor(icst)==2 ) write (ifil,99009) ikmcor(icst,2) , xnorm(2)
99009    format (22x,i4,f12.6)
      enddo
!
   enddo
!   corecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecore
   write (6,99004) ne
99004 format (' ',//,10x,'number of tabulated energies',i5,/)
   write (7,99005) ne , 1
99005 format (/,' number of energies',i5,/,' output format ifmt',i5,/)
!  =====================================================================
!         set up energy table for matrix elements
!  =====================================================================
   neme = min(nemeinp,ne)
   deme = (efermi-dreal(ebot))/dble(neme-1)
   calcint = .false.
   getirrsol = .false.
   do ie = 1 , neme
      eme(ie) = ebot + deme*dble(ie-1)
   enddo
!  =====================================================================
!                    opening of files
!  =====================================================================
   ifilfin = 98
   ifild = ifilcbwf0
   ifilb = ifilcbwf0 + 1
   ifilc = ifilcbwf0 + 2
   open (unit=ifilfin,status='scratch',form='unformatted',access='direct',recl=(4+2*(2+4*nrmax))*2*8)
   open (unit=ifild,status='scratch',form='unformatted',access='direct',recl=(4+2*(2+4*nrmax))*2*8)
   open (unit=ifilb,status='scratch',form='unformatted',access='direct',recl=(4+2*(2+4*nrmax))*2*8)
   open (unit=ifilc,status='scratch',form='unformatted',access='direct',recl=(4+2*(2+4*nrmax))*2*8)
!
!  =====================================================================
!
   de = dreal(etab(2,1)) - dreal(etab(1,1))
   write (6,99014) ne , dreal(etab(1,1)) , dreal(etab(1,1))*ry_ev , dreal(etab(ne,1)) , dreal(etab(ne,1))*ry_ev , de , de*ry_ev
99014 format (' tau:  number of energies read in  ne=',i3,/,'       from   e(1)  =',f8.3,' ry    (',f8.3,' ev)',/,                 &
             &'       to     e(ne) =',f8.3,' ry    (',f8.3,' ev)',/,'       with   de    =',f8.3,' ry    (',f8.3,' ev)')
!
!=======================================================================
!               calculate angular matrix elements
!  angularangularangularangularangularangularangularangularangularangula
!
   write (6,*) 'angular matrix elements will be calculated'
!
   call amecoul(amecig,amecif,iprint,nlmax,nkmmax)
!
   if ( printame ) then
      do il = 1 , 2*nlmax
         write (6,*) 'angular matrix elements for lr=' , il
         call rmatstruct('amecig',amecig(1,1,il),nkmd,nkmmax,0,0,1,1D-8,6)
         call rmatstruct('amecif',amecif(1,1,il),nkmd,nkmmax,0,0,1,1D-8,6)
      enddo
   endif
!  angularangularangularangularangularangularangularangularangularangula
!
! ======================================================================
!               calculate spr-aes-spectrum
! ======================================================================
   write (6,*)
   write (6,*) 'calculating spr-aes spectrum'
   write (6,*)
   call cinit(ncstmax*2*(2*nemax-1),iaes)
!corecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecor
!
   do icst = 1 , ncst
      spag_nextblock_1 = 1
      spag_dispatchloop_1: do
         select case (spag_nextblock_1)
         case (1)
            write (6,*) '       core state ' , icst
            write (6,*) '       -------------'
            call cinit(nkmmax*nkmmax*nkmmax*nememax*nememax,me)
            call cinit(nkmmax*nkmmax*nkmmax,mmed)
            call cinit(nkmmax*nkmmax*nkmmax,mmee)
!     reading me from file if...
            ifilme = 23
            inquire (file='me',exist=readmele)
            if ( readmele .eqv. .true. ) then
               strfile = 'me'
               spag_nextblock_1 = 2
               cycle spag_dispatchloop_1
            endif
            strfile = strme(icst)
            inquire (file=trim(strfile),exist=readmele)
            spag_nextblock_1 = 2
         case (2)
            if ( readmele ) printmele = .not.readmele
            if ( readmele ) open (unit=ifilme,file=trim(strfile))
            if ( readmele ) then
               write (6,*)
               write (6,*) '     reading meles from file !!!!!!!!!!! '
               do i = 1 , 6
                  read (ifilme,*)
               enddo
               iec = 0
               ieb = 0
 5             iec = iec + 1
               if ( iec>neme ) goto 20
 10            ieb = ieb + 1
               do i = 1 , 2
                  read (ifilme,*)
               enddo
               do j = 1 , 1000000
                  read (ifilme,99011,err=15) icstme , ib , ic , id , rme
                  if ( icstme==icst ) me(ib,ic,id,ieb,iec) = rme
               enddo
 15            if ( iec<=neme ) then
                  if ( ieb<neme ) goto 10
                  ieb = 0
                  goto 5
               endif
 20            close (ifilme)
!
               write (6,99015) neme
99015          format (' me: wave functions read in for neme=',i3,' energies')
!
               write (6,99016) dreal(ebot) , dreal(ebot)*ry_ev , efermi , efermi*ry_ev , deme , deme*ry_ev
99016          format (' me:   from   ebot  =',f8.3,' ry    (',f8.3,' ev)',/,'       to     efermi=',f8.3,' ry    (',f8.3,' ev)',/,&
                      &'       with   deme  =',f8.3,' ry    (',f8.3,' ev)')
               close (ifilme)
            else
!=======================================================================
!                calculation of matrix elements
!mememmememecmememmemememmememememememememmememememmememememmememememmem
               write (6,*)
               write (6,*) ' calculating me for (iea,ieb,iec,ied)'
!
               do iec = 1 , neme
                  ec = ebot + deme*(iec-1)
                  call ssite(1,1,ifilc,calcint,getirrsol,ec,p,iprint,nkm,tsst,msst,ssst,mezz,mezj,orbpol)
!
                  do ieb = 1 , neme
                     eb = ebot + deme*(ieb-1)
                     ea = ecor(icst)
                     ee = dreal(eb+ec-ea)
                     ei = dimag(eb)
                     ed = dcmplx(ee,ei)
                     if ( iprint>=0 ) write (6,99010) icst , ieb , iec , iec + ieb - 1
99010                format (28x,'(',3(i3,','),i3,')')
!
                     call ssite(1,1,ifilb,calcint,getirrsol,eb,p,iprint,nkm,tsst,msst,ssst,mezz,mezj,orbpol)
!
                     call ssite(1,1,ifild,calcint,getirrsol,ed,p,iprint,nkm,tsst,msst,ssst,mezz,mezj,orbpol)
!
                     call mecoul(icst,gcor,fcor,mm05cor,nkpcor,ikmcor,lcxray,1,nkmb,ifilb,1,nkmc,ifilc,1,nkmd,ifild,it,me,ieb,iec, &
                               & amecig,amecif,nememax,ncstmax)
                  enddo
               enddo
!     mememmemememmememememememememmememememmememememmememememmeme
!     writing of matrix elements
!     mememmemememmememememememememmememememmememememmememememmeme
               if ( printmele ) then
                  ifilme = 24
                  strfile = strme(icst)
                  open (unit=ifilme,file=trim(strfile))
                  write (6,*)
                  write (6,*) ' matrix elements m(c,v,v,l) output in file me'
                  write (ifilme,'(a)') '    coulomb matrix elements m(c,v,v,l)'
                  write (ifilme,'(10x,a,i2)') 'k_c   = 1,' , ncst
                  write (ifilme,'(10x,a,i3)') 'k_v   = 1,' , nkmb
                  write (ifilme,'(10x,a,i3)') 'k_l   = 1,' , nkmd
                  write (ifilme,*)
                  write (ifilme,'(2a)') '*****************************************' , '*****************************'
                  do iec = 1 , neme
                     do ieb = 1 , neme
                        ec = ebot + deme*dble(iec-1)
                        eb = ebot + deme*dble(ieb-1)
                        write (ifilme,'(a6,f10.6,a6,f10.6)') ' eb = ' , dreal(eb) , ' ec = ' , dreal(ec)
                        write (ifilme,'(6x,a17)') 'k_c k_vb k_vc k_l'
                        do id = 1 , nkmd
                           do ic = 1 , nkmc
                              do ib = 1 , nkmb
                                 if ( cdabs(me(ib,ic,id,ieb,iec))>1.D-14 ) write (ifilme,99011) icst , ib , ic , id ,              &
                                    & me(ib,ic,id,ieb,iec)
                              enddo
                           enddo
                        enddo
                        write (ifilme,'(2a)') '*****************************************' , '*****************************'
                     enddo
                  enddo
                  close (ifilme)
               endif
            endif
!     mememmemememmememememememememmememememmememememmememememmeme
!         if ( (2*ne-1).gt.nemax ) then
!            write (6,*) ' < aes > : increase nemax to ',2*ne - 1
!            stop_regular(routine,' ')
!         end if
            imefin = dimag(etab(1,1))
            nefin = 2*ne - 1
            do ie = nefin , 1 , -1
               etabfin(ie,1) = dcmplx(2*efermi-ecor(icst),imefin) - dcmplx((nefin-ie)*de,0.D0)
            enddo
!     iebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebie
            do ied = 1 , 2*ne - 1
               write (6,*) '  ied=' , ied , '  ed=' , dreal(etabfin(ied,1))
               if ( ied<=ne ) then
                  nebot = 1
                  neint = ied
               elseif ( ied>ne ) then
                  nebot = ied + 1 - ne
                  neint = ne
               endif
               ed = etabfin(ied,1)
!     iecieciecieciecieciecieciecieciecieciecieciecieciecieciecieciecie
               do ieb = nebot , neint
                  iec = ied + 1 - ieb
!
                  if ( (ieb==nebot) .or. (ieb==neint) ) then
                     wgte = de
                  else
                     wgte = 2D0*de
                  endif
                  eb = ebot + dble(ieb-1)*de
                  ec = ebot + dble(iec-1)*de
!-----------------interpolate the calculated meles to the vb-mesh
!     interpolationinterpolationinterpolationinterpolationinterpolation
                  ieme30 = neme
                  ieme40 = neme
                  do ieme = 1 , neme
                     reme = dreal(eme(ieme))
                     if ( (ieme30==neme) .and. (reme>dreal(eb)) ) ieme30 = max(2,ieme-1)
                     if ( (ieme40==neme) .and. (reme>dreal(ec)) ) ieme40 = max(2,ieme-1)
                  enddo
                  ieme30 = min(ieme30,neme-1)
                  ieme40 = min(ieme40,neme-1)
                  p3 = (dreal(eb)-dreal(eme(ieme30)))/deme
                  q4 = (dreal(ec)-dreal(eme(ieme40)))/deme
!
                  wa = q4*(q4-1D0)/2D0
                  wb = p3*(p3-1D0)/2D0
                  wc = 1D0 + p3*q4 - p3**2 - q4**2
                  wd = p3*(p3-2D0*q4+1D0)/2D0
                  we = q4*(q4-2D0*p3+1D0)/2D0
                  wf = p3*q4
                  do ikmb = 1 , nkmb
                     do ikmc = 1 , nkmc
                        do ikmd = 1 , nkmd
                           mmed(ikmb,ikmc,ikmd) = me(ikmb,ikmc,ikmd,ieme30+0,ieme40-1)*wa + me(ikmb,ikmc,ikmd,ieme30-1,ieme40+0)   &
                            & *wb + me(ikmb,ikmc,ikmd,ieme30+0,ieme40+0)*wc + me(ikmb,ikmc,ikmd,ieme30+1,ieme40+0)                 &
                            & *wd + me(ikmb,ikmc,ikmd,ieme30+0,ieme40+1)*we + me(ikmb,ikmc,ikmd,ieme30+1,ieme40+1)*wf
                           mmee(ikmb,ikmc,ikmd) = me(ikmb,ikmc,ikmd,ieme40+0,ieme30-1)*wa + me(ikmb,ikmc,ikmd,ieme40-1,ieme30+0)   &
                            & *wb + me(ikmb,ikmc,ikmd,ieme40+0,ieme30+0)*wc + me(ikmb,ikmc,ikmd,ieme40+1,ieme30+0)                 &
                            & *wd + me(ikmb,ikmc,ikmd,ieme40+0,ieme30+1)*we + me(ikmb,ikmc,ikmd,ieme40+1,ieme30+1)*wf
                        enddo
                     enddo
                  enddo
!     interpolationinterpolationinterpolationinterpolationinterpolation
!ttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
!                   get tau for final states
!                single-scatterer approximation
!     tautautautautautautautautautautautautautautautautautautautautautau
                  eryd = etabfin(ied,1)
                  call cinit(nkmmax*nkmmax*ntmax,tsst)
                  call ssite(1,1,ifilfin,calcint,getirrsol,eryd,p,iprint,nkm,tsst,msst,ssst,mezz,mezj,orbpol)
!
                  call cinit(nkmmax*nkmmax,tautleed)
                  do i = 1 , nkmd
                     do ii = 1 , nkmd
                        tautleed(i,ii) = tsst(i,ii,it)
                     enddo
                  enddo
!     tautautautautautautautautautautautautautautautautautautautautautau
                  read (ifilval,rec=ieb) taub
                  read (ifilval,rec=iec) tauc
!
!-------tild
                  call cinit(nkmmax*nkmmax*nkmmax,mmetildd)
                  call cinit(nkmmax*nkmmax*nkmmax,mmetilde)
                  do ikmb = 1 , nkmb
                     do ikmc = 1 , nkmc
                        do ilam = 1 , nkmd
                           do ilamp = 1 , nkmd
                              mmetildd(ikmb,ikmc,ilam) = mmetildd(ikmb,ikmc,ilam) + dconjg(tautleed(ilamp,ilam))                   &
                               & *mmed(ikmb,ikmc,ilamp)
                              mmetilde(ikmb,ikmc,ilam) = mmetilde(ikmb,ikmc,ilam) + dconjg(tautleed(ilamp,ilam))                   &
                               & *mmee(ikmb,ikmc,ilamp)
                           enddo
                        enddo
                     enddo
                  enddo
!-------tild
                  call cinit(nkmmax*nkmmax*nkmmax,difme)
                  do lam1 = 1 , nkmb
                     do lam2 = 1 , nkmc
                        do lam3 = 1 , nkmd
                           difme(lam1,lam2,lam3) = (mmetildd(lam1,lam2,lam3)-mmetilde(lam2,lam1,lam3))
                        enddo
                     enddo
                  enddo
                  call cinit(2,wtmp)
                  do lamp2 = 1 , nkmb
                     do lamp3 = 1 , nkmb
                        do lamp = 1 , nkmc
                           do lamp1 = 1 , nkmc
                              call cinit(2,waes)
!     sumation over lam(leed)
                              do lam1 = 1 , nkmd
                                 do lam2 = 1 , nkmd
!
                                    if ( (lam1==lam2) .or. (lam2==lammagcpl(lam1)) ) then
                                       waes(1) = waes(1) + cgc(lam1,1)*cgc(lam2,1)*difme(lamp2,lamp,lam1)                          &
                                        & *dconjg(difme(lamp3,lamp1,lam2))
                                       waes(2) = waes(2) + cgc(lam1,2)*cgc(lam2,2)*difme(lamp2,lamp,lam1)                          &
                                        & *dconjg(difme(lamp3,lamp1,lam2))
!
!     sumation over lam(leed)
                                    endif
                                 enddo
                              enddo
                              do ms = 1 , 2
                                 wtmp(ms) = wtmp(ms) + dimag(tauc(lamp,lamp1))*dimag(taub(lamp2,lamp3))*waes(ms)
                              enddo
!
!     sumation over lam(valence)
                           enddo
                        enddo
                     enddo
                  enddo
                  do ms = 1 , 2
                     iaes(icst,ms,ied) = iaes(icst,ms,ied) + wgte*wtmp(ms)
                  enddo
!     iebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebiebie
!             end if
               enddo
!     iecieciecieciecieciecieciecieciecieciecieciecieciecieciecieciecie
               epsd = 16.D0*(dreal(etabfin(ied,1))+c_au**2)/(2.D0*dreal(etabfin(ied,1))+c_au**2)
               do ms = 1 , 2
                  iaes(icst,ms,ied) = iaes(icst,ms,ied)*epsd
               enddo
!
               if ( iprint>0 ) then
                  write (6,*) iaes(icst,1,ied)
                  write (6,*) iaes(icst,2,ied)
               endif
!write spectrum
               write (7,'(1x,e12.5,2x,e12.5,a)') dreal(etabfin(ied,1)) , dreal(etabfin(ied,1)) - dreal(etabfin(nefin,1))
               write (7,'(20x,6(1x,e18.10))') (iaes(icst,is,ied),is=1,2) , iaes(icst,1,ied) + iaes(icst,2,ied)
!
!bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb buildbot
               if ( wrbuildbot .and. ied<=3 ) write (ifilbuildbot,99017) trim(routine) , it , icst , ied , etabfin(ied,1) ,        &
                  & (iaes(icst,is,ied),is=1,2)
!----------------------------------------------------------------formats
99017          format ('# buildbot: ',a,': aes intensity for it =',i5,'  icst =',i3,/,'#',10x,'energy  ied =',i5,' etabfin = ',    &
                     & 2F10.6,/,(1pe22.14))
!bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb buildbot
            enddo
            exit spag_dispatchloop_1
         end select
      enddo spag_dispatchloop_1
!
!corecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecorecore
   enddo
!
   deallocate (iaes,mmed,mmee,rint,taub,tauc,mmetildd,mmetilde)
   deallocate (tautleed,difme,lammagcpl,me,amecif,amecig)
   deallocate (etabfin,eme)
99011 format (i9,i5,i5,i4,2x,2(g21.14))
99012 format (a10,i10)
end subroutine spAES
