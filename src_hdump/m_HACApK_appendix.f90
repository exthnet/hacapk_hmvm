module m_HACApK_hdump
  use m_HACApK_use
  use m_HACApK_base
  use m_HACApK_calc_entry_ij
!  use omp_lib
  implicit integer*4(i-n)
!  integer omp_get_thread_num, omp_get_num_threads
!  external omp_get_thread_num, omp_get_num_threads

contains


  !***HACApK_info_leafmtx
  subroutine HACApK_info_leafmtx(st_leafmtxp,st_ctl,mpinr)
    implicit none
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    type(st_HACApK_lcontrol) :: st_ctl
    integer :: mpinr
    integer*4 :: nd,nlf,ktmax,ip,ltmtx,kt,ndl,ndt,nstrtl,nstrtt,ierr
    character*32 :: fname

    integer :: min_kt, min_ndl, min_ndt
    integer :: max_kt, max_ndl, max_ndt
    integer :: min_long, max_short

    write(*,*)"HACApK_info_leafmtx: begin"

    nd=st_leafmtxp%nd; nlf=st_leafmtxp%nlf; ktmax=st_leafmtxp%ktmax

    ! approximate submatrix
    min_kt = 999999;  min_ndl = 999999;  min_ndt = 999999
    max_kt = 0;  max_ndl = 0;  max_ndt = 0
    do ip=1,nlf
       ltmtx=st_leafmtxp%st_lf(ip)%ltmtx  ; kt    =st_leafmtxp%st_lf(ip)%kt
       ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt
       nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
       if(ltmtx==1)then ! Low-rank matrix
          if(min_ndl > ndl)min_ndl = ndl;      if(min_ndt > ndt)min_ndt = ndt
          if(max_ndl < ndl)max_ndl = ndl;      if(max_ndt < ndt)max_ndt = ndt
          if(min_kt > kt)min_kt = kt;      if(max_kt < kt)max_kt = kt
       elseif(ltmtx==2)then ! Dense matrix
       endif
    enddo
    print*,'approx min ndl=',min_ndl,'ndt=',min_ndt,'kt=',min_kt
    print*,'approx max ndl=',max_ndl,'ndt=',max_ndt,'kt=',max_kt

    ! small dense submatrix
    min_kt = 999999;  min_ndl = 999999;  min_ndt = 999999
    max_kt = 0;  max_ndl = 0;  max_ndt = 0
    min_long = 999999; max_short = 0
    do ip=1,nlf
       ltmtx=st_leafmtxp%st_lf(ip)%ltmtx  ; kt    =st_leafmtxp%st_lf(ip)%kt
       ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt
       nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
       if(ltmtx==1)then ! Low-rank matrix
       elseif(ltmtx==2)then ! Dense matrix
          if(min_ndl > ndl)min_ndl = ndl;      if(min_ndt > ndt)min_ndt = ndt
          if(max_ndl < ndl)max_ndl = ndl;      if(max_ndt < ndt)max_ndt = ndt
          if(ndl<ndt)then
             if(max_short < ndl)max_short = ndl
             if(min_long > ndt)min_long = ndt
          else
             if(max_short < ndt)max_short = ndt
             if(min_long > ndl)min_long = ndl
          endif
       endif
    enddo
    print*,'dense  min ndl=',min_ndl,'ndt=',min_ndt
    print*,'dense  max ndl=',max_ndl,'ndt=',max_ndt
    print*,'dense  max_short=',max_short,'min_long=',min_long
    write(*,*)"HACApK_info_leafmtx: end"
  end subroutine HACApK_info_leafmtx

  !***HACApK_dump_leafmtx
  subroutine HACApK_appendix_dump_leafmtx(st_leafmtxp,st_ctl,mpinr)
    implicit none
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    type(st_HACApK_lcontrol) :: st_ctl
    integer :: mpinr
    integer*4 :: nd,nlf,ktmax,ip,ltmtx,kt,ndl,ndt,nstrtl,nstrtt,ierr
    character*32 :: fname

    write(*,*)"HACApK_appendix_dump_leafmtx: begin"

    write(fname,'(a,i0,a)')'hmatrix_',mpinr,'.bin'
    write(*,*)"file dump to ",fname
    nd=st_leafmtxp%nd; nlf=st_leafmtxp%nlf; ktmax=st_leafmtxp%ktmax
    open( 11, file=fname, action='write', iostat=ierr, form="unformatted" )

    ! for debug, dump only first 2 leavs
    !    write(11) nd,1,ktmax
    !    do ip=1,2

    write(11) nd,nlf,ktmax
    do ip=1,nlf
       ltmtx=st_leafmtxp%st_lf(ip)%ltmtx  ; kt    =st_leafmtxp%st_lf(ip)%kt
       ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt
       nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
       write(11) ltmtx,ndl,ndt,nstrtl,nstrtt,kt
       if(ltmtx==1)then ! Low-rank matrix
          write(11)st_leafmtxp%st_lf(ip)%a1
          write(11)st_leafmtxp%st_lf(ip)%a2
       elseif(ltmtx==2)then ! Dense matrix
          write(11)st_leafmtxp%st_lf(ip)%a1
       endif
    enddo
    close( 11 )
    write(*,*)"HACApK_appendix_dump_leafmtx: end"
  end subroutine HACApK_appendix_dump_leafmtx

  ! modified solve function
  integer function HACApK_appendix_solve(st_leafmtxp,st_bemv,st_ctl,rhs,sol,ztol)
    include 'mpif.h'
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    type(st_HACApK_lcontrol) :: st_ctl
    type(st_HACApK_calc_entry) :: st_bemv
    real*8 :: rhs(st_bemv%nd),sol(st_bemv%nd),ztol
    real*8,pointer :: param(:)
    real*8,dimension(:),allocatable :: u,b,www,ao
    integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)

    lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:); param=>st_ctl%param(:)
    mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1); nthr=lpmd(20)
    param(91)=ztol
    if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_solve start'
    nofc=st_bemv%nd;nffc=1;ndim=3
    nd=nofc*nffc
    if(st_ctl%param(1)>1) write(*,*) 'irank=',mpinr,' lthr=',lthr(0:nthr-1)
    allocate(u(nd),b(nd)); u(:nd)=sol(lod(:nd)); b(:nd)=rhs(lod(:nd))
    if(param(61)==3)then
       !   do il=ndnr_s,ndnr_e
       do il=1,nd
          u(il)=u(il)/st_bemv%ao(lod(il))
          b(il)=b(il)*st_bemv%ao(lod(il))
       enddo
    endif
    if(param(83)>0) then
       allocate(ao(nd))
       do il=1,nd
          zzz=HACApK_entry_ij(il,il,st_bemv)
          ao(il)=1.0d0/zzz
       enddo
       call MPI_Barrier( icomm, ierr )
       ! dump hmatrix
       ! dump info
       print*,'HACApK_info_leafmtx'
       call HACApK_info_leafmtx(st_leafmtxp,st_ctl,mpinr)
       ! dump matrix
       print*,'HACApK_dump_leafmtx'
       st_leafmtxp%nd = nd
       call HACApK_appendix_dump_leafmtx(st_leafmtxp,st_ctl,mpinr)
    end if
  end function HACApK_appendix_solve

  ! modified gensolv function
  integer function HACApK_appendix_gensolv(st_leafmtxp,st_bemv,st_ctl,gmid,rhs,sol,ztol)
    include 'mpif.h'
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    type(st_HACApK_lcontrol) :: st_ctl
    type(st_HACApK_calc_entry) :: st_bemv
    real*8 :: gmid(st_bemv%nd,3),rhs(st_bemv%nd),sol(st_bemv%nd),ztol
    real*8 :: ww(st_bemv%nd),aww(st_bemv%nd),aww1(st_bemv%nd)
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)

    mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
    icomm=st_ctl%lpmd(1)
    lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,gmid,ztol)
    call MPI_Barrier( icomm, ierr )
    ! dump H matrix
    lrtrn=HACApK_appendix_solve(st_leafmtxp,st_bemv,st_ctl,rhs,sol,ztol)
    call MPI_Barrier( icomm, ierr )
    if(st_ctl%param(1)>0 .and. mpinr==0)  write(6,2000) 'lfmtx; time_AX_once  =',(en_measure_time_ax - st_measure_time_ax)/st_ctl%param(99)
9999 continue
    HACApK_appendix_gensolv=lrtrn
  endfunction HACApK_appendix_gensolv

end module m_HACApK_hdump
