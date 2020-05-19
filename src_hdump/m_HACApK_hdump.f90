module m_HACApK_hdump
  use m_HACApK_use
  use m_HACApK_base
  use m_HACApK_calc_entry_ij
  !  use omp_lib
  implicit integer*4(i-n)

contains

  ! ######## ######## ######## ######## ######## ######## ######## ########
  ! HACApK_info_leafmtx
  ! show hmatrix information
  ! ######## ######## ######## ######## ######## ######## ######## ########
  subroutine HACApK_info_leafmtx(st_leafmtxp,st_ctl,mpinr)
    implicit none
    type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    type(st_HACApK_lcontrol),intent(in) :: st_ctl
    integer :: mpinr
    integer*4 :: nd, nlf, ktmax, ip, ltmtx, kt, ndl, ndt, nstrtl, nstrtt, ierr
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

  ! ######## ######## ######## ######## ######## ######## ######## ########
  ! HACApK_hdump_dump_leafmtx
  ! write hmatrix to file (hmatrix_1.bin)
  ! ######## ######## ######## ######## ######## ######## ######## ########
  subroutine HACApK_hdump_dump_leafmtx(st_leafmtxp,st_ctl,mpinr)
    implicit none
    type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    type(st_HACApK_lcontrol),intent(in) :: st_ctl
    integer :: mpinr
    integer*4 :: nd,nlf,ktmax,ip,ltmtx,kt,ndl,ndt,nstrtl,nstrtt,ierr
    character*32 :: fname

    write(*,*)"HACApK_hdump_dump_leafmtx: begin"

    write(fname,'(a,i0,a)')'hmatrix_1.bin'
    write(*,*)"file dump to ",fname
    nd=st_leafmtxp%nd; nlf=st_leafmtxp%nlf; ktmax=st_leafmtxp%ktmax
    open( 11, file=fname, action='write', iostat=ierr, form="unformatted" )

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
    write(*,*)"HACApK_hdump_dump_leafmtx: end"
  end subroutine HACApK_hdump_dump_leafmtx

  ! ######## ######## ######## ######## ######## ######## ######## ########
  ! HACApK_hdump_leafmtx2
  ! write hmatrix to file (hmatrix_2.bin)
  ! ######## ######## ######## ######## ######## ######## ######## ########
  subroutine HACApK_hdump_dump_leafmtx2(st_leafmtxp,st_ctl,mpinr)
    implicit none
    type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    type(st_HACApK_lcontrol),intent(in) :: st_ctl
    integer :: mpinr
    integer*4 :: nd,nlf,ktmax,ip,ltmtx,kt,ndl,ndt,nstrtl,nstrtt,ierr, len
    character*32 :: fname
    integer*4,dimension(:),allocatable :: buf

    write(*,*)"HACApK_hdump_dump_leafmtx2: begin"

    write(fname,'(a,i0,a)')'hmatrix_2.bin'
    write(*,*)"file dump to ",fname
    nd=st_leafmtxp%nd; nlf=st_leafmtxp%nlf; ktmax=st_leafmtxp%ktmax
    open( 11, file=fname, action='write', iostat=ierr, form="unformatted" )

    allocate(buf(nlf))

    len = 0
    do ip=1,nlf
       ltmtx=st_leafmtxp%st_lf(ip)%ltmtx
       if(ltmtx==1)then ! Low-rank matrix
          len = len + st_leafmtxp%st_lf(ip)%ndt * st_leafmtxp%st_lf(ip)%kt
          len = len + st_leafmtxp%st_lf(ip)%ndl * st_leafmtxp%st_lf(ip)%kt
       elseif(ltmtx==2)then ! Dense matrix
          len = len + st_leafmtxp%st_lf(ip)%ndl * st_leafmtxp%st_lf(ip)%ndt
       endif
    enddo

    write(11) nd,nlf,ktmax,len

    do ip=1,nlf
       buf(ip) =st_leafmtxp%st_lf(ip)%ltmtx
    enddo
    write(11) buf
    do ip=1,nlf
       buf(ip) =st_leafmtxp%st_lf(ip)%ndl
    enddo
    write(11) buf
    do ip=1,nlf
       buf(ip) =st_leafmtxp%st_lf(ip)%ndt
    enddo
    write(11) buf
    do ip=1,nlf
       buf(ip) =st_leafmtxp%st_lf(ip)%nstrtl
    enddo
    write(11) buf
    do ip=1,nlf
       buf(ip) =st_leafmtxp%st_lf(ip)%nstrtt
    enddo
    write(11) buf
    do ip=1,nlf
       buf(ip) =st_leafmtxp%st_lf(ip)%kt
    enddo
    write(11) buf
    do ip=1,nlf
       ltmtx=st_leafmtxp%st_lf(ip)%ltmtx
       if(ltmtx==1)then ! Low-rank matrix
          write(11)st_leafmtxp%st_lf(ip)%a1
          write(11)st_leafmtxp%st_lf(ip)%a2
       elseif(ltmtx==2)then ! Dense matrix
          write(11)st_leafmtxp%st_lf(ip)%a1
       endif
    enddo
    close( 11 )
    write(*,*)"HACApK_hdump_dump_leafmtx2: end"
  end subroutine HACApK_hdump_dump_leafmtx2

  ! ######## ######## ######## ######## ######## ######## ######## ########
  ! HACApK_hdump_solve
  ! modified solve function
  ! ######## ######## ######## ######## ######## ######## ######## ########
  subroutine HACApK_hdump_solve(st_leafmtxp, st_bemv, st_ctl, rhs, sol, ztol)
    implicit none
    include 'mpif.h'
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    type(st_HACApK_lcontrol), intent(in) :: st_ctl
    type(st_HACApK_calc_entry), intent(in) :: st_bemv
    real*8, intent(in) :: rhs(st_bemv%nd), sol(st_bemv%nd), ztol
    real*8,pointer :: param(:)
    real*8,dimension(:),allocatable :: u,b,www,ao
    integer*4,pointer :: lpmd(:), lnp(:), lsp(:), lthr(:), lod(:)
    integer :: icomm, ierr
    integer :: il, mpinr, nthr, nofc, nffc, nd
    real*8 :: zzz
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)

    lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp; lthr(0:) => st_ctl%lthr
    lod => st_ctl%lod(:); param=>st_ctl%param(:)
    mpinr=lpmd(3); icomm=lpmd(1); nthr=lpmd(20)
    if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_hdump_solve start'
    nd=st_bemv%nd
    if(st_ctl%param(1)>1) write(*,*) 'irank=',mpinr,' lthr=',lthr(0:nthr-1)
    call MPI_Barrier( icomm, ierr )
    ! dump info
    print*,'HACApK_info_leafmtx'
    call HACApK_info_leafmtx(st_leafmtxp,st_ctl,mpinr)
    ! dump matrix to files
    print*,'HACApK_dump_leafmtx'
    st_leafmtxp%nd = nd
    call HACApK_hdump_dump_leafmtx(st_leafmtxp,st_ctl,mpinr)
    call HACApK_hdump_dump_leafmtx2(st_leafmtxp,st_ctl,mpinr)
  end subroutine HACApK_hdump_solve

  ! ######## ######## ######## ######## ######## ######## ######## ########
  ! HACApK_hdump (modified gensolv function)
  ! ######## ######## ######## ######## ######## ######## ######## ########
  subroutine HACApK_hdump(st_leafmtxp, st_bemv, st_ctl, gmid, rhs, sol, ztol)
    implicit none
    include 'mpif.h'
    type(st_HACApK_leafmtxp), intent(in) :: st_leafmtxp
    type(st_HACApK_lcontrol), intent(in) :: st_ctl
    type(st_HACApK_calc_entry), intent(in) :: st_bemv
    real*8, intent(in) :: gmid(st_bemv%nd,3), rhs(st_bemv%nd), sol(st_bemv%nd), ztol
    integer :: icomm, ierr, ret
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)

    icomm = st_ctl%lpmd(1)
    ret = HACApK_generate(st_leafmtxp,st_bemv,st_ctl,gmid,ztol)
    call MPI_Barrier( icomm, ierr )
    ! dump H matrix
    call HACApK_hdump_solve(st_leafmtxp,st_bemv,st_ctl,rhs,sol,ztol)
    call MPI_Barrier( icomm, ierr )
  end subroutine HACApK_hdump

end module m_HACApK_hdump
