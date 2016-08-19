!浅水方程式系　外部重力波による津波の再現
program tsunami

implicit none

integer,parameter::r8=8
integer::i,j,k,l,m,n

!file
integer::unml=66,ulog=6,uwrite=55,utikei=44
character::file_name*99

! parameter
integer::Nx,Ny
real(kind=r8)::dx,dy
real(kind=r8)::Lx,Ly

! time step 
integer::nend,iscan
real(kind=r8)::dt


! physical constant value
real(kind=r8),parameter::grav=9.8d0
real(kind=r8),parameter::H0=4000.0d0
real(kind=r8):: pi=4.d0*datan(1.d0)

! allocate value
real(kind=r8),allocatable,dimension(:,:)::sshb,ssh,ssha
real(kind=r8),allocatable,dimension(:,:)::h,he
integer      ,allocatable,dimension(:,:)::mask


! tikei 
integer::nisl
integer::i_tikei
character(99)::tikei_file
real(kind=r8),allocatable,dimension(:)::lon(:),lat(:)

real(kind=r8)::dtdx,dtdy

integer,parameter::slv=1
real(kind=r8)::tmp1(-slv:slv,-slv:slv)
real(kind=r8)::tmp2(-slv:slv,-slv:slv)


namelist/param/Nx,Ny,Lx,Ly
namelist/timestep/nend,dt,iscan
namelist/tikei/i_tikei,tikei_file

open(unml,file='./namelist.dat',status='old')
read(unml,param)
write(ulog,param)
read(unml,timestep)
write(ulog,timestep)
read(unml,tikei)
write(ulog,tikei)


allocate(sshb(0:Nx+1,0:Ny+1),ssh(0:Nx+1,0:Ny+1),ssha(0:Nx+1,0:Ny+1),h(0:Nx+1,0:Ny+1),he(0:Nx+1,0:Ny+1),mask(0:Nx+1,0:Ny+1))
allocate(lon(Nx),lat(Ny))

dx = Lx/Nx ; dy = Ly/Ny

!check
write(ulog,*)'dx,dy=',dx,dy

dtdx=dt/dx ; dtdy=dt/dy

call tikei_condition(h,mask,lon,lat)

! call isl_condition(mask,nisl)
! write(ulog,*)'nisl =',nisl

call CFL_condition_check(dt,dx,dy,h)

n=0
call initial_condition(sshb,ssh,ssha)

call write_file(n,ssh)

! time develop +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++

loop_time : do 
   
   n=n+1
   
   do i=1,Nx
      do j=1,Ny

         if(mask(i,j) .ge. 1)then 

            call temporary('DEPTH',tmp1,i,j)
            call temporary('SSH'  ,tmp2,i,j)

            ! ssha(i,j)=2.d0*ssh(i,j)-sshb(i,j) &
                 
            !      + (   ( he(i+1,j)-he(i-1,j) ) * ( ssh(i+1,j)              -ssh(i-1,j) )/4.d0 &  ! x 
            !          + ( he(i,j)             ) * ( ssh(i+1,j)-2.d0*ssh(i,j)+ssh(i-1,j) )      &
            !          ) * grav*(dtdx**2) &
                 
            !      + (   ( he(i,j+1)-he(i,j-1) ) * ( ssh(i,j+1)              -ssh(i,j-1) )/4.d0 &  ! y
            !          + ( he(i,j)             ) * ( ssh(i,j+1)-2.d0*ssh(i,j)+ssh(i,j-1) )      &
            !          ) * grav*(dtdy**2) 

            !! deep bottom 
            tmp1=tmp2+tmp1

            ssha(i,j)=2.d0*ssh(i,j)-sshb(i,j) &
                 
                 + (   ( tmp1(1,0)-tmp1(-1,0) ) * ( tmp2(1,0)               -tmp2(-1,0) )/4.d0 &  ! x 
                     + ( tmp1(0,0)            ) * ( tmp2(1,0)-2.d0*tmp2(0,0)+tmp2(-1,0) )      &
                     ) * grav*(dtdx**2) &
                 
                 + (   ( tmp1(0,1)-tmp1(0,-1) ) * ( tmp2(0,1)               -tmp2(0,-1) )/4.d0 &  ! y
                     + ( tmp1(0,0)            ) * ( tmp2(0,1)-2.d0*tmp2(0,0)+tmp2(0,-1) )      &
                     ) * grav*(dtdy**2) 

         endif

      enddo
   enddo

   !! time step
   sshb(1:Nx,1:Ny) = ssh(1:Nx,1:Ny)
   ssh(1:Nx,1:Ny)  = ssha(1:Nx,1:Ny)

   ! output *** *** *** *** *** *** *** *** ***
   if ( mod(n,iscan) .eq.0 )then

      call write_file(n,ssha)

      write(ulog,*)'check ;',n,maxval(ssha),minval(ssha)

   endif
   ! output *** *** *** *** *** *** *** *** ***

   if ( n .ge. nend)exit loop_time


enddo loop_time

! +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++

contains
  
  subroutine tikei_condition(h,mask,lon,lat)
    real(kind=r8),intent(inout)::h(0:Nx+1,0:Ny+1)
    integer,intent(inout):: mask(0:Nx+1,0:Ny+1)
    real(kind=r8),intent(inout)::lon(1:Nx),lat(1:Ny)

    integer::i,j
    integer::i0
    real(kind=r8)::ratio
    
    
    !! mask = 1 ; ocean
    !! mask = 0 ; land

    !! tikei example
    !!
    !! 00000000000000000000000000
    !! 01111111111111111111111110
    !! 01111111111111111111111110                    
    !! 01111111110000111111111110                    
    !! 01111111110000111111111110                    
    !! 01111111110000111111111110                    
    !! 01111111111111111111111110                    
    !! 01111111111111111111111110
    !! 00000000000000000000000000

    mask=1
    
    mask(0,:)=0 ; mask(Nx+1,:)=0 ; mask(:,0)=0 ; mask(:,Ny+1)=0   
    
    !! block
    ! mask(Nx/2-10:Nx/2+10,Ny/2-10:Ny/2+10)=0
    
       !! thread
    mask(Nx/2,:)=0
    mask(Nx/2,Ny/2-5:Ny/2+5)=1
    
       
       i0=Nx/2
       
       !! step ratio ( 0 < ratio < 1)
       ratio=0.75d0
       
       do i=0,Nx+1
          
          do j=0,Ny+1
          
          if(mask(i,j) .ge. 1 )then 
             !          h(i,j) = H0/2.d0 !H0+dble(i0-i)/(Nx-i0+1)*H0*0.99d0
             h(i,j) = H0*(1.d0 - ratio*1.d0/2.d0*(1.0d0 + tanh( dble(i-i0) ) ) )

          else
!             h(i,j) = H0*(1.d0 - ratio*1.d0/2.d0*(1.0d0 + tanh( dble(i-i0) ) ) )
             h(i,j)=0.0d0
          
          endif

       enddo

    enddo


    if(i_tikei .eq. 1 )then 

       open(utikei,file=tikei_file,action='read',status='old',form='unformatted')
       read(utikei)lat,lon,h(1:Nx,1:Ny)
       close(utikei)
       
       write(ulog,*)'==== === === tikei file is reading . === === ===='
       write(ulog,*)'longitude=',maxval(lon),minval(lon ,mask=lon .gt. 0 )
       write(ulog,*)'latitude =',maxval(lat),minval(lat ,mask=lat .gt. 0 )
       write(ulog,*)'dephth   =',maxval( h ),minval( h  ,mask= h  .gt. 0 )

       mask=1
       do i=1,Nx
          do j=1,Ny
             if(h(i,j).le.10)then
                mask(i,j)=0
                h(i,j)=0.d0
             endif
          enddo
       enddo

       write(ulog,*)'==== === ===  === === ==== === === ==== === === ==== === === ===='

    endif
    
    call tikei_check(h)

  end subroutine tikei_condition

  subroutine initial_condition(sshb,ssh,ssha)
    integer::i,j
    real(kind=r8)::sigma=10.d0
    real(kind=r8)::ssh0=1.d0
    real(kind=r8),intent(inout)::sshb(0:Nx+1,0:Ny+1),ssh(0:Nx+1,0:Ny+1),ssha(0:Nx+1,0:Ny+1)
    integer::i0 , j0

    real(kind=r8)::r

! saka
    i0=Nx/2
    j0=4*Ny/5
    
    ssh=0.d0
    sshb=0.d0
    ssha=0.d0
    
    do i=0,Nx+1
       do j=0,Ny+1
          
          !! gauss
          ssh (i,j) = ssh0*exp(-( (i-i0)**2+(j-j0)**2 ) /sigma**2 )
          sshb(i,j) = ssh0*exp(-( (i-i0)**2+(j-j0)**2 ) /sigma**2 )
          ssha(i,j) = ssh0*exp(-( (i-i0)**2+(j-j0)**2 ) /sigma**2 )
          
          !! plane // y axis
          ! ssh (i,j)=ssh0*exp(-( (i-i0)**2 )/sigma**2 )
          ! sshb(i,j)=ssh0*exp(-( (i-i0)**2 )/sigma**2 )
          ! ssha(i,j)=ssh0*exp(-( (i-i0)**2 )/sigma**2 )

          !! plane // y + x -0.3 = 0

          !! the distance from the line.
          ! r=abs( 1.d0 * i/Nx + 1.d0 * j/Ny -0.3d0 ) / sqrt((1.d0/Nx)**2 + (1.d0/Ny)**2)
          ! ssh (i,j)=ssh0*exp(-( r**2 )/sigma**2 )
          ! sshb(i,j)=ssh0*exp(-( r**2 )/sigma**2 )
          ! ssha(i,j)=ssh0*exp(-( r**2 )/sigma**2 )
          

          if(abs(ssh (i,j)) .lt. 1.0d-10 ) ssh (i,j)=0.0d0
          if(abs(sshb(i,j)) .lt. 1.0d-10 ) sshb(i,j)=0.0d0
          if(abs(ssha(i,j)) .lt. 1.0d-10 ) ssha(i,j)=0.0d0

       enddo
    enddo

  end subroutine initial_condition

  !! macro +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ 

  subroutine temporary(varname,data,ix,iy)
    character(*),intent(in)  :: varname
    real(kind=r8),intent(out):: data(-1:1,-1:1)
    integer,intent(in)    :: ix,iy
    
    integer :: id
    
    
    id = lookup_id(varname)
    
    select case(id)
    case(0)
       data(:,:)=h(ix-slv:ix+slv,iy-slv:iy+slv)
    case(1)
       data(:,:)=ssh(ix-slv:ix+slv,iy-slv:iy+slv)
    end select
    
    call boundary_condition(data,ix,iy)
    
  end subroutine temporary
  
  function lookup_id(varname)
    character(*),intent(in) :: varname
    integer :: lookup_id
    
    select case(varname)
    case('DEPTH')
       lookup_id=0
    case('SSH')
       lookup_id=1
    end select
    
  end function lookup_id
    
  
  subroutine boundary_condition(data,ix,iy)
    integer::i,j
    integer,intent(in)::ix,iy
    real(kind=r8),intent(inout)::data(-1:1,-1:1)

    !! free end
    if(mask(ix-1,iy)*mask(ix+1,iy).eq.0)then
      
       if(mask(ix-1,iy).eq.0)then 
          data(-1,0)=data(1,0)
       elseif(mask(ix+1,iy).eq.0)then 
          data(1,0)=data(-1,0)
       else
          write(ulog,*)'ERROR : narrow region. (i,j)=',ix,iy
          stop
       endif

    endif

    if(mask(ix,iy-1)*mask(ix,iy+1).eq.0)then
      
       if(mask(ix,iy-1).eq.0)then 
          data(0,-1)=data(0,1)
       elseif(mask(ix,iy+1).eq.0)then 
          data(0,1)=data(0,-1)
       else
          write(ulog,*)'ERROR : narrow region. (i,j)=',ix,iy
          stop
       endif

    endif
    
    ! west
    !    pa(0,:)=pa(1,:) ; pa(Nx+1,:)=pa(Nx,:) ;  pa(:,Ny+1)=pa(:,Ny) ;  pa(:,0)=pa(:,1)
    
    !    west
    
    !     do j = 1, Ny   
    !        pa(0,j)=(4.d0*pa(1,j)-pa(2,j))/3.d0
    !     enddo
    
    ! !    east
    
    !     do j = 1, Ny   
    !        pa(Nx+1,j)=(4.d0*pa(Nx,j)-pa(Nx-1,j))/3.d0
    !     enddo
    
    ! !    north
    
    !     do i = 1, Nx   
    !        pa(i,Ny+1)=(4.d0*pa(i,Ny)-pa(i,Ny-1))/3.d0
    !     enddo
    
    ! !    surth
    
    !     do i = 1, Nx   
    !        pa(i,0)=(4.d0*pa(i,1)-pa(i,2))/3.d0
    !     enddo
    
  end subroutine boundary_condition


  subroutine tikei_check(h)
    implicit none
    real(kind=r8),intent(in)::h(0:Nx+1,0:Ny+1)
    integer::i,j
    
    
    !! tikei check === === === === === === === === === === === === === === === === === ===
    write(ulog,*)'  __--^^--__--~~ tikei check __--^^--__--~~ '
    
    do j=0,Ny+1,10
       write(ulog,'(50f6.0)')(h(i,j) ,i=0,Nx+1,10)
    enddo
    
    write(ulog,*)'  __--^^--__--~~ tikei check __--^^--__--~~ '
    !! === === === === === === === === === === === === === === === === === === === === ===
    
  end subroutine tikei_check


  subroutine write_file(n,p)
    integer,intent(in)::n
    real(kind=r8),intent(in)::p(0:Nx+1,0:Ny+1)
    real(kind=r8)           ::pdum(0:Nx+1,0:Ny+1)
    character::file_name*99
    character::file_number*5
    integer::i,j
    
    do i=1,Nx
       do j=1,Ny
          if(mask(i,j).eq.0)then
             pdum(i,j)=-999.d0
          else
             pdum(i,j)=p(i,j)
          endif
       enddo
    enddo

    write(file_number,'(i5.5)')n
    
!saka
!    write(file_name,'(a,a,a)')'./data/ts_',adjustl(file_number),'.dat'
    write(file_name,'(a,a,a)')'./data/ts_saka_',adjustl(file_number),'.dat'
    
    open(uwrite,file=file_name,status='replace',form='unformatted')

    write(uwrite)n*dt,pdum(1:Nx,1:Ny)

    close(uwrite)
    
  end subroutine write_file

  subroutine CFL_condition_check(dt,dx,dy,h)
    integer::i,j
    real(kind=r8),intent(inout)::dt,dx,dy,h(0:Nx+1,0:Ny+1)
    logical::flg=.false.
    
    do i=1,Nx
       do j=1,Ny
          if ( dt/dx .gt. 1.d0/(sqrt(grav*h(i,j))) & 
               .or. dt/dy .gt. 1.d0/(sqrt(grav*h(i,j))) )then 
             
             write(ulog,*)' Warning ! : CFL consition is not satisfied at',i,j
             flg=.true.
             
          endif
       enddo
    enddo
    
    if(flg) stop
    
  end subroutine CFL_condition_check


  subroutine isl_condition(mask,nisl)
    integer::mask(0:Nx+1,0:Ny+1),nisl
    integer::i,j        
    ! nisl=product(mask,mask=mask.gt.0)
    ! nisl=log(dble(nisl))/log(2.d0)

    nisl=0
    do i=0,Nx+1
       do j=0,Ny+1
          if(mask(i,j).eq.2)nisl=nisl+1
       enddo
    enddo
    
  end subroutine isl_condition


  
end program tsunami

