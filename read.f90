program read_data

implicit none

integer,parameter::uread=44,ulog=6,uwrite=66,r8=8,unml=33

integer::i,j,k
integer::dlen
integer,parameter::nmax=1000000

real(kind=r8)::dlon(nmax),dlat(nmax),dh(nmax)

integer::Nx,Ny
real(kind=r8)::Lx,Ly
real(kind=r8),allocatable,dimension(:)::lon(:),lat(:)
real(kind=r8),allocatable,dimension(:)::h(:,:)

namelist/param/Nx,Ny,Lx,Ly

open(unml,file='./namelist.dat',status='old')
read(unml,param)
write(ulog,param)

!! initialize
dlon=0.d0 ; dlat=0.d0 ; dlen=0 ; dh=0.d0
allocate(lon(Nx),lat(Ny),h(Nx,Ny))


open(uread,file='./jodc-depth500mesh-20160813165833.txt',action='read',status='old')

!! reading
do i=1,nmax
   read(uread,*,end=100)k,dlat(i),dlon(i),dh(i)
   dlen=dlen+1
enddo
100 write(ulog,*)'length=',dlen

close(uread)

write(ulog,*)'longitude=',maxval(dlon),minval(dlon ,mask=dlon .gt. 0 )
write(ulog,*)'latitude =',maxval(dlat),minval(dlat ,mask=dlat .gt. 0 )

do j=1,Ny
   lat(j)=( minval(dlat,mask=dlat .gt. 0)*(j-1) + maxval(dlat)*(Ny-j+1) )/Ny
enddo

! !check 
! write(ulog,*)lat

do i=1,Nx
   lon(i)=( minval(dlon,mask=dlon .gt. 0)*(Nx-i+1) + maxval(dlon)*(i-1) )/Nx
enddo

!check 
!write(ulog,*)lon

h=0.d0

loop_lon : do i=1,Nx
   loop_lat : do j=1,Ny

      if (mod(i,50) .eq. 0 .and. mod(j,50) .eq. 0 ) write(ulog,*)'check : ',i,j

      loop_k : do k=1,dlen
         
         if( lat(j) .ge. dlat(k) .and.  lat(j+1) .lt. dlat(k) &
              .and. &
             lon(i) .le. dlon(k) .and.  lon(i+1) .gt. dlon(k) )then
            
            h(i,j)=dh(k)

            !write(ulog,*)k
            
            exit loop_k

         endif

      enddo loop_k
  
 enddo loop_lat

! if (mod(i,50) .eq. 0)write(ulog,*)'check : ',h(i,:)

enddo loop_lon



!!! lat ; lon ; h
open(uwrite,file='tikei_hokan.dat',form='unformatted',status='replace')

write(uwrite)lat(1:Ny),lon(1:Nx),h(1:Nx,1:Ny)

! do j=1,Ny
!    write(uwrite,*)(h(i,j),i=1,Nx)
! enddo

close(uwrite)
   
end program read_data

