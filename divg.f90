subroutine divg 
use param
use local_arrays, only: vx,vy,vz,dph
use mpi_param, only: kstart,kend
implicit none
integer :: jc,jp,kc,kp,ic,ip
real    :: usdtal,dqcap   

usdtal = 1.d0/(dt*al)

do kc=kstart,kend
  kp=kc+1
  do jc=1,n2m
    jp=jpv(jc)
      do ic=1,n1m
        ip=ipv(ic)
        dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
              +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
              +(vz(ic,jc,kp)-vz(ic,jc,kc))*dx3
        dph(ic,jc,kc)=dqcap*usdtal
      enddo
   enddo
enddo
end

subroutine divgck(qmax)
use param
use local_arrays, only: vy,vz,vx
use mpih
use mpi_param, only: kstart,kend
implicit none
real,intent(out) :: qmax
integer :: jc,kc,kp,jp,ic,ip
real    :: dqcap,dvol,my_qmax
  

dvol=1.d0/(dx1*dx2*dx3)
qmax=0.d0                                                     

do kc=kstart,kend
  kp=kc+1
  do jc=1,n2m
    jp=jpv(jc)
      do ic=1,n1m
        ip=ipv(ic)
        dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
              +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
              +(vz(ic,jc,kp)-vz(ic,jc,kc))*dx3
        qmax = max(abs(dqcap),qmax)          
enddo
enddo
enddo

call MpiAllMaxRealScalar(qmax)

end

subroutine vmaxv 
use param
use local_arrays, only: vy,vz,vx
use mpi_param, only: kstart,kend
use mpih
implicit none
integer :: jc,kc,ic

 vmax=-100.d0

 do kc=kstart,kend
  do jc=1,n2m
    do ic=1,n1m
     vmax(1) = max(vmax(1),abs(vx(ic,jc,kc)))
     vmax(2) = max(vmax(2),abs(vy(ic,jc,kc)))
     vmax(3) = max(vmax(3),abs(vz(ic,jc,kc)))
   enddo
  enddo
 enddo

call MpiAllMaxRealScalar(vmax(1))
call MpiAllMaxRealScalar(vmax(2))
call MpiAllMaxRealScalar(vmax(3))
end


subroutine divgloc 
use param
use local_arrays, only: vy,vz,vx
use mpih
use mpi_param, only: kstart,kend
implicit none
integer :: jc,kc,kp,jp,ic,ip
real    :: dqcap

resid = 1.0d0-2
  
do kc=kstart,kend
  kp=kc+1
  do jc=1,n2m
    jp=jpv(jc)
      do ic=1,n1m
        ip=ipv(ic)
        dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
              +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
              +(vz(ic,jc,kp)-vz(ic,jc,kc))*dx3
        if (abs(dqcap).gt.resid) then
           write(*,*) ic,jc,kc,myid
        endif
enddo
enddo
enddo

return     
end
