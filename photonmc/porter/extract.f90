implicit real*8(a-h,o-z) 
character(len=20) :: word 
double precision :: x, y, z, temp,grad, srad 
integer :: np 

call system('wc -l 17by17by34_a0.1_phi16_theta8.dat') 
write(*,*) 'enter number of rows' 
read(*,*) np 
open(11,file='17by17by34_a0.1_phi16_theta8.dat',status='old') 
open(12,file='centerline_srad.dat',status='unknown') 
read(11,*)
 
do i=1, np 
  read(11,*) x, y, z, temp,grad, srad  
  if(x.eq.1.0d0.and.y.eq.1.0d0) write(12,*) x, y, z, temp,grad, srad 
end do 

end  
 
