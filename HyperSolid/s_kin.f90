subroutine s_kin(nee,ndof,acc,Mass,km,fkin,ka)
implicit none
integer :: ndof,nee,i,j
real(8),intent(in) :: Mass(ndof,ndof)
real(8) :: fkin(ndof),acc(ndof),km(ndof,ndof),ka(nee,nee)
ka(1:ndof,1:ndof)=ka(1:ndof,1:ndof)+km(:,:)
!write(*,*) shape(Mass),shape(acc)
fkin(:)=0.0d0
!write(*,*) fkin
do i=1,ndof
!write(*,*) 'i= ',i
	do j=1,ndof
		!write(*,*) 'j= ', j
		!write(*,*) 'Mass= ', Mass(i,j)
		!write(*,*) 'acc= ',acc(j)
		!write(*,*) 'prod= ',Mass(i,j)*acc(j)
		fkin(i)=fkin(i)+Mass(i,j)*acc(j)
	enddo
enddo
write(4,*) '------acceleration-----'
do i=1,ndof
write(4,*) acc(i)
enddo
write(4,*) '-'
write(21,*) '-------Mass Matrix---------'
do i=1,ndof
	write(21,*) Mass(i,:)
enddo
end subroutine s_kin