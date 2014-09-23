subroutine s_gdof(ndof,ng,tmpgdof,gdof,gx)
implicit none
integer :: i,b,ndof,tmpgdof(ndof),ng,gdof(ng)
real(8) :: gx(ng)
i=1
do b=1,ndof
	if (tmpgdof(b) .ne.  0) then
		gdof(i)=tmpgdof(b)
		i=i+1
	endif
enddo
gx(:)=0.0d0
end subroutine s_gdof