subroutine s_ess(nsd,nen,ne,ndof,eldof,ien,ng,tmpgdof,fext)
implicit none
integer :: file,ien(eldof,ne),ne,ng,stat,dof,nen,tmpgdof(ndof),nsd,ndof,eldof,ng_true(ndof)
integer,allocatable :: gdof(:)
integer :: a(nen),b,el,bc,p,q,i,j,c
real(8) :: fext(ndof),bf(nsd)
ng=0
ng_true(:)=0
file=11
open(file, FILE="sbc_solid.in", STATUS="old",action="read")
tmpgdof(:)=0
fext(:)=0.0d0
do while (i .lt. huge(1))
read(11,*,end=150) el, a(1:nen), bc
	do b=1,nen
		if (b .gt. 2) then
			c=b-2
		else
			c=b+2
		endif
		if (a(b) .eq. 1) then
			if (bc .gt. 10000) then
				if (bc .eq. 10110) then
					ng=ng+1
					dof=nsd*(c-1)+1
					p=ien(dof,el)
					tmpgdof(ng)=p
					ng=ng+1
					dof=nsd*(c-1)+2
					q=ien(dof,el)
					tmpgdof(ng)=q
				else if (bc .eq. 10100) then
					ng=ng+1
					dof=nsd*(c-1)+1
					p=ien(dof,el)
					tmpgdof(ng)=p
				else if (bc .eq. 10010) then
					ng=ng+1
					dof=nsd*(c-1)+2
					q=ien(dof,el)
					tmpgdof(ng)=q
				endif
			else
					dof=nsd*(c-1)+1
					p=ien(dof,el)
					fext(p)=fext(p)+0.0d0
			endif
		endif	
	enddo
enddo
!remove nonunique entries
150 do i=1,ndof
		do j=1,ndof
			if ((tmpgdof(i) .eq. tmpgdof(j)) .and. (i .ne. j))then
				tmpgdof(j)=0
			endif
		enddo
	enddo
where (tmpgdof .ne. 0)
	ng_true=1
end where

ng=sum(ng_true)
write(*,*) 'ng= ', ng
end subroutine s_ess