subroutine s_ess(ien_sbc)
use solid_variables
implicit none
integer :: file,tmpgdof(ndof_solid),ng_true(ndof_solid),ng
integer :: a(nen_solid),b,el,bc,p,q,i,j,c,dof
integer :: ien_sbc(ne_sbc_1,nen_solid+2)
ng=0
ng_true(:)=0
file=11
!open(file, FILE="sbc_solid.in", STATUS="old",action="read")
tmpgdof(:)=0
!write(*,*) 'ne_sbc_1= ',ne_sbc
do i=1,ne_sbc_1
!read(11,*,end=150) el, a(1:nen_solid), bc
el=ien_sbc(i,1)
a(1:nen_solid)=ien_sbc(i,2:nen_solid+1)
bc=ien_sbc(i,6)
if (bc .eq. 999) then
	bc=10110
else if (bc .eq. -999) then
	bc=10000
end if 
	do b=1,nen_solid
		if (b .gt. 2) then
			c=b-2
		else
			c=b+2
		endif
		if (a(b) .eq. 1) then
			if (bc .gt. 10000) then
				if ((bc .eq. 10110)) then
					ng=ng+1
					dof=nsd_solid*(c-1)+1
					p=lm_solid(dof,el)
					tmpgdof(ng)=p
					ng=ng+1
					dof=nsd_solid*(c-1)+2
					q=lm_solid(dof,el)
					tmpgdof(ng)=q
				else if (bc .eq. 10100) then
					ng=ng+1
					dof=nsd_solid*(c-1)+1
					p=lm_solid(dof,el)
					tmpgdof(ng)=p
				else if (bc .eq. 10010) then
					ng=ng+1
					dof=nsd_solid*(c-1)+2
					q=lm_solid(dof,el)
					tmpgdof(ng)=q
				endif
			endif
		endif	
	enddo
enddo
!remove nonunique entries
150 do i=1,ndof_solid
		do j=1,ndof_solid
			if ((tmpgdof(i) .eq. tmpgdof(j)) .and. (i .ne. j))then
				tmpgdof(j)=0
			endif
		enddo
	enddo
where (tmpgdof .ne. 0)
	ng_true=1
end where

nsol_ebc=sum(ng_true)
neq_solid=ndof_solid+nsol_ebc
allocate(IPIV(neq_solid))
allocate(gdof(nsol_ebc))
allocate(gx(nsol_ebc))
i=1
do b=1,ndof_solid
	if (tmpgdof(b) .ne.  0) then
		gdof(i)=tmpgdof(b)
		i=i+1
	endif
enddo
end subroutine s_ess
