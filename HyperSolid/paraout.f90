subroutine paraout(t,dnew,vel,solid_con,sel,fext)
use solid_variables 
implicit none
integer :: i,t,solid_con(ne_solid,nen_solid)
real(8) :: dnew(ndof_solid)
real(8) :: vel(ndof_solid)
real(8) :: sel(3,ne_solid)
real(8) :: fext(ndof_solid)
character(len=12) :: flnm
character(len=5) :: x1
character(len=5) :: x2

write (x1,20) t
flnm='pt1'//x1//'.vtk'
write(*,*) flnm
open(unit=22,file=flnm)
write(22,80) '# vtk DataFile Version 2.0'
write(22,80)'solid geometry'
write(22,80)'ASCII'
write(22,80)'DATASET POLYDATA'
write(22,90)'POINTS',nn_solid,'double'
do i=1,nn_solid
	write(22,70) dnew(nsd_solid*(i-1)+1)+xref_solid(nsd_solid*(i-1)+1), dnew(nsd_solid*(i-1)+2)+xref_solid(nsd_solid*(i-1)+2),0.0d0
enddo
write(22,95)'POLYGONS',ne_solid,ne_solid*5
do i=1,ne_solid
	write(22,30) 4,solid_con(i,1)-1,solid_con(i,2)-1,solid_con(i,3)-1,solid_con(i,4)-1
enddo
write(22,96)'CELL_DATA',ne_solid
write(22,97)'SCALARS','sxx','float',1
write(22,80)'LOOKUP_TABLE default'
do i=1,ne_solid
	write(22,98) sel(1,i)
enddo
write(22,97)'SCALARS','syy','float',1
write(22,80)'LOOKUP_TABLE default'
do i=1,ne_solid
	write(22,98) sel(2,i)
enddo
write(22,97)'SCALARS','sxy','float',1
write(22,80)'LOOKUP_TABLE default'
do i=1,ne_solid
	write(22,98) sel(3,i)
enddo
write(22,96)'POINT_DATA',nn_solid
write(22,97)'VECTORS','vel','float'
do i=1,ndof_solid,nsd_solid
	write(22,70) vel(i),vel(i+1),0.0d0
enddo
write(22,97)'VECTORS','fext','float'
do i=1,ndof_solid,nsd_solid
	write(22,70) fext(i),fext(i+1),0.0d0
enddo
70 format(f10.4,f10.4,f10.4) 
20 format(I5.5) 
30 format(I5,I5,I5,I5,I5) 
80 format(A)
90 format(A,x,i5,x,A)
95 format(A,x,i5,x,i5)
96 format(A,x,i5)
97 format(A,x,A,x,A,x,i5)
98 format(f10.1)
99 format(A,x,A,x,A)
end subroutine paraout