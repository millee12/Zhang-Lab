! this module is used to apply boundary conditions
module form
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formd(ds,rngface,ien)
  use global_constants
  use fluid_variables, only: nn,ne,ndf,nsd,nen,neface,nrng,nnface,mapping,bc,bv,etype,ic,static,udf,vdf,wdf, maxconn

  implicit none

  integer  rngface(neface,ne),ien(nen,ne)
  real(8) :: ds(ndf,nn)
  integer :: idf, inl, iec, irng, ieface, inface, inn
  real(8) :: eps1,eps2
  real(8) :: hs(nrng+1,nn), h(nrng+1,nn)

  eps1 = -1000000.0 
  eps2 = -10000.0 
  h(:,:) = 0.0d0
  ds(:,:) = eps1
  
  do ieface=1,neface
     do inface=1,nnface
        inl = mapping(ieface,inface,etype)
        do iec=1,ne
           irng = rngface(ieface,iec)
           if (irng.ne.0) h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0
        enddo
     enddo
  enddo

  hs = h

  do irng=1,nrng
     do inn=1,nn
        do idf=1,ndf
           if (hs(irng,inn).gt.1.0e-8) then
              if (bc(idf,irng) .gt. 0) then
                ds(idf,inn) = bv(idf,irng)
			  endif
           endif
        enddo
     enddo
  enddo


  do inn=1,nn
     do idf=1,ndf
        if(ds(idf,inn).lt.eps2) then
			ds(idf,inn) = ic(idf)
		endif
     enddo
  enddo

  if(static) ds(:,:)=0.0

  return
end subroutine formd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formid(ids, rngface, ien)
  use fluid_variables

  implicit none

  integer :: ids(ndf,nn),rngface(neface,ne),ien(nen,ne)
  integer :: idf, inl, iec, irng, ieface, inface, inn,k
  real(8) :: epsr,epsl
  real(8) :: ds(ndf,nn),d(ndf,nn)

  d(1:ndf,1:nn) = 0.0d0

  epsr = 0.0001
  epsl = 0.000001

  do ieface=1,neface
     do inface=1,nnface
        inl = mapping(ieface,inface,etype)
        do iec=1,ne
           irng = rngface(ieface,iec)
           if(irng.ne.0) then
              do idf = 1,ndf
                 if(d(idf,ien(inl,iec)).lt.epsr) then
					d(idf,ien(inl,iec)) = bc(idf,irng)+epsl
					
				  endif
              enddo
           endif
        enddo
     enddo
  enddo

  ds = d
  ids = ds

  if(static) ids(:,:) = 1

  nq = 0

  do inn=1,nn
     do idf=1,ndf
        if(ids(idf,inn).eq.0) then
           nq = nq + 1
           ids(idf,inn) = nq
        else
           ids(idf,inn) = 0
        endif
     enddo
  enddo

  return
end subroutine formid

end module form
