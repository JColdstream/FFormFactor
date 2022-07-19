program structure_factor
use omp_lib

implicit none
include 'form_factor.inc'
integer(8):: iter

frame = 1
call start

do while ( trajstatus .eq. 0 .and. (lastframe .eq. 0 .or. frame .le. lastframe) )
  !if (mod(frame, 50) .eq. 0) write(*,*) frame
  write(*, *) frame
  call readheader
  call readcoordinates
  if (calc_com .eq. 1) then
    call initialcog
    call iteratecog
    call cogstructure
  elseif (calc_com .eq. 0) then
    call structure
  else
    write(*, *) 'Centre of mass calculation must take an input value of 0 or 1.'
    call exit
  endif
  frame = frame+1
enddo

! test .xyz file to check molecule rebuilding
open(unit = 110, file='test.xyz', status='unknown')
do iter = 1, nanalyse
  write(110, *) 'C', atom(:, iter)
enddo

write(*, *) 'Final frame: ', frame-1
call structure_output

contains

subroutine start

call get_command_argument(1, inputfile)
if (len_trim(inputfile) == 0) then
  write(*,*) 'Provide an input file.'
  call exit
endif

call readinput

! open lammps trajectory file
open(11, file = trajfile, status='old')
! get ntotal
call readheader

! allocate arrays
call readscatteringlengths
call nanalyse_setup
allocate(atom(3, nanalyse))
allocate(atomtype(nanalyse))
call wavevector_setup
rewind(11)

call skipframe
frame = frame+nskip

end subroutine start

! reads analysis.input input file
subroutine readinput
open(10, file = inputfile, status='old')
read(10, *)
read(10, *)
read(10, '(a)') trajfile
read(10, *) scatteringlengths
read(10, *) calc_com
read(10, *) nskip
read(10, *) lastframe
read(10, *) ntypes
read(10, '(a)') strtypesanalyse
typesanalyse = string_to_integers(strtypesanalyse, ",")
read(10, *) logqmin
read(10, *) logqmax
read(10, *) nq
write(*,*) trajfile


end subroutine readinput

! reads the scattering lengths from file given in the input file.
subroutine readscatteringlengths
  integer(sp):: i, j

  allocate(b(ntypes))
  ! reads scattering lengths input
  if (scatteringlengths .eq. '0') then
    do i = 1, ntypes
      b(i) = 1.0_dp
    enddo
  else
    open(12, file = scatteringlengths, status='old')
    read(12, *)
    read(12, *)
    do i = 1, ntypes
      read(12, *) j, b(j)
    enddo
  endif
end subroutine readscatteringlengths

! skips frames that don't need to be analysed
subroutine skipframe
integer:: i
do i = 1, nskip*(ntotal+9)
  read(11, *)
enddo
write(*,*) "SKIPPED FRAMES : ", nskip
write(*,*)
end subroutine skipframe

! reads in box length, timestep and natoms from the LAMMPS trajectory headers
subroutine readheader
read(11, *, iostat = trajstatus) headertext(1)   
read(11, *, iostat = trajstatus) timestep
read(11, *, iostat = trajstatus) headertext(2)
read(11, *, iostat = trajstatus) ntotal
read(11, *, iostat = trajstatus) headertext(3)
read(11, *, iostat = trajstatus) xmin, xmax
read(11, *, iostat = trajstatus) ymin, ymax
read(11, *, iostat = trajstatus) zmin, zmax
read(11, *, iostat = trajstatus) headertext(4)

!if(natom  .ne. ntotal) stop 'mismatched number of atoms'

lx = xmax-xmin
ly = ymax-ymin
lz = zmax-zmin

! Edit pbc conditions if it needs to be used for a non-cubic box.
if ((lx .ne. ly) .or. (ly .ne. lz)) then
  write(*, *) 'Script is not suitable for non-cubic boxes.'
  call exit
endif

!write(*,*) timestep
!write(*,*) xmin, xmax
!write(*,*) ymin, ymax
!write(*,*) zmin, zmax

end subroutine readheader

subroutine nanalyse_setup
  integer(sp):: i, atom_type, nindex
  nanalyse = 0
  do i = 1, ntotal
    read(11, *) nindex, atom_type
    if (atom_type .gt. ntypes) then
      write(*,*) "Number of atom types mismatch." 
      call exit
    endif
    if ( any(typesanalyse .eq. atom_type) ) then
      scattering_sum = scattering_sum+b(atom_type)
      nanalyse = nanalyse+1
    endif
  enddo
  ! array to keep track of which atoms we need to calculate
  allocate(molanalyse(nanalyse))
end subroutine nanalyse_setup

subroutine readcoordinates
  integer(sp):: i, j, nindex, temp_atom_type
  real(dp):: tempx, tempy, tempz
  j = 1
  do i = 1, ntotal
    read(11, *, iostat = trajstatus)  nindex, temp_atom_type, tempx, tempy, tempz

    if (any(typesanalyse .eq. temp_atom_type)) then
      atom(:, j) = (/tempx, tempy, tempz/)
      atomtype(j) = temp_atom_type
      j = j+1
    endif
  enddo

 ! zero coordinates
 atom(1, :) = atom(1, :)-xmin
 atom(2, :) = atom(2, :)-ymin
 atom(3, :) = atom(3, :)-zmin
end subroutine readcoordinates

subroutine wavevector_setup
  integer(8):: iq
  real(dp):: logq
  allocate(q(nq+1))
  allocate(p(nq+1))
  q = 0.0_dp
  p = 0.0_dp
  do iq = 0, nq
    logq = logqmin + (logqmax-logqmin)*iq/nq
    q(iq+1) = q(iq+1) + 10**logq
  enddo
end subroutine wavevector_setup

subroutine initialcog
  integer(8):: mol, i

  cog = 0.0_dp
  do mol = 1, nanalyse
    do i = 1, 3
      cog(i) = cog(i) + atom(i, mol) 
    enddo
  enddo
  cog = cog/nanalyse
end subroutine initialcog

subroutine iteratecog
  integer(8):: i, mol
  real(dp), dimension(3):: temp_cog, dxyz, cog_diff
  real(dp):: cog_diff_mag 
  logical:: final_cog

  final_cog = .False.

  do while (final_cog .eqv. .False.)
    temp_cog = cog
    cog = 0.0_dp
    do mol = 1, nanalyse
      do i = 1, 3
        ! dxyz(i) = dxyz(i) + (lx*anint(dxyz(i)/lx))
        cog(i) = cog(i) + atom(i, mol) - lx*anint((atom(i, mol)-temp_cog(i))/lx)
      enddo
    enddo
    cog = cog/nanalyse
    do i = 1, 3
      if (cog(i) .gt. lx) then
        cog(i) = cog(i) - lx
      elseif (cog(i) .lt. 0.0) then
        cog(i) = cog(i) + lx
      endif
    enddo
    ! write(*,*) cog
    cog_diff = cog-temp_cog
    ! check how close the new centre of geometry is to the previous iteration
    cog_diff_mag = 0.0_dp
    do i = 1, 3
      cog_diff_mag = cog_diff_mag+cog_diff(i)**2
    enddo
    if (cog_diff_mag .eq. 0.0_dp) then
      final_cog = .True.
    endif
  enddo
end subroutine iteratecog

subroutine cogstructure
  integer(8):: i, iq, mol, moli, molj
  real(dp):: tempq, drsq, dr, drcog, qrij
  real(dp), dimension(3):: dxyz
  rg = 0.0_dp
  do mol = 1, nanalyse
    do i = 1, 3
      drcog = atom(i, mol) - cog(i)
      if ( abs(drcog) .gt. lx/2 ) then
        atom(i, mol) = atom(i, mol) - lx*anint(drcog/lx)
      endif
      rg = rg + (atom(i, mol) - cog(i))**2
    enddo
  enddo
  write(*, *) 'Rg =', dsqrt(rg/nanalyse)

  !$OMP PARALLEL DO &
  !$OMP SCHEDULE(DYNAMIC) &
  !$OMP DEFAULT(NONE) &
  !$OMP SHARED(nq, q, lx, atom, nanalyse, b, atomtype, rg) &
  !$OMP PRIVATE(iq, moli, molj, drsq, dr, qrij, dxyz, i, tempq) &
  !$OMP REDUCTION(+:p)
  do iq = 1, nq+1
    tempq = 0.0_dp
    do moli = 1, nanalyse
      do molj = 1, nanalyse
        if (moli .eq. molj) then
          tempq = tempq+b(atomtype(moli))**2
        else
          drsq = 0.0_dp
          dxyz = atom(:, molj) - atom(:, moli)
          do i = 1, 3
            drsq = drsq+dxyz(i)**2
          enddo
          dr = dsqrt(drsq)
          qrij = q(iq)*dr
          tempq = tempq+b(atomtype(moli))*b(atomtype(molj))*sin(qrij)/qrij
        endif
      enddo
    enddo
    p(iq) = p(iq) + tempq
  enddo
  !$OMP END PARALLEL DO
  pqnrm = pqnrm+1.0_dp
          
end subroutine cogstructure

subroutine structure
  integer(8):: i, iq, moli, molj
  real(dp):: tempq, drsq, dr
  real(dp), dimension(3):: dxyz
  ! real(dp), external:: ddot
  !$OMP PARALLEL DO &
  !$OMP SCHEDULE(GUIDED) &
  !$OMP DEFAULT(NONE) &
  !$OMP SHARED(nq, q, lx, atom, nanalyse, b, atomtype) &
  !$OMP PRIVATE(iq, moli, molj, drsq, dr, dxyz, i, tempq) &
  !$OMP REDUCTION(+:p)
    do iq = 1, nq+1
      tempq = 0.0_dp
      do moli = 1, nanalyse
        do molj = moli, nanalyse 
          if (moli .eq. molj) then
            tempq = tempq+b(atomtype(moli))**2
          else
            dxyz = atom(:, molj) - atom(:, moli)
            drsq = 0.0_dp
            ! magnitude of the vector
            do i = 1, 3
              if ( abs(dxyz(i)) .gt. lx/2 ) then
                dxyz(i) = dxyz(i) - lx*anint(dxyz(i)/lx)
              endif
              drsq = drsq+dxyz(i)**2 
            enddo
            dr = dsqrt(drsq)
            tempq = tempq+b(atomtype(moli))*b(atomtype(molj)) * sin(q(iq)*dr)/q(iq)*dr
          endif
        enddo
      enddo
      p(iq) = p(iq) + tempq
    enddo 
  !$OMP END PARALLEL DO
  ! keep track of how many times we have done the structure factor sum
  ! sk will be normalised in the output
  pqnrm = pqnrm+1.0
end subroutine structure

subroutine structure_output
  integer(4):: iq
  ! Writes q and p(q) to file
  open(unit = 101, file='formfactor.out', status='unknown')
  do iq = 1, nq+1
    write(101, *) q(iq), p(iq)/dble(scattering_sum)**2/pqnrm
  enddo
end subroutine structure_output

function string_to_integers(str, sep) result(a)
    integer, allocatable:: a(:)
    character(*):: str
    character:: sep
    integer:: i, n_sep

    n_sep = 0
    do i = 2, len_trim(str)
      if (str(i:i)==sep .and. str(i-1:i-1)/=sep) then
        n_sep = n_sep+1
        str(i:i) = ','
       end if
    end do
    allocate(a(n_sep+1))
    read(str, *) a
end function

end program structure_factor
