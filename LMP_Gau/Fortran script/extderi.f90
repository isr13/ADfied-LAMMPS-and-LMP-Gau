!This script is edited by Indu Sekhar Roy (RWTH Aachen) for LMP_Gau integration. The original script was for gau_xtb and was written by Dr. Tian Lu at Beijing Kein Research Center for Natural Sciences (www.keinsci.com)

program extderi
implicit real*8 (a-h,o-z)
character*80 outname,c80,arg2,arg3
real*8 polar(6)
real*8,allocatable :: ddip(:),hess(:,:)
logical alive

call getarg(1,outname)
call getarg(2,arg2)
call getarg(3,arg3)
read(arg2,*) natm
read(arg3,*) ider !1=only force 2=also Hessian

open(11,file=outname,status="replace")

inquire(file="log.lammps",exist=alive)
if (alive.eqv..false.) then
	write(*,*) "Error: log.lammps is not existed in current folder!"
	stop
end if
open(10,file="log.lammps",status="old")

do while(.true.)
	read(10,"(a)") c80
	if (index(c80,"PotEng")/=0) exit
end do
read(10,*) !skipping the readin at step 0
read(10,*) stp, ene !this will read the energy at step 1
write(11,"(4D20.12)") (ene/627.5095),0D0,0D0,0D0
close(10)

inquire(file="dump.all",exist=alive)
if (alive.eqv..false.) then
	write(*,*) "Error: dump.all is not existed in current folder!"
	stop
end if
open(10,file="dump.all",status="old")
do iatm=1,9
	read(10,*)
end do
do iatm=1,natm
	read(10,*) fx,fy,fz
	write(11,"(3D20.12)") fx*(-0.52917721092/627.5095),fy*(-0.52917721092/627.5095),fz*(-0.52917721092/627.5095) !(1.88973*627.5095) (-0.529177249)/627.5095)
end do
close(10)

if (ider==2) then
	polar=0
	write(11,"(3D20.12)") polar
	allocate(ddip(9*natm))
	ddip=0
	write(11,"(3D20.12)") ddip
	allocate(hess(3*natm,3*natm))
	inquire(file="Hessian_AD.txt",exist=alive)
	if (alive.eqv..false.) then
		write(*,*) "Error: Hessian_AD.txt is not existed in current folder!"
		stop
	end if
	open(10,file="Hessian_AD.txt",status="old")
	!read(10,*) commented as there is no comment at the first line in Hessian_AD.txt
	read(10,*) ((hess(i,j),j=1,3*natm),i=1,3*natm)
	write(11,"(3D20.12)") (((hess(i,j)*(0.52917721092*0.52917721092/627.5095)),j=1,i),i=1,3*natm)!(1.88973*1.88973*627.5095) *0.529177249*0.529177249/((627.5095)
	close(10)
end if
close(11)
end program
