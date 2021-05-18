program optim_pp
real, parameter :: pi=3.1415926
real, parameter :: g=9.81
real, parameter :: patm0=101325
real, parameter :: t0=288
real, parameter :: ro0=1.225
real, parameter :: gkz=9.986e14
real, parameter :: rpl=6.371e6
real :: x,k,Rt,T
real :: vk,po,pk,lamo,muk
real :: wa,j10,dj1,dj2
real :: h,his,vis,vtek,dvdmu,dhdmu
real :: patm,th,roh,gh,ah,lamoh
real :: mut,mah,st,tet,cx
real :: r,faer,fgr,fvis
real :: dv,dvo,poo,pko,lamoo,muko,ho,vo
real :: kj1,dzeta,pia,Ac
real :: j10v,gamt,qgg
integer :: delta
open (unit=10,file='opt_ish.txt')
open (unit=11,file='opt_rez.txt')
read (10,*) x
read (10,*) k
read (10,*) Rt
read (10,*) T
read (10,*) j1
read (10,*) gamt
read (10,*) delta
read (10,*) st
close (unit=10)
vk=11.19*sqrt(tan(pi*x/(180*222.4))*tan(pi/180*(45-x/444.8)))*1000
dvo=vk
do po=5000,10000,500
do pk=60,120,10.0
wa=sqrt(2*k/(k-1)*Rt*T*g)*sqrt(1-(0.75/pk)**((k-1)/k))
!j10=sqrt(2*k/(k-1)*Rt*T*g)*sqrt(1-(0.75/pk)**((k-1)/k))+11.22*sqrt(g*Rt*T/k*((k+1)/2)**((k+1)/(k-1)))*(0.75/pk-1/pk)
!j10=sqrt(2*k/(k-1)*Rt*T*g)*sqrt(1-(0.75/pk)**((k-1)/(k+1)))
!wa=sqrt(2*k/(k-1)*Rt*T*g)*sqrt(1-(0.75/pk)**((k-1)/k))+11.22*sqrt(g*Rt*T)/sqrt(k*(2/(k-1)**((k+1)/(k-1))))*(0.75/pk-1/pk)
pia=0.75/pk
Ac=sqrt(g*k/(Rt*T))*sqrt((2/(k+1))**((k+1)/(k-1)))
do lamo=0.1,0.4,0.05
do muk=0,0.5,0.05
mut=1.0
vtek=0
h=0
140 if (h.lt.11e3) then
	patm=patm0*(1-h/44300)**5.256
	th=t0*(1-h/44300)
	roh=ro0*(1-h/44300)**4.256
	else
		patm=22690*exp(-1*(h-11000)/6340)
		th=216.5
		roh=0.365*exp(-1*(h-11000)/6340)
	end if
gh=gkz/(rpl+h)**2.0
lamoh=lamo*gh/g
if (mut.gt.0.95) then 
	tet=pi/2
end if
if ((mut.gt.0.45).and.(mut.lt.0.95)) then
	tet=(4*(pi/2-pi/4)*(mut-0.45)**2.0)+pi/4
end if
if (mut.lt.0.45) then 
	tet=pi/4
end if
ah=20*sqrt(th)
mah=vtek/ah
if (mah.lt.0.8) then 
	cx=0.29
end if
if ((mah.gt.0.8).and.(mah.lt.1.068)) then 
	cx=mah-0.51
end if
if (mah.gt.1.068) then 
	cx=0.091+0.5/mah
end if


kj1=sqrt((1-pia*((k-1)/k))/(1-(patm/(pk*10**5.0))**((k-1)/k)))
dzeta=(2/(k+1))**(1/(k-1))*sqrt((k-1)/(k+1))/sqrt(pia**(2/k)-pia**((k+1)/k))
qgg=0.7*pk/(gamt*1.0e4)
j10=(kj1*j1+dzeta/(0.98*Ac))*0.98*0.96*0.96/(1+qgg*delta)
!write(*,'(7(1x,f10.3))') j10,j10v
!read(*,*)


dj1=j10/wa
dj2=dj1-1
r=wa*dj1/mut
faer=wa*lamoh*cx*roh*vtek**2.0/(2*g*po*mut)
fgr=wa*lamoh*sin(tet)
fvis=wa*dj2*patm/(patm0*mut)
dvdmu=r-faer-fgr-fvis
vtek=vtek+dvdmu*st
dhdmu=wa*lamoh*vtek*sin(tet)/gh
h=h+dhdmu*st

if (mut.gt.muk) then
	mut=mut-st
	goto 140
	else
	goto 150
end if
150 vis=vtek
his=h
call sr
end do
end do
end do
end do
write (11,"('Потребное значение скорости, м/с =',2(2x,f10.2))") vk
write (11,"('Оптимальное значение скорости, м/с =',2(2x,f10.2))") vo
write (11,"('Расчетное значение удельного импульса J10, с =',2(2x,f10.2))") j10o/9.81
write (11,"('Расчетное значение скорости истечения Wa, m/с =',2(2x,f10.2))") wa
write (11,"('Оптимальное значение По, кг/м2 =',2(2x,f10.0))") poo
write (11,"('Оптимальное значение pk, атм =',2(2x,f10.0))") pko
write (11,"('Оптимальное значение lam =',2(2x,f10.3))") lamoo
write (11,"('Оптимальное значение muk =',2(2x,f10.3))") 1-muko


close (unit=11)
muko=1-muko
j10=j10o/9.81
call rdtt (k,Ac,gamt,j10,wa,poo,pko,lamoo,muko)


stop
 
CONTAINS
subroutine sr
dv=abs(vk-vis)
if (dv.gt.dvo) then 
goto 200
end if
j10o=j10
dvo=dv
poo=po
pko=pk
lamoo=lamo
muko=muk
ho=his
vo=vis
200 end subroutine sr
end program optim_pp 


