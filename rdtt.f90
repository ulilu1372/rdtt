subroutine rdtt (k,a,gamt,j10,wa,po,pk,lamo,muk)
real, parameter :: g=9.81
real, parameter :: pi=3.1415926
real :: k,a,gamt,j10,wa,po,pk,lamo,muk
real :: m0,mgch,m0sum,alfak,alfad,ettat
real :: mpu0,mou0,mpo0,mho0
real :: alfac,alfadn,alfas,alfadnt,alfast
real :: alfapu,alfaou,alfapo,alfaho
real :: fiz,alfatz,alfatd,alfatg
real :: B,alftz,sigmak,nu,t0,tn
real :: pkm,epse,ku,kappa,kappap,nk
real :: kkm,delkotn,deltotn,kt,dotn,u,u1,kzk
real :: L,L11,L12,Lk,Lk13
real :: kL,ka,kzm,gampr,gamk,gamtp
real :: kdn,ks,kdnt,kst
real :: dm,mt,r0,gps,fkr,dkr,ksi,fa,da,daotn,md,mk
real :: delk,dk,tauk,delt,d,d1
real :: dpp,d1pp,d1ppotn
real :: mc,mdn,ms,mdnt,mst,mct
real :: ldn,ldnp,ldnpotn,dv,ldnz,ldnzotn,ds,ef,laotn,la,lu,fi,ls,ldkr,dvh
real :: lpo,dellso,lho,lso,lgch,lgchotn,lgchk,beta,dkon,gamgch,xcd,xcm,lc,lcotn
real :: lp,x1,x3,x4,x5,x6,xcm0,xcmk,mpu,mou,mpo,mho

integer :: utoplenie
integer :: ugch





open (unit=14,file='IshDan.txt')
read (14,*) utoplenie
read (14,*) ugch
read (14,*) mgch
read (14,*) sigmak
read (14,*) gamk
read (14,*) B
read (14,*) t0
read (14,*) u1
read (14,*) nu
read (14,*) ku
read (14,*) nk
read (14,*) dotn
read (14,*) lotn
read (14,*) ka
read (14,*) ksi
read (14,*) betas
read (14,*) ef
read (14,*) fi
read (14,*) lgchotn
read (14,*) yotn
read (14,*) lcotn
close (unit=14)

sigmak=sigmak*1.0e6
betas=betas*pi/180



!utoplenie=1				!функция включения, для утопленного сопла 1, для неутопленного 0
!ugch=1					!функция включения, для конической гч 0, для цилиндро-конической 1
!mgch=1000				!масса головной части
!sigmak=1250e6			!предел прочности материала стенки, Па
!gamk=7800				!плотность материала КС, кг/м2, для стали 7800, для композита 1800
!B=200					!термохимическая константа
!t0=35					!начальная температура (градусы цельсия)
!u1=3.34					!единичная скорость горения, мм/с*(см2/кг)^nu
!nu=0.27					!степень в законе горения
!ku=1.06					!технологический разброс скорости горения, 1.03 - 1.06
!nk=1.3					!запас прочности, для сталей 1,3 и для композитов 2,5
!dotn=0.3				!относительный диаментр канала, 0.25 - 0.3
!lotn=2.5				!относительная длина заряда, 3.0 - 6.0
!ka=1.1					!коэффициент, учитывающий массу элементов крепления, для сталей 1.1 и для композитов 1.2 
!ksi=9					!степень расширения сопла по площади
!betas=17*pi/180			!угол полураствора сопла 12 - 20 град
!ef=0.3					!степень утопления сопла 0 - 0.3
!fi=0.98					!коэффициент при расчете сопла, 0.3 - 1.0
!lgchotn=3				!удлинение конической головной части, 1.5 - 3.0
!yotn=0.1				!запас статической устойчивости гч, 0.03 - 0.1
!lcotn=1.0				!относительная длина цилиндрической части гч Lц/Lгчк 0.5 - 1.0

mpu0=50					!масса приборов управления
mou0=25					!масса органов управления
mpo0=20					!масса приборного отсека
mho0=35					!масса хвостового отсека
gamtp=1500				!плотность тзп, кг/м2
tn=20					!нормальная температура заряда (градусы цельсия)
kappa=120				!параметр Победоносцева, 110 - 120
kappap=90				!пороговое значение параметра Победоносцева
kt=0.015				!скорость уноса тзп, см/с
kdn=3.5*1.0e-3			!коэффициент днища, кгс/см2
ks=0.065				!коэффициент сопла, кгс/см2
kdnt=2.5*1.0e-5			!коэффициент тзп днища, кгс/см2*с
kst=2*1.0e-5			!коэффициент тзп сопла, 1/с
alfapu=6.67*1.0e-3		!массовый коэффициент приборов управления
alfaou=7.5*1.0e-3		!массовый коэффициент органов управления
alfapo=2.41*1.0e-3		!массовый коэфициент приборного отсека
alfaho=4.04*1.0e-3		!массовый коэффициент хвостового отсека
alfatz=0.0025			!массовый коэффициент топлива заправки
alfatd=0.0045			!массовый коэффициент достартового расхода топлива
alfatg=0.025			!массовый коэффициент гарантированного запаса топлива

!термодинамический коэффициент топлива
alftz=B/(B-(t0-tn))
!эразионный коэффициент твердого топлива
epse=sqrt(kappa/kappap)
!максимальное давление в камере сгорания, атм
pkm=pk*(alftz*epse*ku)**(1/(1-nu))
!относительная толщина стенки корпуса
delkotn=nk*pkm*1.0e5/(2*sigmak)
!относительный диаметр камеры РДТТ
kkm=1-2*delkotn
!примерная скорость горения тт, см/с
u=0.1*u1*pk**nu
!относительная толщина термопокрытия
deltotn=kt*(1-dotn)/(2*u)
!коэффициент заполнения камеры
kzk=1/(1+2*deltotn)
!примерный горящий свод, см
e1=0.1*u1*pk**nu*muk*lamo*wa/g
!определение длины твердотопливного заряда (канальной шашки)
L11=po*muk/(kkm**2.0*(kzk**2.0-po*g*100/(kkm**2.0*lamo*gamt*u*wa*kappa))*gamt)
L12=2.0/3.0*(2*e1/(100*(1-dotn)))*(kzk**2.0*(2.0/3.0-dotn**2.0)*(1/(kzk**2.0-(po*g*100)/(kkm**2.0*lamo*gamt*u*wa*kappa)))-1)
L=L11-L12
!определение длины камеры сгорания
Lk13=2*e1/(100*(1-dotn))*(nk*pkm*1.0e5/sigmak+kt/u*(1-dotn))
Lk=L+Lk13
!относительная длина цилиндрической части камеры
kL=1-2.0/(3.0*lotn)
!относительный диаметр заряда
kzm=(1-2*delkotn)/(1+2*deltotn)
!приведенная плотность материала камеры
gampr=gamk+kzm*kt*(1-dotn)/(2*u)*gamtp/(delkotn*ka)
!массовый коэффициент цилиндрической части кс
alfac=4*kL*ka*L*delkotn*gampr/(po*muk)
!массовый коэффициент днища
alfadn=kdn*1.0e4/(po*muk)
!массовый коэффициент сопла
alfas=ks*10e4/(0.98*a*pkm*1.0e4*muk*lamo*j10)
!массовый коэффициент тзп днища
alfadnt=kdnt*1.0e4*lamo*j10/po
!массовый коэффициент тзп сопла
alfast=kst*j10
!массовый коэффициент двигателя
alfad=alfac+2*alfadn+alfas+2*alfadnt+alfast
!массовый коэффициент конструкции
alfak=alfapu+alfaou+alfapo+alfaho
!коэффициент заправки
fiz=1+alfatz+alfatd+alfatg
!относительная масса топлива заправки
ettat=fiz*muk/(1+alfatd*muk)

!write(*,'(7(1x,f10.4))') atan(2.8/5)*180/pi
!read(*,*)

!суммарная масса деталей кc, независящая от проектных параметров
m0sum=mpu0+mou0+mpo0+mho0	
!стартовая масса ракеты
m0=(mgch+m0sum)/(1-alfak-(1+alfad)*ettat)
!диаметр миделева сечения
dm=2*sqrt(m0/(pi*po))
!масса топлива
mt=ettat*m0
!тяга РДТТ
r0=m0/lamo
!расход продуктов сгорания
gps=m0*g/(lamo*wa)
!площадь критического сечения сопла
fkr=m0*g/(0.98*wa*lamo*a*pk*1.0e4)
dkr=sqrt(4*fkr/pi)
!площадь среза сопла
fa=ksi*fkr
da=sqrt(4*fa/pi)
daotn=da/dm
if (daotn.gt.0.75) then
	print*,'      ahtung           '
	print*,'slishkom zdorovoe soplo'
	read(*,*)
end if
!масса дигателя
md=alfad*mt
!масса конструкции изделия
mk=alfak*m0+m0sum
!толщина стенки камеры сгорания
delk=nk*pkm*1.0e5*dm/(2*(sigmak+nk*pkm*1.0e5))
!внутренний диаметр камеры
dk=dm-2*delk
!примерное время работы двигателя
tauk=muk*lamo*j10
!толщина тзп
delt=kt*tauk/100
!диаметр заряда
d=dm-2*delk-2*delt
dpp=2*muk*lamo*wa*u/(100*g*(1-sqrt(100*po*g/(gamt*u*lamo*wa*kappa*kzk**2.0))))
!диаметр канала
d1=d-2*e1/100
d1ppotn=sqrt(100*po*g/(gamt*u*lamo*wa*kappa*kzk**2.0))
d1pp=d1ppotn*dpp
!масса цилиндрической части камеры
mc=alfac*mt
!масса одного днища камеры
mdn=alfadn*mt
!масса сопла
ms=alfas*mt
!масса тзп одного днища
mdnt=alfadnt*mt
!масса тзп сопла
mst=alfast*mt
!масса тзп цилиндрической части камеры
mct=pi*kzk*dk*L*delt*gamtp
if (utoplenie==0) then
!длинна металлического эллиптического днища
ldn=0.33*dm
ldnp=ldn
ldnz=ldn
!длина сопла
ls=dkr*(sqrt(ksi)-1)/(2*tan(betas))
else
!диаметр под воспламенитель
dv=0.2*dk
!длина воспламенителя
lvos=0.1*dm
!относительная длина переднего днища
ldnpotn=0.3+0.05*(dv/dk)
!длина переднего днища
ldnp=ldnpotn*dm
!диаметр отверстия под сопло
ds=dkr*(1.5+ef*(sqrt(ksi)-1.5))
!относительная длина заднего днища
ldnzotn=0.3+0.05*(ds/dk)
!длина заднего днища
ldnz=ldnzotn*dm
laotn=3.2*fi*(sqrt(ksi))**(0.829+0.298*k**2.0)
!длина сопла
ls=(1-ef)*laotn*dkr/2.0
!длина утопленной части сопла
lu=ef*laotn*dkr/2.0
!длина сверхзвуковой части сопла
la=laotn*dkr/2.0
!длина докритической части сопла
ldkr=0.7*dkr
!диаметр входной части сопла
dvh=1.5*dkr
end if
!длина приборного отсека, в диапозоне 0,8 - 1,2
lpo=1.0*dm
!длина переходного отсека, в диапозоне 0,1 - 0,2
dellso=0.15*dm
!длина хвостового отсека
lho=ldnz+ls
!длина соединительного отсека
lso=ldnp+dellso
!расчет длины головной части

!длина головной части
lgch=lgchotn*dm

if (ugch==0) then
!длина конической части гч
lgchk=4.0/3.0*lgch*(yotn+2.0/3.0)
!угол конуса гч
beta=atan(1/(2*lgchotn))
dkon=lgchk/lgchotn
lc=0.0
!приведенная плотность головной части
gamgch=12*mgch*lgchotn**2.0/(pi*lgchk**3.0)
!расстояние до центра давления гч
xcd=2.0*lgch/3.0
!расстояние до центра масс гч
xcm=3.0*lgchk/4.0	
	else
	lgchk=lgch/(1+lcotn)
	beta=atan((1+lcotn)/(2*lgchotn))
	gamgch=12*mgch*(1+lcotn)/(pi*lgchotn*dm**3.0)
	lc=lgch-lgchk
end if

if (gamgch.lt.1500) then
	print*,'           ahtung             '
	print*,'plotnost gch menshe 1500 kg/m3'
	read(*,*)
end if

if (gamgch.gt.2500) then
	print*,'           ahtung             '
	print*,'plotnost gch bolshe 2500 kg/m3'
	read(*,*)
end if

!полная длина изделия
lp=lgch+lpo+dellso+lk+ls

!налагаемые ограничения
if (e1/(100*d).gt.0.5) then
	print*,'         ahtung              '
	print*,'ochen bolshoy goryaschiy svod'
	read(*,*)
end if

if (d1/dkr.lt.1.01) then
	print*,'     ahtung         '
	print*,'slishkom uzkiy kanal'
	read(*,*)
end if

if (lp/dm.gt.12) then
	print*,'        ahtung            '
	print*,'slishkom bolshoe udlinenie'
	read(*,*)
end if

!определение положения центра масс

!координата цм конической части гч
x1=3.0/4.0*lgchk
!координата цм приборного отсека
x3=lgch+0.5*lpo
!координата цм камеры
x4=lgch+lpo+dellso+0.5*Lk
!координата цм хвостового отсека
x5=lgch+lpo+dellso+Lk+0.5*lho
!координата цм сопла
x6=lgch+lpo+dellso+Lk+ls*(1-(da**2.0/4+0.5*da*dkr+3.0*dkr**2.0/4.0)/(4*(da**2.0+da*dkr/4.0+dkr**2.0/4.0)))

!масса приборов управления, приборного отсека, хвостового отсека и органов управления
mpu=alfapu*m0+mpu0
mpo=alfapo*m0+mpo0
mho=alfaho*m0+mho0
mou=alfaou*m0+mou0

!координата цм изделия до старта
xcm0=(x1*mgch+x3*(mpu+mpo)+x4*(md+mt)+x5*(mho+mou)+x6*ms)/m0
!координата цм изделия в конце АУТ
xcmk=(x1*mgch+x3*(mpu+mpo)+x4*md+x5*(mho+mou)+x6*ms)/(m0-mt)

!write(*,'(7(1x,f10.4))') gamgch
!read(*,*)

open (unit=12,file='trh.txt')
write (12,"('Стартовая масса изделия, кг =',2(2x,f10.3))") m0
write (12,"('Масса твердого топлива, кг =',2(2x,f10.3))") mt
write (12,"('Тяга двигателя, кг =',2(2x,f12.3))") r0
write (12,"('Расход продуктов сгорания, кг/с =',2(2x,f10.3))") gps
write (12,"('Время работы двигателя, с =',2(2x,f10.3))") tauk
write (12,"('Расчетный горящий свод, м =',2(2x,f10.3))") e1/100
write (12,"('Диамертр миделева сечения, м =',2(2x,f10.3))") dm
write (12,"('Полная длина изделия, м =',2(2x,f10.3))") lp
write (12,"('Относительное удлинение =',2(2x,f10.3))") lp/dm
close (unit=12)

open (unit=13,file='mgh.txt')
write (13,"('Стартовая масса изделия, кг =',2(2x,f10.3))") m0
write (13,"('Масса твердого топлива, кг =',2(2x,f10.3))") mt
write (13,"('Масса головной части, кг =',2(2x,f10.3))") mgch
write (13,"('Приведенная плотность гч =',2(2x,f10.3))") gamgch
write (13,"('Масса приборов управления, кг =',2(2x,f10.3))") mpu
write (13,"('Масса приборного отсека, кг =',2(2x,f10.3))") mpo
write (13,"('Масса пустого двигателя, кг =',2(2x,f10.3))") md
write (13,"('Масса снаряженного двигателя, кг =',2(2x,f10.3))") md+mt
write (13,"('Масса конструкции носителя, кг =',2(2x,f10.3))") mk
write (13,"('Масса сопла, кг =',2(2x,f10.3))") ms
write (13,"('Масса хвостового отсека, кг =',2(2x,f10.3))") mho
write (13,"('Масса органов управления, кг =',2(2x,f10.3))") mou
write (13,"('Диаметр миделева сечения, м =',2(2x,f10.3))") dm
write (13,"('Наружный диаметр заряда, м =',2(2x,f10.3))") d
write (13,"('Расчетный диаметр канала, м =',2(2x,f10.3))") d1
write (13,"('Диаметр критического сечения, м =',2(2x,f10.3))") dkr
write (13,"('Диаметр среза сопла, м =',2(2x,f10.3))") da
write (13,"('Толщина стенки камеры сгорания, м =',2(2x,f10.3))") delk
write (13,"('Толщина ТЗП, м =',2(2x,f10.3))") delt
write (13,"('Полная длина изделия, м =',2(2x,f10.3))") lp
write (13,"('Относительное удлинение =',2(2x,f10.3))") lp/dm
write (13,"('Длина головной части, м =',2(2x,f10.3))") lgch
write (13,"('Длина конуса головной части, м =',2(2x,f10.3))") lgchk
write (13,"('Длина приборного отсека, м =',2(2x,f10.3))") lpo
write (13,"('Длина соединительного отсека, м =',2(2x,f10.3))") lso
write (13,"('Расстояние между приборным отсеком и передним днищем, м =',2(2x,f10.3))") dellso
write (13,"('Длина переднего днища, м =',2(2x,f10.3))") ldnp
write (13,"('Длина заряда, м =',2(2x,f10.3))") L
write (13,"('Длина камеры сгорания, м =',2(2x,f10.3))") Lk
write (13,"('Длина заднего днища, м =',2(2x,f10.3))") ldnz
write (13,"('Длина сопла, м =',2(2x,f10.3))") ls
write (13,"('Длина сверхзвуковой части утопленного сопла, м =',2(2x,f10.3))") la
write (13,"('Длина докритической части утопленного сопла, м =',2(2x,f10.3))") ldkr
write (13,"('Длина утопленной части сопла, м =',2(2x,f10.3))") lu
write (13,"('Диаметр входной части утопленного сопла, м =',2(2x,f10.3))") dvh
close (unit=13)

200 end subroutine rdtt

















