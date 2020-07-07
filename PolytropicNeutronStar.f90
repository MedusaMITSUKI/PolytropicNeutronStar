program polytrope
  implicit none
  ! 定数設定 (SI単位系)
  double precision,parameter :: pi = atan(1.0d0) * 4.0d0  ! 円周率
  double precision,parameter :: G = 6.67408d-11           ! 万有引力定数
  double precision,parameter :: c = 2.99792458d+8         ! 光速定数
  double precision,parameter :: hbar = 1.0545718d-34      ! プランク定数
  double precision,parameter :: m_n = 1.674927471d-27     ! 核子の質量
  double precision,parameter :: Msun = 1.98892d+30        ! 太陽質量

  ! 状態方程式(EOS)の定数設定
  double precision,parameter :: Gamma0 = 5.0d0 / 3.0d0    ! simeq ポリトロピック指数
  double precision,parameter :: K0 = ( 3.0d0 * pi ** 2.0d0 ) ** ( 2.0d0 / 3.0d0 ) * hbar ** 2.0d0 & 
    / ( 5.0d0 * m_n ** ( 8.0d0 / 3.0d0 ) )                ! 比例定数
  double precision,parameter :: Gamma1 = 3.0d0            ! simeq ポリトロピック指数 (2 ~ 4の値に変更してみよう！）
  double precision,parameter :: rho1 = 5d+17              ! 低密度・高密度領域の境目の密度
  double precision,parameter :: P1 = K0 * rho1 ** Gamma0  ! 低密度・高密度領域の境目のEOS
  double precision,parameter :: K1 = P1 / rho1 ** Gamma1  ! 比例定数

  ! 初期値/変数設定
  double precision,parameter :: dr = 1.0d0                ! 星の初期半径
  double precision :: r                                   ! 星の半径
  double precision :: m                                   ! 星の質量
  double precision :: P                                   ! 星の圧力
  double precision :: rhoc                                ! 星の中心密度

  ! Runge-Kutta (4th Order)用変数設定
  double precision,parameter :: nmin = 17.5d0             ! 中心密度の冪(最小)
  double precision,parameter :: nmax = 20.0d0             ! 中心密度の冪(最大)
  double precision :: n                                   ! 中心密度の冪
  double precision :: dn                                  ! 中心密度の冪の変化量
  double precision :: k(4)                                ! 係数
  double precision :: j(4)                                ! 係数
  integer :: i                                            ! ループカウンタ
  integer, parameter :: imin = 1                          ! ループカウンタ(最小)
  integer, parameter :: imax = 1000                       ! ループカウンタ(最大)
  integer, parameter :: istep = 1                         ! ループカウンタステップ

  ! 計算結果ファイル設定
  open (1, file='polytrope.dat', status='unknown', form='formatted')

  do i = imin, imax, istep

    ! 初期値設定
    dn = (nmax - nmin) * (dble(i) / dble(imax))
    n = nmin + dn
    rhoc = 10.0d0 ** n ! 星の中心密度
    r = 1.0d-6         ! 初期半径
    P = pressure(rhoc) ! 圧力
    m = 0d0            ! 星の初期質量

    ! Runge-Kutta
    do while ( P > 0d0 )
      k(1) = dr * dPdr( r, P, m )
      k(2) = dr * dPdr( r + dr / 2.0d0, P + k(1) / 2.0d0, m ) 
      k(3) = dr * dPdr( r + dr / 2.0d0, P + k(2) / 2.0d0, m ) 
      k(4) = dr * dPdr( r + dr, P + k(3), m ) 

      j(1) = dr * dmdr( r, m, P )
      j(2) = dr * dmdr( r + dr / 2.0d0, m + j(1) / 2.0d0, P ) 
      j(3) = dr * dmdr( r + dr / 2.0d0, m + j(2) / 2.0d0, P ) 
      j(4) = dr * dmdr( r + dr, m + j(3), P )  

      m = m + ( j(1) + 2.0d0 * j(2) + 2.0d0 * j(3) + j(4) ) / 6.0d0 
      P = P + ( k(1) + 2.0d0 * k(2) + 2.0d0 * k(3) + k(4) ) / 6.0d0
      r = r + dr
    end do

    ! 計算結果書き込み
    !! 出力は M [Solar mass], R [km], rhoc [kg m^{-3}]
    write(1,*) m/Msun, r*1.0d-3, rhoc
  end do

  close(1)

  ! グラフ描画
  call setting_tov_gnuplot()
  call execute_command_line('gnuplot "polytrope.plt"')

contains

  ! 状態方程式 = (圧力 P)
  double precision function pressure (rho)
    double precision,intent(in) :: rho
    if (rho < rho1) then
      pressure = K0 * rho ** Gamma0
    else 
      pressure = K1 * rho ** Gamma1
    end if
    return
  end function pressure

  ! 状態方程式 (密度 rho)
  double precision function density (P)
    double precision,intent(in) :: P
    if ( P < P1 ) then
      density = ( P / K0 ) ** ( 1.0d0 / Gamma0 )
    else
      density = ( P / K1 ) ** ( 1.0d0 / Gamma1 ) 
    end if
    return
  end function density

  ! TOV方程式
  double precision function dPdr (r, P, m)
    double precision,intent(in) :: r, P, m
    double precision :: rho
    rho = density(P)

    dPdr = -G * ( rho + P / c ** 2.0d0 ) * ( m + 4.0d0 * pi * r ** 3.0d0 * P / c ** 2.0d0 )
    dPdr = dPdr / ( r ** 2.0d0 * ( 1.0d0 - 2.0d0 * G * m / ( r * c ** 2.0d0 ) ) )
    return
  end function dPdr

  ! TOV方程式 
  double precision function dmdr (r, m, P)
    double precision,intent(in) :: r, m, P
    double precision :: rho
    rho = density(P)

    dmdr = 4.0d0 * pi * r ** 2.0d0 * rho
    return
  end function dmdr

  ! グラフ描画 (Gnuplot)
  subroutine setting_tov_gnuplot()
    open (5, file = 'polytrope.plt', status = 'unknown', form = 'formatted')
    write (5, '(a)') 'set terminal qt enhanced font "Arial, 20" dashed dashlength 1'
    write (5, '(a)') 'set format "%g"'
    write (5, '(a)') 'set xlabel "Radius [km]"'
    write (5, '(a)') 'set ylabel "Mass [Solar mass]"'
    write (5, '(a)') 'set mxtics 2'
    write (5, '(a)') 'set mytics 2'

    ! M, Rが0に近い値の場合はグラフから省く
    write (5, '(a)') 'plot "polytrope.dat" using ($2 <= 1E-3? 1/0 : $2):($1 <= 1E-3? 1/0 : $1) with lines lw 3 notitle"'
    write (5, '(a)') 'pause -1'
  end subroutine setting_tov_gnuplot

end program polytrope
