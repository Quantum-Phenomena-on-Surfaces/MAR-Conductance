module declarations

  implicit none

! PARAMETERS

  integer, parameter :: kk = SELECTED_INT_KIND (10)
  integer, parameter :: q = SELECTED_REAL_KIND(10)
  integer, parameter :: qs = SELECTED_REAL_KIND(5)
  integer, parameter :: qc = SELECTED_REAL_KIND(10)
  real (q), parameter :: pi_d = 3.14159265358979323846_q
  real (q), parameter :: sqrtpi_d = 1.77245385090551602729_q
  real (q), parameter :: Hartree = 27.2116_q
  complex (qc), parameter :: ui = (0._q,1._q)
  complex (qc), parameter :: ur = (1._q,0._q)
  complex (qc), parameter :: zero = (0._q,0._q)

! code conventions and parameters
  integer, parameter :: sg_R=1,sg_L=-1 !if you touch this, be careful with
! sums over Floquet indices!!
! L electrode is the first index and right electrode the second index


   character (100) :: Name_output
   character (100) :: DName_output

! numbers
  integer :: i, j, n, m, n1, i_V
  integer :: N_Bias, N_Floquet, N_e, i_e

  real (q) :: S,  Delta_R, Delta_L, J_, U, Temp
  real (q) :: Gamma_0,  Dynes_L, Dynes_R
  real (q) :: Bias_ini, Bias_fin, step_BIAS, V, omega
  real (q) ::  theta_i, phi_i, ed_up, ed_down, phi_L, phi_R 
  real (q) :: curr_0_L, curr_1, curr_2, e_x, phi_, bandwidth, curr_avant

!arrays
  complex (qc), allocatable :: S_l (:,:,:,:,:), S_r (:,:,:,:,:), S_a (:,:,:,:,:)
  complex (qc), allocatable :: S_M (:,:,:,:,:)
  complex (qc), allocatable :: G_l (:,:,:,:,:), G_r (:,:,:,:,:), G_a (:,:,:,:,:)
  complex (qc), allocatable :: gr (:,:,:,:,:), ga (:,:,:,:,:)

end module declarations
