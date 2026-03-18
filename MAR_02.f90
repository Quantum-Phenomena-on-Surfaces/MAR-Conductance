program currentSC
!
!
! This program computes the time-dependent SC current
! and its DC average
!
!
! couplings are exponential in energy to help
! convergence in Floquet indices
!        Gamma_(i_e)=Gamma_0*exp(-(e/(4*Bandwidth))**2)
!
! Coded by N. Lorente
! Juin 2020
!
Use declarations
Use GreenFunctions
Use SelfEnergies
use Tools
Use current
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Entrée (le fichier doit exister)
open (1,file='MAR.in',status='old')

    read (1,*) S ! Molecular spin
    read (1, *) theta_i, phi_i ! angles of spin in deg
    read (1, *) Gamma_0!  interelectrode coupling in eV
    read (1, *) Bandwidth!  Bandwidth in eV
    read (1, *) Delta_L! Left SC gap in eV (this is the sample)
    read (1, *) Delta_R! Right SC gap in eV
    read (1, *) Dynes_L! Dynes' Broadening in eV
    read (1, *) Dynes_R
    read (1, *) J_ ! exchange coupling in eV
    read (1, *) U ! potential scattering in eV
    read (1,*) Bias_ini ! in volts
    read (1,*) Bias_fin ! in volts
    read (1,*) N_Bias ! sampling in bias
    read (1,*) N_Floquet ! Cutoff in Floquet summation
!! we extend to negative numbers too
         N_Floquet=2*N_Floquet + 1 !including zero
    read (1,*) N_e ! Sampling in energy integration
    read (1,*) Temp ! in Kelvin
    read (1,*) phi_L, phi_R !phases in deg
    read (1,*) Name_output ! character
    read (1,*) DName_output ! character

close (1)
open (2, file=Name_output)
open (3, file=DName_output)
      write (3,*) '#Bias (mV), Conductance in Go' 

 

! change units

      !Angles from deg to rad:
       theta_i = theta_i*pi_d/180._q; phi_i=phi_i*pi_d/180._q
       phi_L = phi_L*pi_d/180._q; phi_R = phi_R*pi_d/180._q
       phi_=phi_L-phi_R
               !Energies from eV to Hartree
               Bandwidth=Bandwidth/Hartree ! a.u.
               Gamma_0=Gamma_0/Hartree !a.u.
               Delta_L=Delta_L/Hartree;Delta_R=Delta_R/Hartree
               Dynes_L=Dynes_L/Hartree;Dynes_R=Dynes_R/Hartree
               J_ = J_ /Hartree; U= U/Hartree ! a.u.
               Bias_ini=Bias_ini/Hartree; Bias_fin=Bias_fin/Hartree
          !Temp: Kelvin to Hartree
           Temp = Temp * 25.8_q/(1000*Hartree*300._q)
!Loop on BIAS goes from Bias_ini to Bias_fin in N_Bias steps:
   step_Bias=(Bias_fin-Bias_ini)/(N_Bias-1)

!fin d'entrée
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ALLOCATIONS
! zero-order impurity Green's function, advanced and retarded:
      allocate (gr(N_e,4,4,N_Floquet,N_Floquet), ga(N_e,4,4,N_Floquet,N_Floquet))
! Lesser than self-energy
      allocate (S_l (N_e,4,4,N_Floquet+2,N_Floquet+2))
      allocate (S_r (N_e,4,4,N_Floquet+2,N_Floquet+2))
      allocate (S_a (N_e,4,4,N_Floquet+2,N_Floquet+2))
      allocate (S_M (N_e,4,4,N_Floquet+2,N_Floquet+2))
      allocate (G_r (N_e,4,4,N_Floquet,N_Floquet))
      allocate (G_a (N_e,4,4,N_Floquet,N_Floquet))
      allocate (G_l (N_e,4,4,N_Floquet,N_Floquet))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       print *, ' '
       print *, ' CurrentSC to compute time-dependent SC currents'
       print *, ' '
       print *, ' CONVENTION: bias drop sign::'
       print *, '  sign_L is negative (L is the sample)'
! Hardwired in declarations
       print *, ' '
       print *, ' '
       print *, 'Transmission=', (2*Gamma_0/Bandwidth)/(1._q+Gamma_0/(2.*Bandwidth))**2
       print *, ' '
            if (N_Floquet < 3) then
                print *, ' '
                print *, 'Number of Floquet iterations must be at least 3, but it is=', N_Floquet
                print *, ' STOP '
                print *, ' '
                stop
            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TEST print out BCS DOS
!          omega=Bias_fin
!     call g_d_0(S, J_, U, N_Floquet, Delta_R, Dynes_R,    &
!       &       omega, N_e, theta_i, phi_i, gr, ga, Bandwidth)
!
!          n=(N_Floquet-1)/2+1
!          print *, n
!!
! do i_e=1, N_e
!
!          e_x=-0.5*omega+(i_e-1)*omega/(N_e-1)
!!
!          write (124,'(2g14.4)') e_x*Hartree, -aimag (gr(i_e,1,1,n,n))/pi_d 
! enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! START LOOP on bias
      curr_avant=0.
BIAS: do i_V=1, N_Bias
           V=Bias_ini+(i_V-1)*step_Bias
           omega=V ! un peu con quoi
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set up free left electrode Green's function
     call g_d_0(N_Floquet, Delta_L, Dynes_L,    &
       &       omega, N_e, theta_i, phi_i, gr, ga, Bandwidth)

! Set up Sigma_< (S_l)

     call Sigma_less(omega, Temp, Dynes_L, Delta_L, &
       &    Gamma_0, Bandwidth, phi_, S_l, N_e, N_Floquet)

!           print *, 'out of Sigma_less'

! Set up Sigma advanced and retarded (S_a and S_r)

    call Sigma (omega, Temp, Dynes_L, Delta_L, theta_i, phi_i, &
       &    Gamma_0, Bandwidth, phi_, S_a, S_r, & 
       &    S, U, J_, S_M, N_e, N_Floquet)

!           print *, 'out of Sigma'

! Set up full left electrode Green's functions

     call Green (omega, Temp, S_M, S_l, S_r, S_a, N_e, N_Floquet, gr, ga, G_l, G_a)


!           print *, 'out of Green'


! Compute Current in Floquet space
      
     call curr (omega, N_e, N_Floquet, G_a, S_l, G_l, S_r, curr_0_L, curr_1, curr_2)




! Print out different contributions of the current and conductance

      print *, i_V,'Bias=', V*Hartree, 'Current =', 0.5*curr_0_L/pi_d
      write (2,'(4e14.4)') V*Hartree*1000, 0.5*curr_0_L/pi_d, -0.5*curr_1/pi_d, -0.5*curr_2/pi_d 
! results in atomic units
      if (i_V>1) then
      write (3,'(2e14.4)') V*Hartree*1000, 0.5*(curr_0_L-curr_avant)/step_Bias  !Results in G0
! MODIFIED June 3rd 2022. I was missing a 1/2pi_d from the Floquet integral in the current
!
      endif
      curr_avant=curr_0_L

      call analyses (omega, N_e, N_Floquet, G_a, S_r) ! Plot different contributions to analyse the data

      enddo BIAS

     print *, 'Transmission=', (2*Gamma_0/Bandwidth)/(1._q+Gamma_0/(2.*Bandwidth))**2
     write (2,*) 'Transmission=', (2*Gamma_0/Bandwidth)/(1._q+Gamma_0/(2.*Bandwidth))**2
! END LOOP on bias
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program currentSC
