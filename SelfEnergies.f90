module SelfEnergies
Use declarations

CONTAINS

! Self less
     subroutine Sigma_less(omega, Temp, Dynes_L, Delta_L, &
       &    Gamma_0, Bandwidth, phi_, S_l, N_e, N_Floquet)

     implicit none
     integer :: N_e, N_Floquet, i_e, n, np
     real (q) :: e_x, e, omega, Temp, Dynes_L
     real (q) :: Delta_L,  Gamma_0, Bandwidth, phi_
     complex (qc), intent (out):: S_l (:,:,:,:,:)
     real (q), dimension (N_e) :: Gamma_
      
     S_l=(0._q,0._q)

! Diagonal in Floquet
         do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)
         Gamma_(i_e)=Gamma_0*exp(-(2*e/Bandwidth)**2)

      do n=1, N_Floquet
         np =n-1-(N_Floquet-1)/2

!Achtung
!! sg_L=-1 so that R electrode is at positive bias

         e_x= e+(np-sg_L)*omega
         S_l(i_e,1,1,n,n)= &
          &  ui*fermi(e_x,Temp)*aimag((e_x+ui*Dynes_L)/sqrt(-(e_x+ui*Dynes_L)**2+Delta_L**2))
         S_l(i_e,2,2,n,n)= &
          &  ui*fermi(e_x,Temp)*aimag((e_x+ui*Dynes_L)/sqrt(-(e_x+ui*Dynes_L)**2+Delta_L**2))

         e_x= e+(np+sg_L)*omega
         S_l(i_e,3,3,n,n)= &
          &  ui*fermi(e_x,Temp)*aimag((e_x+ui*Dynes_L)/sqrt(-(e_x+ui*Dynes_L)**2+Delta_L**2))
         S_l(i_e,4,4,n,n)= &
          &  ui*fermi(e_x,Temp)*aimag((e_x+ui*Dynes_L)/sqrt(-(e_x+ui*Dynes_L)**2+Delta_L**2))

      enddo
         enddo

! Non-diagonal

         do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)
      do n=1, N_Floquet+2*sg_L
! because sg_L is negative, not marvellous...
         np =n-1-(N_Floquet-1)/2

         e_x= e+(np-sg_L)*omega
         S_l(i_e,1,3,n,n-2*sg_L)= &
          &  -ui*fermi(e_x,Temp)*Delta_L*exp(ui*phi_)  &
          &  *aimag(1._q/sqrt(-(e_x+ui*Dynes_L)**2+Delta_L**2))
         S_l(i_e,2,4,n,n-2*sg_L)= &
          &  -ui*fermi(e_x,Temp)*Delta_L*exp(ui*phi_)  &
          &  *aimag(1._q/sqrt(-(e_x+ui*Dynes_L)**2+Delta_L**2))

      enddo
         enddo
!
         do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)
      do n=1-2*sg_L, N_Floquet
! because sg_L is negative, not marvellous...

         np =n-1-(N_Floquet-1)/2

         e_x= e+(np+sg_L)*omega
         S_l(i_e,3,1,n,n+2*sg_L)= &
          &  -ui*fermi(e_x,Temp)*Delta_L*exp(-ui*phi_)  &
          &  *aimag(1._q/sqrt(-(e_x+ui*Dynes_L)**2+Delta_L**2))
         S_l(i_e,4,2,n,n+2*sg_L)= &
          &  -ui*fermi(e_x,Temp)*Delta_L*exp(-ui*phi_)  &
          &  *aimag(1._q/sqrt(-(e_x+ui*Dynes_L)**2+Delta_L**2))

      enddo
         enddo


         do i_e=1, N_e
         S_l (i_e,:,:,:,:)=Gamma_(i_e)*S_l (i_e,:,:,:,:)
         enddo
                        return

     end subroutine
!
! Self retarded and advanced
!
     subroutine Sigma (omega, Temp, Dynes_L, Delta_L, theta_i, phi_i, &
       &    Gamma_0, Bandwidth, phi_, S_a, S_r, &
       &    S, U, J, S_M, N_e, N_Floquet)
     implicit none
     integer :: N_e, N_Floquet, i_e, n, np
     real (q) :: S, J, U
     real (q) :: e_x, e, omega, Temp, Dynes_L
     real (q) :: Delta_L, Delta_R, Gamma_0, Bandwidth, phi_, theta_i, phi_i 
     complex (qc), intent (out):: S_a (:,:,:,:,:)
     complex (qc), intent (out):: S_r (:,:,:,:,:)
     complex (qc), intent (out):: S_M (:,:,:,:,:)
     real (q), dimension (N_e) :: Gamma_

       S_r = (0._q,0._q)
       S_a = (0._q,0._q)
       S_M = (0._q,0._q)

! Diagonal
         do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)
         Gamma_(i_e)=Gamma_0*exp(-(e/(4*Bandwidth))**2)
      do n=1, N_Floquet
         np =n-1-(N_Floquet-1)/2

         e_x= e+(np-sg_L)*omega
         S_r (i_e,1,1,n,n) =  &
             (-e_x-ui*Dynes_L)/sqrt(Delta_L**2-(e_x+ui*Dynes_L)**2)
         S_a (i_e,1,1,n,n) =  &
             (-e_x+ui*Dynes_L)/sqrt(Delta_L**2-(e_x-ui*Dynes_L)**2)
         S_r (i_e,2,2,n,n) = &
             (-e_x-ui*Dynes_L)/sqrt(Delta_L**2-(e_x+ui*Dynes_L)**2)
         S_a (i_e,2,2,n,n) =  &
             (-e_x+ui*Dynes_L)/sqrt(Delta_L**2-(e_x-ui*Dynes_L)**2)

         e_x= e+(np+sg_L)*omega
         S_r (i_e,3,3,n,n) = &
             (-e_x-ui*Dynes_L)/sqrt(Delta_L**2-(e_x+ui*Dynes_L)**2)
         S_a (i_e,3,3,n,n) =  &
             (-e_x+ui*Dynes_L)/sqrt(Delta_L**2-(e_x-ui*Dynes_L)**2)
         S_r (i_e,4,4,n,n) = &
             (-e_x-ui*Dynes_L)/sqrt(Delta_L**2-(e_x+ui*Dynes_L)**2)
         S_a (i_e,4,4,n,n) =  &
             (-e_x+ui*Dynes_L)/sqrt(Delta_L**2-(e_x-ui*Dynes_L)**2)

      enddo
         enddo

! Non-diagonal
         do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)
      do n=1, N_Floquet+2*sg_L
         np =n-1-(N_Floquet-1)/2

         e_x= e+(np-sg_L)*omega
         S_r (i_e,1,3,n,n-2*sg_L) = &
             Delta_L*exp(ui*phi_)/sqrt(Delta_L**2-(e_x+ui*Dynes_L)**2)
         S_a (i_e,1,3,n,n-2*sg_L) =  &
             Delta_L*exp(ui*phi_)/sqrt(Delta_L**2-(e_x-ui*Dynes_L)**2)
         S_r (i_e,2,4,n,n-2*sg_L) = &
             Delta_L*exp(ui*phi_)/sqrt(Delta_L**2-(e_x+ui*Dynes_L)**2)
         S_a (i_e,2,4,n,n-2*sg_L) =  &
             Delta_L*exp(ui*phi_)/sqrt(Delta_L**2-(e_x-ui*Dynes_L)**2)

      enddo
         enddo

         do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)
      do n=1-2*sg_L, N_Floquet
         np =n-1-(N_Floquet-1)/2

         e_x= e+(np+sg_L)*omega
         S_r (i_e,3,1,n,n+2*sg_L) = &
             Delta_L*exp(-ui*phi_)/sqrt(Delta_L**2-(e_x+ui*Dynes_L)**2)
         S_a (i_e,3,1,n,n+2*sg_L) =  &
             Delta_L*exp(-ui*phi_)/sqrt(Delta_L**2-(e_x-ui*Dynes_L)**2)
         S_r (i_e,4,2,n,n+2*sg_L) = &
             Delta_L*exp(-ui*phi_)/sqrt(Delta_L**2-(e_x+ui*Dynes_L)**2)
         S_a (i_e,4,2,n,n+2*sg_L) =  &
             Delta_L*exp(-ui*phi_)/sqrt(Delta_L**2-(e_x-ui*Dynes_L)**2)

      enddo
         enddo


         do i_e=1, N_e
         S_r (i_e,:,:,:,:)=0.5*Gamma_ (i_e)*S_r (i_e,:,:,:,:)
         S_a (i_e,:,:,:,:)=0.5*Gamma_(i_e)*S_a (i_e,:,:,:,:)
         enddo

! Magnetic contribution (due to its locality), it
! is the same regardless of _r or _a  _l
! Since it is a potential-like term, there is no _l
         do n= 1, N_Floquet
         S_M (:,1,1,n,n) = 0.5*J*S*cos(theta_i)+U
         S_M (:,1,2,n,n) = 0.5*J*S*sin(theta_i)*exp(-ui*phi_i)
         S_M (:,2,1,n,n) = 0.5*J*S*sin(theta_i)*exp(ui*phi_i)
         S_M (:,2,2,n,n) = -0.5*J*S*cos(theta_i)+U
         S_M (:,3,3,n,n) = 0.5*J*S*cos(theta_i)-U
         S_M (:,3,4,n,n) = 0.5*J*S*sin(theta_i)*exp(ui*phi_i)
         S_M (:,4,3,n,n) = 0.5*J*S*sin(theta_i)*exp(-ui*phi_i)
         S_M (:,4,4,n,n) = -0.5*J*S*cos(theta_i)-U
         enddo


                        return

     end subroutine
!
! Fermi occupation function
!
    function fermi (e, T)
      real (q) :: fermi, e, T
      if ( e > 0._q ) then
          Fermi = exp(-e/T)/(exp(-e/T)+1._q)
       else
         Fermi = 1._q/(exp(e/T)+1._q)
       endif
      return
    end function fermi
end module SelfEnergies
