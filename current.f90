module current
Use declarations
CONTAINS
     subroutine curr(omega, N_e, N_Floquet, G_a, S_l, G_l, S_r, curr_0_L, curr_1, curr_2)
        implicit none
        integer :: i_e, n, n1, N_e, N_FLoquet
        real (q):: curr_0_L, curr_1, curr_2, omega
        complex (qc), intent (in)::G_a(:,:,:,:,:), G_l (:,:,:,:,:)
        complex (qc), intent (in)::S_l (:,:,:,:,:), S_r (:,:,:,:,:)

        
!
! Only compute the first Floquet 
! component (m=0 corresponds to the dc current)

! Trapeze-rule integration on energies

        
! Left electrode

        i_e=1
     curr_1 = 0

     do n=1, N_Floquet
     do n1=1, N_Floquet
     curr_1 = curr_1 + real( &
      &  sum(S_r(i_e,1,:,n,n1)*G_l(i_e,:,1,n1,n)) &
      & +sum(S_r(i_e,2,:,n,n1)*G_l(i_e,:,2,n1,n)) &
      & -sum(S_r(i_e,3,:,n,n1)*G_l(i_e,:,3,n1,n)) &
      & -sum(S_r(i_e,4,:,n,n1)*G_l(i_e,:,4,n1,n)) )
     enddo
     enddo


     do i_e = 2, N_e-1
     do n=1, N_Floquet
     do n1=1, N_Floquet

     curr_1 =curr_1 +2*real( &
      &  sum(S_r(i_e,1,:,n,n1)*G_l(i_e,:,1,n1,n)) &
      & +sum(S_r(i_e,2,:,n,n1)*G_l(i_e,:,2,n1,n)) &
      & -sum(S_r(i_e,3,:,n,n1)*G_l(i_e,:,3,n1,n)) &
      & -sum(S_r(i_e,4,:,n,n1)*G_l(i_e,:,4,n1,n))) 
 
     enddo
     enddo
     enddo

       i_e=N_e

     do n=1, N_Floquet
     do n1=1, N_Floquet
     curr_1 =curr_1 +real( &
      &  sum(S_r(i_e,1,:,n,n1)*G_l(i_e,:,1,n1,n)) &
      & +sum(S_r(i_e,2,:,n,n1)*G_l(i_e,:,2,n1,n)) &
      & -sum(S_r(i_e,3,:,n,n1)*G_l(i_e,:,3,n1,n)) &
      & -sum(S_r(i_e,4,:,n,n1)*G_l(i_e,:,4,n1,n))) 
     enddo
     enddo

        i_e=1
     curr_2 = 0

     do n=1, N_Floquet
     do n1=1, N_Floquet
     curr_2 = curr_2 + real( &
      &  sum(S_l(i_e,1,:,n,n1)*G_a(i_e,:,1,n1,n)) &
      & +sum(S_l(i_e,2,:,n,n1)*G_a(i_e,:,2,n1,n)) &
      & -sum(S_l(i_e,3,:,n,n1)*G_a(i_e,:,3,n1,n)) &
      & -sum(S_l(i_e,4,:,n,n1)*G_a(i_e,:,4,n1,n)))
     enddo
     enddo


     do i_e = 2, N_e-1
     do n=1, N_Floquet
     do n1=1, N_Floquet

     curr_2 =curr_2 +2*real( & 
      &  sum(S_l(i_e,1,:,n,n1)*G_a(i_e,:,1,n1,n)) &
      & +sum(S_l(i_e,2,:,n,n1)*G_a(i_e,:,2,n1,n)) &
      & -sum(S_l(i_e,3,:,n,n1)*G_a(i_e,:,3,n1,n)) &
      & -sum(S_l(i_e,4,:,n,n1)*G_a(i_e,:,4,n1,n)))

     enddo
     enddo
     enddo

       i_e=N_e

     do n=1, N_Floquet
     do n1=1, N_Floquet
     curr_2 =curr_2 +real( &
      &  sum(S_l(i_e,1,:,n,n1)*G_a(i_e,:,1,n1,n)) &
      & +sum(S_l(i_e,2,:,n,n1)*G_a(i_e,:,2,n1,n)) &
      & -sum(S_l(i_e,3,:,n,n1)*G_a(i_e,:,3,n1,n)) &
      & -sum(S_l(i_e,4,:,n,n1)*G_a(i_e,:,4,n1,n)))
     enddo
     enddo

     curr_1= curr_1 * abs(omega)/(N_e-1)
     curr_2= curr_2 * abs(omega)/(N_e-1)

     curr_0_L = 0.5*(curr_1+curr_2)

                return
     end subroutine curr
end module current
