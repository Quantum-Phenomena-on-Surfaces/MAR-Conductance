module GreenFunctions
Use declarations
Use algebra
Use selfenergies

CONTAINS


   subroutine g_d_0 (N_Floquet, Delta, Dynes,    &
       &       omega, N_e, theta_i, phi_i, gr, ga, Bandwidth)

   implicit none
   integer :: i_e, N_e, N_Floquet, n, np
   real (q) :: Bandwidth
   real (q) :: omega, theta_i, phi_i, Delta, Dynes
   complex (qc) :: e, den
   complex (qc), intent (out):: gr(:,:,:,:,:), ga(:,:,:,:,:)


   gr = (0._q,0._q)
   ga = (0._q,0._q)

    do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)+ui*Dynes
      do n=1, N_Floquet
         np =n-1-(N_Floquet-1)/2

        den =-1._q/(Bandwidth*sqrt(Delta**2-(e+np*omega)**2))

      gr(i_e,1,1,n,n)=den*(e+np*omega)
      gr(i_e,1,3,n,n)=den*(-Delta)
      gr(i_e,2,2,n,n)=den*(e+np*omega)
      gr(i_e,2,4,n,n)=den*(-Delta)
      gr(i_e,3,1,n,n)=den*(-Delta)
      gr(i_e,3,3,n,n)=den*(e+np*omega)
      gr(i_e,4,2,n,n)=den*(-Delta)
      gr(i_e,4,4,n,n)=den*(e+np*omega)


      enddo
   enddo


    do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)-ui*Dynes
      do n=1, N_Floquet
         np =n-1-(N_Floquet-1)/2

        den =-1._q/(Bandwidth*sqrt(Delta**2-(e+np*omega)**2))

      ga(i_e,1,1,n,n)=den*(e+np*omega)
      ga(i_e,1,3,n,n)=den*(-Delta)
      ga(i_e,2,2,n,n)=den*(e+np*omega)
      ga(i_e,2,4,n,n)=den*(-Delta)
      ga(i_e,3,1,n,n)=den*(-Delta)
      ga(i_e,3,3,n,n)=den*(e+np*omega)
      ga(i_e,4,2,n,n)=den*(-Delta)
      ga(i_e,4,4,n,n)=den*(e+np*omega)


      enddo
   enddo

     return
   end subroutine

! Full Greens function

   subroutine Green (omega, T, S_M, S_l, S_r, S_a, N_e, N_Floquet, gr, ga, G_l, G_a) 
!  subroutine Green (omega, T, S_M, S_l, S_r, S_a, N_e, N_Floquet, gr, ga, G_l, G_r, G_a) 
     implicit none
     integer :: i_e, N_e, N_Floquet, n, m, i, j, n1, m1, i1, j1
     integer :: np
     real (q) :: e, T, ff, omega
     complex (qc), intent (in):: gr(:,:,:,:,:), ga(:,:,:,:,:), S_l (:,:,:,:,:)
     complex (qc), intent (in):: S_r (:,:,:,:,:), S_a (:,:,:,:,:), S_M (:,:,:,:,:)
     complex (qc), intent (out):: G_a(:,:,:,:,:), G_l (:,:,:,:,:)
!    complex (qc), intent (out):: G_r(:,:,:,:,:), G_a(:,:,:,:,:), G_l (:,:,:,:,:)
!    complex (qc) :: gg(4*N_Floquet,4*N_Floquet)
     complex (qc) :: aa(4*N_Floquet,4*N_Floquet)
     complex (qc) :: ab(4*N_Floquet,4*N_Floquet)
!    complex (qc) :: ssr(4*N_Floquet,4*N_Floquet)
     complex (qc) :: ssm(4*N_Floquet,4*N_Floquet)
     complex (qc) :: ssa(4*N_Floquet,4*N_Floquet)
!    complex (qc) :: GGr(4*N_Floquet,4*N_Floquet)
     complex (qc) :: GGa(4*N_Floquet,4*N_Floquet)
     complex (qc) :: gl (4, 4, N_Floquet, N_Floquet)
     complex (qc) :: Amb (4, 4, N_Floquet, N_Floquet)
     complex (qc) :: Ama (4, 4, N_Floquet, N_Floquet)
     complex (qc) :: Mat_2 (4, 4, N_Floquet, N_Floquet)
     complex (qc) :: Mat_3 (4, 4, N_Floquet, N_Floquet)


! Retarded and Advanced
! Transform to a matrix


     G_l = zero
     G_a = zero
     aa = zero
     ssa = zero
     ssm = zero 
     gl = zero

     do i_e = 1, N_e

! We compute gl from gr and ga

       e=-0.5*omega+(i_e-1)*omega/(N_e-1)

       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet

          np = n-1-(N_Floquet-1)/2
          ff = fermi (e+np*omega, T)

        gl (i,j,n,n) = ff * (ga(i_e,i,j,n,n)-gr(i_e,i,j,n,n))

       enddo
       enddo
       enddo

       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet
        aa (n+(i-1)*N_Floquet,m+(j-1)*N_Floquet)= &
     &           ga(i_e,i,j,n,m)
        ssa(n+(i-1)*N_Floquet,m+(j-1)*N_Floquet)= &
     &           S_a (i_e,i,j,n,m)
        ssm(n+(i-1)*N_Floquet,m+(j-1)*N_Floquet)= &
     &           S_M (i_e,i,j,n,m)
       enddo
       enddo
       enddo
       enddo

! G^r  Retarded
!    call inversion (gg)
!    GGr = gg-ssr-ssm
!    call inversion (GGr)

! G^a Advanced
     call inversion (aa)
     GGa = aa-ssa-ssm 
     call inversion (GGa)


       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet
!         G_r ( i_e, i, j, n, m) = GGr (n+(i-1)*N_Floquet,m+(j-1)*N_Floquet) 
          G_a ( i_e, i, j, n, m) = GGa (n+(i-1)*N_Floquet,m+(j-1)*N_Floquet) 
       enddo
       enddo
       enddo
       enddo
!
! G^<   Less than
!
! We solve it using:
!         G^<=inv[1-gr*Sigma(2)] * (g< + g<*Sigma(1)*G^a + g^r*Sigma^<*G^a)
! Sigma (1) = S_M + S_a
! Sigma (2) = S_M + S_r
! We start calculating the first matrix
! then inverted it
! Then we calculate the sum of the three terms
! Then we multiply the inverse of the first matrix
! times the sum of the three terms
!           Easy peasy


! We create a matrix Amb=inv[1-gr*Sigma(2)]
! we initialize Amb to the identity matrix in both Nambu+Floquet space
       Amb = zero
     
       do i=1, 4
       do n=1, N_Floquet
          Amb ( i, i, n, n) =  ur
       enddo
       enddo
! product of g*Sigma, matrix form
       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet

          do j1=1, 4
          do m1=1, N_Floquet
          Amb ( i, j, n, m) = Amb ( i, j, n, m) - gr (i_e, i, j1, n, m1)*  &
     &        (S_r (i_e, j1, j, m1, m)+S_M (i_e, j1, j, m1, m))
          enddo
          enddo

       enddo
       enddo
       enddo
       enddo
 
       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet
          ab (n+(i-1)*N_Floquet,m+(j-1)*N_Floquet) = Amb(i,j,n,m)
       enddo
       enddo
       enddo
       enddo

       call inversion (ab) 

! Finally 
! We create a matrix Amb=inv[1-gr*Sigma(2)]
       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet
          Amb (i,j,n,m) = ab(n+(i-1)*N_Floquet,m+(j-1)*N_Floquet)
       enddo
       enddo
       enddo
       enddo

!Sum of the three matrices
!1st matrix is g^<
!2nd matrix is g^< * Sigma(1) * G^a

          Ama = zero

       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet

          do j1=1, 4
          do m1=1, N_Floquet
          Ama ( i, j, n, m) = Ama ( i, j, n, m) + gl (i, j1, n, m1)*  &
     &        (S_a (i_e, j1, j, m1, m)+S_M (i_e, j1, j, m1, m))
          enddo
          enddo

       enddo
       enddo
       enddo
       enddo

          Mat_2 = zero

       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet

          do j1=1, 4
          do m1=1, N_Floquet
          Mat_2 ( i, j, n, m) = Mat_2 ( i, j, n, m) + Ama (i, j1, n, m1)*  &
     &        G_a (i_e, j1, j, m1, m)
          enddo
          enddo

       enddo
       enddo
       enddo
       enddo

!3rd matrix is g^r * S^< * G^a
          Ama = zero
       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet

          do j1=1, 4
          do m1=1, N_Floquet
          Ama ( i, j, n, m) = Ama ( i, j, n, m) + gr (i_e, i, j1, n, m1)*  &
     &        S_l (i_e, j1, j, m1, m)
          enddo
          enddo

       enddo
       enddo
       enddo
       enddo

          Mat_3 = zero
       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet

          do j1=1, 4
          do m1=1, N_Floquet
          Mat_3 ( i, j, n, m) = Mat_3 ( i, j, n, m) + Ama (i, j1, n, m1)*  &
     &        G_a (i_e, j1, j, m1, m)
          enddo
          enddo

       enddo
       enddo
       enddo
       enddo

! G^<
       do i=1, 4
       do j=1, 4
       do n=1, N_Floquet
       do m=1, N_Floquet

          do j1=1, 4
          do m1=1, N_Floquet
          G_l (i_e, i, j, n, m) = G_l (i_e, i, j, n, m) + Amb (i, j1, n, m1)*  &
     &        (gl(j1, j, m1, m) + Mat_2 (j1, j, m1, m) + Mat_3 (j1, j, m1, m))
          enddo
          enddo

       enddo
       enddo
       enddo
       enddo

      

    enddo

   end subroutine
end module GreenFunctions
