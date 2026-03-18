Module Tools
Use declarations
CONTAINS
!
      subroutine analyses (omega, N_e, N_Floquet, G_a, S_r) ! Plot different contributions to analyse the data
         integer :: i_e, n, n1, N_e, N_FLoquet
         real (q):: e, omega
         complex (qc), intent (in):: G_a (:,:,:,:,:),  S_r (:,:,:,:,:)

         open (4, file='PDOS0.dat')
         open (8, file='PDOSn.dat')
         open (9, file='Selfn.dat')

         nzero = (N_Floquet+1)/2


         do i_e=1, N_e
         e=-0.5*omega+(i_e-1)*omega/(N_e-1)

             write (4,'(4g14.4)') omega*HARTREE*1000,e*HARTREE*1000, aimag (G_a (i_e,1,1,nzero,nzero)),  &
                   aimag (G_a (i_e,4,4,nzero,nzero))
             write (8,'(8g14.4)') omega*HARTREE*1000,e*HARTREE*1000, aimag (G_a (i_e,1,1,nzero+1,nzero+1)),  & 
                   aimag (G_a (i_e,4,4,nzero+1,nzero+1)), aimag (G_a (i_e,1,1,nzero+2,nzero+2)),  & 
                   aimag (G_a (i_e,4,4,nzero+2,nzero+2)), aimag (G_a (i_e,1,1,nzero-1,nzero-1)),  &
                   aimag (G_a (i_e,4,4,nzero-1,nzero-1))   
             write (9,'(8g14.4)') omega*HARTREE*1000,e*HARTREE*1000, aimag (S_r (i_e,1,1,nzero,nzero)),  & 
                   aimag (S_r (i_e,4,4,nzero,nzero)), aimag (S_r (i_e,1,1,nzero+1,nzero+1)), & 
                   aimag (S_r (i_e,4,4,nzero+1,nzero+1)), aimag (S_r (i_e,1,1,nzero+2,nzero+2)), & 
                   aimag (S_r (i_e,4,4,nzero+2,nzero+2))
          
         enddo

             write (4,'(4g14.4)')
             write (8,'(8g14.4)')
             write (9,'(8g14.4)')

         return

      end subroutine analyses 

end Module Tools
