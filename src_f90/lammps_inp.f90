! To generate LAMMPS input file for SIC (VEC-LiSTFSI) systems
! Use in conjunction with lmp_params.f90 and ran_numbers.f90
PROGRAM LAMMPSINP

  USE PARAMS

  IMPLICIT NONE
  
  LOGICAL :: input_coor = .false.
  REAL :: bondl, bondlsq
  INTEGER :: ierror,narg

  bondlsq = 0

  CALL SYSTEM_CLOCK(S)

  OPEN (unit = outfile, file = "lmp_input.dat", status ="replace",&
       & action="write",iostat=ierror)
  
  IF(ierror /= 0) STOP "Cannot open lmp_input.txt"
  
  CALL SANITY_CHECKS()
  CALL CREATEFILE() 
  CALL COMPUTE_BOX()
  CALL ALLOCATE_ARRAYS()
  CALL INPCOR()
  CALL SET_IMGFLAGS()
  CALL LMP_COORD()
  CALL DEALLOCATE_ARRAYS()
  CALL CROSS_CHECK_NUM_BOND_TYPES()
  CALL CLEAN_AND_CLOSE_ALL()

END PROGRAM LAMMPSINP

!--------------------------------------------------------------------

SUBROUTINE LMP_COORD()

  USE PARAMS
  
  IMPLICIT NONE
  
  INTEGER :: i,j, k, ierror
  INTEGER ::  bondid, anglid, dihdid
  REAL ::  massval
  i = 1
  
  PRINT *, "Writing LAMMPS Datafile .. "

20 FORMAT(5X,I0,2X,A)
22 FORMAT(5X,I0,2X,A)
24 FORMAT(5X,I0,2X,F14.6,2X,A)
  
  OPEN (unit=10, file = datafile, status="replace",action=&
       &"write",iostat = ierror)
  
  IF(ierror /= 0) STOP "Failed to open datafile"
  
  IF(frac_unpoly == 0) THEN 
     WRITE (10,*) "Data for CG p(mVEC-r-nLiMTSFI) simulations w/o unpo&
          &lymerized VEC monomers"
  ELSE
     WRITE (10,*) "Data for CG p(mVEC-r-nLiMTSFI) simulations with unp&
          &olymerized VEC monomers"
  END IF

  WRITE (10,*) 
  WRITE (10,20) totblobs, "atoms"

  IF(numbondtypes /= 0) THEN
     WRITE (10,20) nbonds, "bonds"
  ELSE
     WRITE (10,20) 0, "bonds"
  END IF

  IF(numangltypes /= 0) THEN
     WRITE (10,20) nangls, "angles"
  ELSE
     WRITE (10,20) 0, "angles"
  END IF

  IF(numdihdtypes /= 0) THEN
     WRITE (10,20) ndihds, "dihedrals"
  ELSE
     WRITE (10,20) 0, "dihedrals"
  END IF

  WRITE (10,20) 0, "impropers"

 
  WRITE (10,20) numatomtypes, "atom types"
  WRITE (10,20) SUM(bflag_arr(:,1)), "bond types"
  WRITE (10,22) numangltypes, "angle types"
  WRITE (10,22) numdihdtypes, "dihedral types"
  WRITE (10,22) 0, "improper types"

  WRITE (10,*)
  WRITE (10,24) 0, boxl_x, "xlo xhi"
  WRITE (10,24) 0, boxl_y, "ylo yhi"
  WRITE (10,24) 0, boxl_z, "zlo zhi"
  WRITE (10,*)
  WRITE (10,*) "Masses"
  WRITE (10,*)
  
  ! Writing Masses
  
  WRITE(10,'(I0,1X,F14.8)') 1, 12.4 ! VEC without C=O
  WRITE(10,'(I0,1X,F14.8)') 2, 4.0  ! C=O of VEC
  WRITE(10,'(I0,1X,F14.8)') 3, 48.7 ! MTFSI

  IF(frac_unpoly .NE. 0.0) THEN
     WRITE(10,'(I0,1X,F14.8)') 4, 12.4 ! unpolymerized VEC w/o C=O 
     WRITE(10,'(I0,1X,F14.8)') 5, 4.0  ! C=O of unpolymerized VEC
  END IF

  WRITE(10,'(I0,1X,F14.8)') 6, 1.0  ! Li

  ! Writing atomic coordinates
  
  WRITE (10,*) 
  WRITE (10,*) "Atoms"
  WRITE (10,*)

  DO i = 1,totblobs

     IF(unwrapped .EQV. .false.) THEN !Wrap it and write ix,iy,iz
        WRITE(10,'(3(I0,1X),4(F14.6,1X),3(I0,1X))') aidvals(i,1),&
             & aidvals(i,2), aidvals(i,3), charge(i), rxyz(i,1),&
             & rxyz(i,2),rxyz(i,3), ixyz(i,1), ixyz(i,2), ixyz(i,3)

     ELSE ! Ideally ixyz should be zero
        WRITE(10,'(3(I0,1X),4(F14.6,1X))') aidvals(i,1),aidvals(i,2),&
             & aidvals(i,3), charge(i), rxyz(i,1) + boxl_x*ixyz(i,1),&
             & rxyz(i,2) + boxl_y*ixyz(i,2),rxyz(i,3)+ boxl_z*ixyz(i&
             &,3)
        
     END IF
    
  END DO

  IF(numbondtypes /= 0) THEN

     ! Writing Bond Details  
     
     WRITE (10,*)
     WRITE (10,*) "Bonds"
     WRITE (10,*)
     
     DO i = 1,nbonds
        
        WRITE(10,'(4(I0,2X))') topo_bond(i,1), topo_bond(i,2),&
             & topo_bond(i,3), topo_bond(i,4)

     END DO
        
  END IF

  CLOSE(unit = 10)

END SUBROUTINE LMP_COORD

!----------------------------------------------------------------------

SUBROUTINE SET_IMGFLAGS()

  USE PARAMS

  INTEGER :: i,j,k
  REAL :: rx, ry, rz

  ixyz = 0 ! If I don't assign, errors may creep

  ! Polymerized molecules
  DO i = 1,N_poly
     
     ! First blob
     k = (i-1)*blob_per_ch + 1

     ixyz(k,1) = 0
     ixyz(k,2) = 0
     ixyz(k,3) = 0
     
     DO j = 1,blob_per_ch-1
        
        k = (i-1)*blob_per_ch + j
        
        ! Branched VEC blobs
        IF (aidvals(k+1,3) > 1 .OR. aidvals(k+1,3) .LE. CG_per_mon)&
             & THEN

           rx = rxyz(k,1) - rxyz(k+1,1)
           ry = rxyz(k,2) - rxyz(k+1,2)
           rz = rxyz(k,3) - rxyz(k+1,3)
           CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
           CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
           CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))

        ELSE
           ! VEC blob in main chain or the anion blob
           IF (aidvals(k,3) == 1 .OR. aidvals(k,3) == CG_per_mon+1)&
                & THEN

              rx = rxyz(k,1) - rxyz(k+1,1)
              ry = rxyz(k,2) - rxyz(k+1,2)
              rz = rxyz(k,3) - rxyz(k+1,3)
              CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
              CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
              CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))

            
           ELSEIF (aidvals(k,3) == CG_per_mon) THEN
              ! Image flags for blobs in main chain right after the
              ! branch ends
              rx = rxyz(k-CG_per_mon+1,1) - rxyz(k+1,1)
              ry = rxyz(k-CG_per_mon+1,2) - rxyz(k+1,2)
              rz = rxyz(k-CG_per_mon+1,3) - rxyz(k+1,3)
              CALL IMGFLAGS(rx,ixyz(k-CG_per_mon+1,1),boxl_x,ixyz(k+1,1))
              CALL IMGFLAGS(ry,ixyz(k-CG_per_mon+1,2),boxl_y,ixyz(k+1,2))
              CALL IMGFLAGS(rz,ixyz(k-CG_per_mon+1,3),boxl_z,ixyz(k+1,3))
                           
           END IF

        END IF
             
     END DO
     
  END DO

  ! Unpolymerized VEC molecules
  DO i = 1, T_unpoly_VEC
     
     ! First blob of the VEC unpolymerized molecule
     k = N_poly*blob_per_ch + (i-1)*CG_per_mon + 1

     ixyz(k,1) = 0
     ixyz(k,2) = 0
     ixyz(k,3) = 0

     DO j = 1, CG_per_mon-1

        k = N_poly*blob_per_ch + (i-1)*CG_per_mon + j
        
        rx = rxyz(k,1) - rxyz(k+1,1)
        ry = rxyz(k,2) - rxyz(k+1,2)
        rz = rxyz(k,3) - rxyz(k+1,3)
        CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
        CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
        CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))

     END DO

  END DO

END SUBROUTINE SET_IMGFLAGS

!----------------------------------------------------------------------

SUBROUTINE IMGFLAGS(dist, img, boxl, imgout)
    
  USE PARAMS
  IMPLICIT NONE
  
  REAL, INTENT(IN) :: dist,boxl
  INTEGER, INTENT(IN) :: img
  INTEGER, INTENT(OUT) :: imgout
  INTEGER :: nx
  
  IF(dist > boxl/2) THEN
     
     nx = img + 1
     
  ELSEIF(dist < -boxl/2) THEN
     
     nx = img - 1
     
  ELSE
     
     nx = img
     
  END IF
  
  imgout = nx
  
END SUBROUTINE IMGFLAGS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_BOX()

  USE PARAMS

  IMPLICIT NONE
 
  WRITE(outfile,*) "Simulation inputs ...."

  volbox = REAL(totblobs)/REAL(density)
  boxl_x = volbox**(1.0/3.0)
  boxl_y = volbox**(1.0/3.0)
  boxl_z = volbox**(1.0/3.0)

  WRITE(outfile,*) "----------System level details-------------------"
  WRITE(outfile,*) "Total # of particles/blobs: ", totblobs
  WRITE(outfile,*) "# of CG blobs per polymer (VEC) monomer: ",&
       & CG_per_mon
  WRITE(outfile,*) "Total # of polymer chains: ", N_poly
  WRITE(outfile,*) "Total # of molecules (chains+unpolymerized mons): &
       &", nmols
  WRITE(outfile,*) "Total # of anion blobs: ", N_anions
  WRITE(outfile,*) "Total # of cations: ", N_cations
  WRITE(outfile,*) "Fraction of unpolymerized-mers: ", frac_unpoly
  WRITE(outfile,*) "Total # of VEC monomers: ", T_VEC_mons
  WRITE(outfile,*) "Total # of polymerized VEC mons: ", T_poly_VEC
  WRITE(outfile,*) "Total # of unpolymerized-mers: ", T_unpoly_VEC
  
  WRITE(outfile,*) "----------Chain level details--------------------"
  WRITE(outfile,*) "Ratio between anions and VEC mons per chain: ",&
       & an_poly_rat
  WRITE(outfile,*) "# of VEC monomers per chain: ", VEC_per_ch
  WRITE(outfile,*) "# of blobs per chain: ", blob_per_ch
  WRITE(outfile,*) "# of anion blobs per chain: ", ideal_an_per_ch

  WRITE(outfile,*) "----------Box details----------------------------"
  WRITE(outfile,*) "LX/LY/LZ: ", boxl_x, boxl_y, boxl_z
  WRITE(outfile,*) "Box volume: ", volbox
  WRITE(outfile,*) "Density: ", density

  WRITE(outfile,*) "----------Topology details-----------------------"
  WRITE(outfile,*) "Number of atomtypes: ", numatomtypes
  WRITE(outfile,*) "Number of bondtypes: ", numbondtypes
  WRITE(outfile,*) "Number of angletypes: ", numangltypes
  
END SUBROUTINE COMPUTE_BOX

!--------------------------------------------------------------------

SUBROUTINE INPCOR()
  
  USE PARAMS

  IMPLICIT NONE
  
  INTEGER :: i,j,k,u,v,ierror
  INTEGER :: an_per_ch, cgcnt, poly_mon
  INTEGER :: bid_start, bondid_new, bondtemp
  REAL :: theta, phi
  REAL :: csum
  LOGICAL :: arr_bound
  
  CALL RAN_INIT(S,X)
  
  WRITE(outfile,*) "Random Initial Configuration : NRRW"

  ! Create polymeric chains

  WRITE(outfile,*) "Generating chain configurations .. "

  PRINT *, "------Generating mixed VEC-anion chains------------"
  i = 1; bondid_new = 0; bflag_arr = 0
  
  ! Polymerized chains
  DO WHILE (i .LE. N_poly)

     print *, "Chain ID: ", i
     bondtemp = bondid_new; bid_start = bondid_new + 1
     an_per_ch = 0
     cgcnt = 1
     j = 1
     poly_mon = 1

     ! First blob is VEC
     k      = (i-1)*blob_per_ch + 1     
     theta  = math_pi*RAN1(X)
     phi    = 2*math_pi*RAN1(X)

     CALL CREATE_FIRST_VEC_MONOMER(i,theta,phi,cgcnt,k,bondtemp&
          &,bid_start)

     j = 2
     DO WHILE (j .LE. M_poly) !j runs over MONOMERS and not BLOBS

        IF(is_ion_sep == 0) THEN ! If anions are part of the backbone

           ! Second monomer onwards can be STSFI (anion)
           IF (j .NE. M_poly .AND. RAN1(X) .LE.  an_poly_rat) THEN ! can
              ! be attached to blob-1 only

              theta     = math_pi*RAN1(X)
              phi       = 2*math_pi*RAN1(X)
              
              bondtemp  = bondtemp + 1

              ! Check if more blobs/bonds are created than the array
              ! bounds
              CALL CHECK_ARRAY_BOUNDS(i,k+1,bondtemp,arr_bound)

              IF (arr_bound .EQV. .TRUE.) THEN
                 CALL CREATE_ANION_IN_MIXEDCHAIN(i,cgcnt,theta,phi,k&
                      &,bondtemp,an_per_ch)
              ELSE
                 j = M_poly+ideal_an_per_ch
              END IF
                 
           ELSE
              
              ! CG Blobs
              cgcnt = 1
              DO WHILE (cgcnt .LE. CG_per_mon)
                           
                 k = (i-1)*blob_per_ch + CG_per_mon*poly_mon + cgcnt +&
                      & an_per_ch
                 bondtemp     = bondtemp + 1
                 ! Check if more blobs/bonds are created than the array
                 ! bounds
                 CALL CHECK_ARRAY_BOUNDS(i,k,bondtemp,arr_bound)
                 
                 IF (arr_bound .EQV. .TRUE.) THEN
           
                    theta     = math_pi*RAN1(X)
                    phi       = 2*math_pi*RAN1(X)
                    CALL CREATE_VEC_IN_MIXEDCHAIN(i,cgcnt,theta,phi&
                         &,k,bondtemp,an_per_ch)
                    cgcnt = cgcnt + 1

                 ELSE

                    j = M_poly
                    cgcnt = CG_per_mon+1

                 END IF
                 
              END DO

              poly_mon = poly_mon + 1
              
           END IF

           j = j + 1

        ELSE ! if anions are separate from the backbone

           ! CG Blobs
           cgcnt = 1
           DO WHILE (cgcnt .LE. CG_per_mon)
              
              k = (i-1)*blob_per_ch + CG_per_mon*poly_mon + cgcnt
              bondtemp     = bondtemp + 1
              ! Check if more blobs/bonds are created than the array
              ! bounds
              CALL CHECK_ARRAY_BOUNDS(i,k,bondtemp,arr_bound)
              
              IF (arr_bound .EQV. .TRUE.) THEN
                 
                 theta     = math_pi*RAN1(X)
                 phi       = 2*math_pi*RAN1(X)
                 CALL CREATE_PURE_VEC_CHAINS(i,cgcnt,theta,phi,k,&
                      &,bondtemp)
                 cgcnt = cgcnt + 1
                 
              ELSE
                 
                 j = M_poly
                 cgcnt = CG_per_mon+1
                 
              END IF
             
           END DO

           poly_mon = poly_mon + 1
           j = j + 1

        END IF
        
     END DO

     IF (an_per_ch == ideal_an_per_ch .AND. is_ion_sep == 0) THEN

        bondid_new = bondtemp
        CALL CHECK_BOND_TYPES(i, bid_start, bondid_new)
        i = i + 1
        
     END IF
     
  END DO
  PRINT *, "------Generated mixed VEC-anion chains------------"
  

  PRINT *, "---------Polymerized chain data-------------------"
  PRINT *, "Ideal number of polymer blobs: ", N_poly*blob_per_ch
  PRINT *, "Total number of polymer blobs: ", k
  PRINT *, "NCations, NAnions ", N_cations, N_anions
  PRINT *, "Number of bond types in poly VEC: ", btype_poly_VEC
  PRINT *, "------Generated polymerized chains----------------"


  btype_poly_VEC = SUM(bflag_arr(:,1)) ! # of btypes in poly_VEC

  ! Create unpolymerized monomers
  IF (frac_unpoly > 0.0) THEN
     PRINT *, "------Generating unpolymerized chains-------------"
     CALL CREATE_UNPOLYMERIZED_VEC_MONOMERS(i,cgcnt,k,bondtemp,bid_start)
     PRINT *, "------Generated unpolymerized chains--------------"
  END IF
  
  ! Create pure polyanions if necessary
  IF(is_ion_sep) THEN
     PRINT *, "------Generating pure polyanions---------------"
     CALL CREATE_PURE_POLY_ANIONS(i,cgcnt,k,bondtemp)
     i = i + N_poly - 1
     PRINT *, "------Generated pure polyanions----------------"
  END IF


  PRINT *, "------Generating lithium cations-------------------"     
  CALL CREATE_LITHIUM_CATIONS(k)
  PRINT *, "------Generated lithium cations--------------------"     

  IF (unwrapped .EQV. .false.) THEN
     ! PBC
     
     DO i = 1,totblobs
        
        rxyz(i,1) = rxyz(i,1) - boxl_x*floor(rxyz(i,1)/boxl_x)
        rxyz(i,2) = rxyz(i,2) - boxl_y*floor(rxyz(i,2)/boxl_y)
        rxyz(i,3) = rxyz(i,3) - boxl_z*floor(rxyz(i,3)/boxl_z)
        
     END DO

  END IF

  ! Check charge neutrality
  csum = 0.0
  DO i = 1, totblobs

     csum = csum + charge(i)
     
  END DO

  IF(abs(csum) .GT. 0.00000001) THEN

     OPEN(unit = 93,file="charge.txt",action="write",status="replace")

     DO i = 1,totblobs
        WRITE(93,*) i,aidvals(i,3), charge(i)
     END DO

     WRITE(outfile,*) "System not charge neutral", csum
     PRINT *, "ERROR: System not charge neutral", csum
     CLOSE(93)
     STOP
     
  ELSE

     WRITE(outfile,*) "Good Charge Neutrality ", csum
     PRINT *, "Good Charge Neutrality ", csum

  END IF
 
END SUBROUTINE INPCOR

!--------------------------------------------------------------------

SUBROUTINE CHECK_ARRAY_BOUNDS(chid,aid,bid,arr_bound)

  USE PARAMS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: chid,aid, bid
  LOGICAL, INTENT(OUT) :: arr_bound

  arr_bound = .TRUE. ! Within array bounds

  IF (aid > chid*blob_per_ch .OR. bid > chid*(blob_per_ch-1)) THEN

     arr_bound = .FALSE.

  END IF

END SUBROUTINE CHECK_ARRAY_BOUNDS
!--------------------------------------------------------------------

SUBROUTINE CREATE_FIRST_VEC_MONOMER(i,theta,phi,cgcnt,k,bondtemp)

  USE PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i
  REAL, INTENT(INOUT) :: theta, phi
  INTEGER, INTENT(INOUT) :: k, bondtemp, cgcnt
  LOGICAL :: arr_bound

  rxyz(k,1) = RAN1(X)*boxl_x
  rxyz(k,2) = RAN1(X)*boxl_y
  rxyz(k,3) = RAN1(X)*boxl_z
  
  ! Create atom types here
  aidvals(k,1) = k ! Atom ID
  aidvals(k,2) = i ! Mol ID
  aidvals(k,3) = 1 ! Atom type
  charge(k)    = (0.5*(1+CG_per_mon)-REAL(cgcnt))*2*charge_poly
     
  DO cgcnt = 2, CG_per_mon

     k = (i-1)*blob_per_ch + cgcnt
     bondtemp     = bondtemp + 1
     CALL CHECK_ARRAY_BOUNDS(i,k,bondtemp,arr_bound)
     
     theta     = math_pi*RAN1(X)
     phi       = 2*math_pi*RAN1(X)
     rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
     rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
     rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
     aidvals(k,1) = k
     aidvals(k,2) = i
     aidvals(k,3) = cgcnt
     charge(k)    = (0.5*(1+CG_per_mon)-REAL(cgcnt))*2*charge_poly 
     CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1),aidvals(k,1)&
          &,160)
        
  END DO
  
END SUBROUTINE CREATE_FIRST_VEC_MONOMER

!--------------------------------------------------------------------

SUBROUTINE CREATE_ANION_IN_MIXEDCHAIN(i,cgcnt,theta,phi,k,bondtemp&
     &,an_per_ch)

  USE PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i, cgcnt
  REAL, INTENT(IN) :: theta, phi
  INTEGER, INTENT(INOUT) :: k,bondtemp,an_per_ch
  
  ! Generate anions bonded to the backbone

  aidvals(k+1,1) = k + 1
  aidvals(k+1,2) = i
  aidvals(k+1,3) = CG_per_mon + 1
  charge(k+1)    = -1.0
  
  IF(aidvals(k,3) == CG_per_mon + 1) THEN !anion-anion
     rxyz(k+1,1) = rxyz(k,1) + r0init*sin(theta)*cos(phi)
     rxyz(k+1,2) = rxyz(k,2) + r0init*sin(theta)*sin(phi)
     rxyz(k+1,3) = rxyz(k,3) + r0init*cos(theta)
     CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k+1,1),aidvals(k,1),100)
     
  ELSEIF(aidvals(k,3) == 1) THEN ! Wrong choice
     
     PRINT *, "WRONG CHOICE" , k
     STOP
     
  ELSE ! poly_mon - anion bond
     rxyz(k+1,1) = rxyz(k+1-CG_per_mon,1) + r0init*sin(theta)*cos(phi)
     rxyz(k+1,2) = rxyz(k+1-CG_per_mon,2) + r0init*sin(theta)*sin(phi)
     rxyz(k+1,3) = rxyz(k+1-CG_per_mon,3) + r0init*cos(theta)
     CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k+1,1),aidvals(k+1&
          &-CG_per_mon,1),100)
     
  END IF
     
  an_per_ch = an_per_ch + 1 ! update an_per_ch
  k = k + 1 ! pseudo update if anion-anion bond is formed
      
END SUBROUTINE CREATE_ANION_IN_MIXEDCHAIN

!--------------------------------------------------------------------

SUBROUTINE CREATE_VEC_IN_MIXEDCHAIN(i,cgcnt,theta,phi,k,bondtemp)

  USE PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i, cgcnt
  REAL, INTENT(IN) :: theta, phi
  INTEGER, INTENT(IN) :: k,bondtemp

  !CG Blobs
  aidvals(k,1) = k
  aidvals(k,2) = i
  aidvals(k,3) = cgcnt
  charge(k)    = (0.5*(1+CG_per_mon)-REAL(cgcnt))*2*charge_poly
  
  IF (cgcnt == 1) THEN
     IF(aidvals(k-1,3) == CG_per_mon + 1) THEN !poly_mon-anion
        rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
        rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
        rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
        CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1),aidvals(k,1)&
             &,200)
     ELSE !poly_mon_(k-1) - poly_mon_k
        rxyz(k,1) = rxyz(k-CG_per_mon,1) + r0init*sin(theta)*cos(phi)
        rxyz(k,2) = rxyz(k-CG_per_mon,2) + r0init*sin(theta)*sin(phi)
        rxyz(k,3) = rxyz(k-CG_per_mon,3) + r0init*cos(theta)
        CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-CG_per_mon,1)&
             &,aidvals(k,1),200)
     END IF
  ELSE !poly_mon_k - poly_mon_k (blob-blob connection)
     rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
     rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
     rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
     CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1),aidvals(k,1),200)
  END IF

END SUBROUTINE CREATE_VEC_IN_MIXEDCHAIN
  
!--------------------------------------------------------------------

SUBROUTINE CREATE_PURE_VEC_CHAINS(i,cgcnt,theta,phi,k,bondtemp)

  USE PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i, cgcnt
  REAL, INTENT(IN) :: theta, phi
  INTEGER, INTENT(IN) :: k, bondtemp

  !CG Blobs
  aidvals(k,1) = k
  aidvals(k,2) = i
  aidvals(k,3) = cgcnt
  charge(k)    = (0.5*(1+CG_per_mon)-REAL(cgcnt))*2*charge_poly
  
  IF (cgcnt == 1) THEN
     rxyz(k,1) = rxyz(k-CG_per_mon,1) + r0init*sin(theta)*cos(phi)
     rxyz(k,2) = rxyz(k-CG_per_mon,2) + r0init*sin(theta)*sin(phi)
     rxyz(k,3) = rxyz(k-CG_per_mon,3) + r0init*cos(theta)
     CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-CG_per_mon,1),aidvals(k&
          &,1),200)
  ELSE !poly_mon_k - poly_mon_k (blob-blob connection)
     rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
     rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
     rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
     CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1),aidvals(k,1),200)
  END IF

END SUBROUTINE CREATE_PURE_VEC_CHAINS
  
!--------------------------------------------------------------------

SUBROUTINE CREATE_PURE_POLY_ANIONS(i,cgcnt,k,bondtemp)

  USE PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i
  INTEGER, INTENT(INOUT) :: cgcnt,k,bondtemp
  INTEGER :: polyan_id, dum_id
  REAL :: theta, phi
  
  !Polyanion blobs
  DO polyan_id = i, i + N_poly-1

     k = k + 1 
     aidvals(k,1) = k
     aidvals(k,2) = polyan_id
     aidvals(k,3) = CG_per_mon + 1
     charge(k)    = -1

     theta     = math_pi*RAN1(X)
     phi       = 2*math_pi*RAN1(X)
     rxyz(k,1) = RAN1(X)*boxl_x
     rxyz(k,2) = RAN1(X)*boxl_y
     rxyz(k,3) = RAN1(X)*boxl_z

     DO dum_id = 2,ideal_an_per_ch
        k = k + 1
        aidvals(k,1) = k
        aidvals(k,2) = polyan_id
        aidvals(k,3) = cgcnt
        charge(k)    = -1
        
        rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
        rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
        rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
        CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1),aidvals(k,1)&
             &,200)

     END DO

  END DO

  
END SUBROUTINE CREATE_PURE_POLY_ANIONS
  
!--------------------------------------------------------------------

SUBROUTINE CREATE_UNPOLYMERIZED_VEC_MONOMERS(i,cgcnt,k,bondtemp&
     &,bid_start)

  USE PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i
  INTEGER, INTENT(INOUT) :: cgcnt,k,bondtemp,bid_start
  REAL :: theta, phi
  INTEGER :: un_chid

  ! Create unpolymerized polymer (VEC) monomers
  DO un_chid = i, i+T_unpoly_VEC-1

     ! Head VEC blob
     k = k + 1; cgcnt = 1
     theta     = math_pi*RAN1(X)
     phi       = 2*math_pi*RAN1(X)
     rxyz(k,1) = RAN1(X)*boxl_x
     rxyz(k,2) = RAN1(X)*boxl_y
     rxyz(k,3) = RAN1(X)*boxl_z

     ! Create atom types here
     aidvals(k,1) = k ! Atom ID
     aidvals(k,2) = un_chid ! Mol ID
     aidvals(k,3) = CG_per_mon + 2 ! Atom type
     charge(k)    = (0.5*(1+CG_per_mon)-REAL(cgcnt))*2*charge_poly   
     
     IF(un_chid == i) bid_start = bondtemp + 1
     DO cgcnt = 2, CG_per_mon

        k = k + 1
        bondtemp  = bondtemp + 1
        
        theta     = math_pi*RAN1(X)
        phi       = 2*math_pi*RAN1(X)
        rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
        rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
        rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
        aidvals(k,1) = k
        aidvals(k,2) = un_chid
        aidvals(k,3) = CG_per_mon + cgcnt + 1 !cgcnt is from 2; so +1
        charge(k)    = (0.5*(1+CG_per_mon)-REAL(cgcnt))*2*charge_poly 
        CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1),aidvals(k,1)&
             &,260)
        
     END DO

     ! For consistency in bflag_arr
     IF(un_chid == i) CALL CHECK_BOND_TYPES(un_chid,bid_start&
          &,bondtemp)

  END DO

END SUBROUTINE CREATE_UNPOLYMERIZED_VEC_MONOMERS

!--------------------------------------------------------------------

! Create lithium cations
SUBROUTINE CREATE_LITHIUM_CATIONS(k)

  USE PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: k
  INTEGER :: li_id
  
  DO li_id = k+1, k + N_cations
  
     rxyz(li_id,1) = RAN1(X)*boxl_x
     rxyz(li_id,2) = RAN1(X)*boxl_y
     rxyz(li_id,3) = RAN1(X)*boxl_z
     aidvals(li_id,1) = k+1
     aidvals(li_id,2) = N_poly + T_unpoly_VEC + 1
     aidvals(li_id,3) = CG_per_mon*(1 + CEILING(frac_unpoly)) + 2
     charge(li_id)    = 1.0
        
  END DO

END SUBROUTINE CREATE_LITHIUM_CATIONS

!--------------------------------------------------------------------

SUBROUTINE ASSIGN_BOND_TOPO(bid,aid1,aid2,iderr)

  USE PARAMS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: bid, aid1, aid2,iderr
  INTEGER :: btype
  INTEGER :: atype1, atype2

  atype1 = aidvals(aid1,3)
  atype2 = aidvals(aid2,3)
  
  IF (aidvals(aid1,2) .LE. N_poly) THEN !Polymerized molecules

     IF (atype1 == atype2) THEN

        IF (atype1 == 1) THEN
           btype = 1 ! VEC_head-VEC_head
        ELSEIF (atype1 == CG_per_mon + 1) THEN
           btype = CG_per_mon + 2 ! Anion-Anion
        ELSE
           PRINT *, "Bond error ", bid, aid1, aid2,atype1,atype2,iderr
        END IF

     ELSEIF ((atype1 == 1 .OR. atype1 == CG_per_mon+1) .AND. (atype2 ==&
          & 1 .OR. atype2 == CG_per_mon+1) ) THEN

        btype = CG_per_mon + 1 ! VEC-Anion

     ELSEIF((atype1 .NE. CG_per_mon+1) .AND. (atype2 .NE. CG_per_mon+1))&
          & THEN

        btype = MAX(atype1, atype2) ! Branched mon - branched mon

     ELSE

        PRINT *, "Bond error ", bid, aid1, aid2,atype1,atype2,iderr

     END IF

  ELSE !Unpolymerized VEC molecules

     ! Accounts for anion-anion bonding not formed
     btype = min(atype1,atype2) + btype_poly_VEC - (CG_per_mon + 1)

  END IF

  topo_bond(bid,1) = bid
  topo_bond(bid,2) = btype
  topo_bond(bid,3) = aid1
  topo_bond(bid,4) = aid2

END SUBROUTINE ASSIGN_BOND_TOPO
!--------------------------------------------------------------------

SUBROUTINE CREATEFILE() !CREATEFILE(narg)

  USE PARAMS
  IMPLICIT NONE

  CHARACTER (LEN = 6) :: nmon_char
  CHARACTER (LEN = 7) :: prefix = "VECdata"
  CHARACTER (LEN = 4) :: frat_char
 
  WRITE(nmon_char,"(I0)") M_poly
  WRITE(frat_char,"(F4.2)") an_poly_rat

  datafile = prefix//'_'//trim(adjustl(nmon_char))//'_'&
       &//trim(adjustl(frat_char))//'.dat'
     
  WRITE(outfile,*) "Data file generated for simulation is ",&
       & trim(datafile)

  PRINT *, "Data file generated for simulation is ", trim(datafile)
  
END SUBROUTINE CREATEFILE

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARRAYS()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus
 
  ALLOCATE(rxyz(1:totblobs,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate rxyz"
  ALLOCATE(charge(1:totblobs), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate charge"
  ALLOCATE(aidvals(1:totblobs,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate aidvals"
  ALLOCATE(ixyz(1:totblobs,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate ixyz"
  ALLOCATE(topo_bond(1:2*nbonds,4), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate topo_bonds"

END SUBROUTINE ALLOCATE_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE CHECK_BOND_TYPES(chid,bid1,bid2)

  USE PARAMS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: chid, bid1, bid2
  INTEGER :: bchk, btyp_cnt

  ! Check if btype is present for angle generation by VMD
  DO bchk = bid1, bid2
     DO btyp_cnt = 1, numbondtypes
        IF(topo_bond(bchk,2) == btyp_cnt .AND. bflag_arr(btyp_cnt,1) &
             &== 0) THEN
           bflag_arr(btyp_cnt,1) = 1
           bflag_arr(btyp_cnt,2) = bchk ! First instance - bid
           bflag_arr(btyp_cnt,3) = chid ! First instance - chid
        END IF
     END DO
  END DO

END SUBROUTINE CHECK_BOND_TYPES
!--------------------------------------------------------------------

SUBROUTINE CROSS_CHECK_NUM_BOND_TYPES()
  
  USE PARAMS
  IMPLICIT NONE
  INTEGER :: i

  WRITE(outfile,*) "----------For debug----------------------------"
  DO i = 1, numbondtypes
     WRITE(outfile,*) i, bflag_arr(i,1), bflag_arr(i,2)
  END DO
 
  IF(SUM(bflag_arr(:,1)) .NE. numbondtypes) THEN
     WRITE(outfile,*) "Some bond types (mostly anion-anion) are not fo&
          &rmed. See above!"
  END IF

END SUBROUTINE CROSS_CHECK_NUM_BOND_TYPES
!--------------------------------------------------------------------

SUBROUTINE SANITY_CHECKS()

  USE PARAMS
  IMPLICIT NONE

  IF(ideal_an_per_ch < 1) THEN
     
     PRINT *, "ERROR: an_poly_rat in lmp_params.f90 is too small"
     PRINT *, "Ideal anions per chain: ", ideal_an_per_ch
     PRINT *, "Consider changing ideal_an_per_ch"
     STOP

  END IF
     
END SUBROUTINE SANITY_CHECKS
!--------------------------------------------------------------------

SUBROUTINE CLEAN_AND_CLOSE_ALL()

  USE PARAMS
  IMPLICIT NONE

  WRITE(outfile,*) "Successfully Written Datafile .."
  CLOSE (unit = outfile)

END SUBROUTINE CLEAN_AND_CLOSE_ALL
!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE PARAMS
  IMPLICIT NONE

  DEALLOCATE(rxyz)
  DEALLOCATE(ixyz)
  DEALLOCATE(charge)
  DEALLOCATE(aidvals)

END SUBROUTINE DEALLOCATE_ARRAYS

!--------------------------------------------------------------------
