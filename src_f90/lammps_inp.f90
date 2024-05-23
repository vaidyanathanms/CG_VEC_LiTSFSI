! To generate LAMMPS input file for polyelectrolyte simulations
! Use in conjunction with lmp_params.f90 and ran_numbers.f90
PROGRAM LAMMPSINP

  USE PARAMS

  IMPLICIT NONE
  
  LOGICAL :: input_coor = .false.
  REAL :: bondl, bondlsq
  INTEGER :: ierror,narg

  bondlsq = 0

  CALL SYSTEM_CLOCK(S)

!!$  narg = IARGC()
  
!!$  IF(narg == 1) THEN

  OPEN (unit = outfile, file = "lmp_input.dat", status ="replace",&
       & action="write",iostat=ierror)
  
  IF(ierror /= 0) STOP "Cannot open lmp_input.txt"
  
  CALL CREATEFILE() !CREATEFILE(narg)
  CALL COMPUTE_BOX()
  CALL ALLOCATE_ARRAYS()
  CALL INPCOR()
  CALL LMP_COORD()
  CALL DEALLOCATE_ARRAYS()
  WRITE(outfile,*) "Successfully Written Datafile .."
  CLOSE (unit = outfile)

!!$  END IF

END PROGRAM LAMMPSINP

!--------------------------------------------------------------------

SUBROUTINE LMP_COORD()

  USE PARAMS
  
  IMPLICIT NONE
  
  INTEGER :: i,j, k, ierror
  INTEGER ::  bondid, anglid, dihdid
  REAL ::  massval, rx, ry, rz
  i = 1
  
  PRINT *, "Writing LAMMPS Datafile .. "

20 FORMAT(5X,I0,2X,A)
22 FORMAT(5X,I0,2X,A)
24 FORMAT(5X,I0,2X,F14.6,2X,A)
  
  OPEN (unit=10, file = datafile, status="replace",action=&
       &"write",iostat = ierror)
  
  IF(ierror /= 0) STOP "Failed to open datafile"
     
  WRITE (10,*) "Data for CG p(mVEC-r-nLISTSFI) simulations"
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
  WRITE (10,20) numbondtypes, "bond types"
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
  
!!$Fetching image information

  ixyz = 0 ! If I don't assign, errors may creep
  
  DO i = 1,N_poly
     
     k = (i-1)*blob_per_ch + 1

     ixyz(k,1) = 0
     ixyz(k,2) = 0
     ixyz(k,3) = 0
     
     DO j = 1,blob_per_ch-1
        
        k = (i-1)*blob_per_ch + j
        
        IF (aidvals(k+1,3) == 2) THEN

           rx = rxyz(k,1) - rxyz(k+1,1)
           ry = rxyz(k,2) - rxyz(k+1,2)
           rz = rxyz(k,3) - rxyz(k+1,3)
           CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
           CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
           CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))

        ELSE

           IF (aidvals(k,3) == 1 .OR. aidvals(k,3) == CG_per_mon+1) THEN

              rx = rxyz(k,1) - rxyz(k+1,1)
              ry = rxyz(k,2) - rxyz(k+1,2)
              rz = rxyz(k,3) - rxyz(k+1,3)
              CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
              CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
              CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))

             
           ELSEIF (aidvals(k,2) == 2) THEN

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

  ! Writing Masses

  DO i = 1,numatomtypes

     massval = 1.000
     WRITE(10,'(I0,1X,F14.8)') i, massval

  END DO
  
  ! Writing atomic corrdinates
  
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

        IF(ixyz(i,1) .NE. 0 .OR. ixyz(i,2) .NE. 0 .OR. ixyz(i,3) .NE.&
             & 0) THEN
           PRINT *, "Long bonds found", i, ixyz(i,1), ixyz(i,2),&
                & ixyz(i,3)
           STOP
        END IF
        
     END IF

     
  END DO

  IF(numbondtypes /= 0) THEN

     ! Writing Bond Details  
     
     bondid = 0
     WRITE (10,*)
     WRITE (10,*) "Bonds"
     WRITE (10,*)
     
     DO i = 1,nbonds
        
        WRITE(10,'(4(I0,2X))') topo_bond(i,1), topo_bond(i,2),&
             & topo_bond(i,3), topo_bond(i,4)

     END DO
        
  END IF

!!$  IF(numangltypes /= 0) THEN
!!$
!!$     ! Writing Angle Details
!!$
!!$     anglid = 0
!!$     WRITE (10,*)
!!$     WRITE (10,*) "Angles"
!!$     WRITE (10,*)
!!$     
!!$     DO i = 1,N_brush
!!$        
!!$        DO j = 1,M_brush-2
!!$           
!!$           anglid = anglid + 1
!!$           k = (i-1)*M_brush + j           
!!$           WRITE(10,'(5(I0,2X))') anglid, angltype, aidvals(k,1)&
!!$                &,aidvals(k+1,1),aidvals(k+2,1)
!!$           
!!$        END DO
!!$        
!!$     END DO
!!$
!!$     DO i = 1,N
!!$        
!!$        DO j = 1,M-2
!!$           
!!$           anglid = anglid + 1
!!$           k = (i-1)*M + j + N_brush*M_brush
!!$           WRITE(10,'(5(I0,2X))') anglid, angltype, aidvals(k,1)&
!!$                &,aidvals(k+1,1),aidvals(k+2,1)
!!$           
!!$        END DO
!!$        
!!$     END DO
!!$
!!$  END IF
!!$
!!$  IF(numdihdtypes /= 0) THEN
!!$
!!$     ! Writing Dihedral Details
!!$     
!!$     dihdid = 0
!!$     WRITE (10,*)
!!$     WRITE (10,*) "Dihedrals"
!!$     WRITE (10,*)
!!$     
!!$     DO i = 1,N_brush
!!$        
!!$        DO j = 1,M_brush-3
!!$           
!!$           dihdid = dihdid + 1
!!$           k = (i-1)*M_brush + j 
!!$
!!$           WRITE(10,'(6(I0,2X))') dihdid, dihdtype, aidvals(k,1)&
!!$                &,aidvals(k+1,1), aidvals(k+2,1), aidvals(k+3,1)
!!$           
!!$        END DO
!!$        
!!$     END DO
!!$
!!$     DO i = 1,N
!!$        
!!$        DO j = 1,M-3
!!$           
!!$           dihdid = dihdid + 1
!!$           k = (i-1)*M + j + N_brush*M_brush
!!$
!!$           WRITE(10,'(6(I0,2X))') dihdid, dihdtype, aidvals(k,1)&
!!$                &,aidvals(k+1,1), aidvals(k+2,1), aidvals(k+3,1)
!!$           
!!$        END DO
!!$        
!!$     END DO
!!$
!!$  END IF

  CLOSE(unit = 10)

END SUBROUTINE LMP_COORD

!--------------------------------------------------------------------  

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

  volbox = REAL(totpart)/REAL(density)
  boxl_x = volbox**(1.0/3.0)
  boxl_y = volbox**(1.0/3.0)
  boxl_z = volbox**(1.0/3.0)

  ! Extra 2.0 factor above so that diagonal is approximately equal to
  !  sqrt(area/num_chains)

  WRITE(outfile,*) "----------System level details-------------------"
  WRITE(outfile,*) "Total particles: ", totpart
  WRITE(outfile,*) "Total blobs: ", totblobs
  WRITE(outfile,*) "Total # of polymer chains: ", N_poly
  WRITE(outfile,*) "Total # of anion blobs: ", N_anions
  WRITE(outfile,*) "Total # of cations: ", N_cations
  
  WRITE(outfile,*) "----------Chain level details--------------------"
  WRITE(outfile,*) "# of CG blob per polymer monomer: ", CG_per_mon
  WRITE(outfile,*) "Ratio between anions and VEC mons per chain: ",&
       & an_poly_rat
  WRITE(outfile,*) "# of blobs per chain: ", ideal_an_per_ch +&
       & CG_per_mon*M_poly
  WRITE(outfile,*) "# of anion blobs per chain: ", N_anions

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
  INTEGER :: an_per_ch,bondid_new, bondtemp, cgcnt, poly_mon
  REAL, PARAMETER :: r0init  = 0.97
  REAL, PARAMETER :: r0sq3   = r0init/sqrt(3.0)
  REAL, PARAMETER :: rmaxsq  = r0init*r0init
  REAL, PARAMETER :: math_pi = 3.14159265359
  REAL :: theta, phi
  REAL :: rx, ry, rz,rval
  LOGICAL :: in_box
  REAL :: X(3)
  REAL :: csum
  LOGICAL :: arr_bound
  
!!$  CALL RAN_INIT(S,X)
  CALL RANDOM_INIT(.TRUE., .TRUE.)
  
  WRITE(outfile,*) "Random Initial Configuration : NRRW"

! Create polymeric chains

  WRITE(outfile,*) "Generating chain configurations .. "

  i = 1; bondid_new = 0
  
  DO WHILE (i .LE. N_poly)

     bondtemp = bondid_new
     an_per_ch = 0
     cgcnt = 1
     j = 1
     poly_mon = 1
     
     ! First blob is VEC
     k = (i-1)*blob_per_ch + 1     
!!$     theta     = math_pi*RAN1(X)
!!$     phi       = 2*math_pi*RAN1(X)
     theta     = math_pi*RAND()
     phi       = 2*math_pi*RAND()
     rxyz(k,1) = r0init*sin(theta)*cos(phi)
     rxyz(k,2) = r0init*sin(theta)*sin(phi)
     rxyz(k,3) = r0init*cos(theta)

     ! Create atom types here
     aidvals(k,1) = k ! Atom ID
     aidvals(k,2) = i ! Mol ID
     aidvals(k,3) = 1 ! Atom type
     charge(k)    = (REAL(cgcnt)-0.5*(1+CG_per_mon))*2*charge_poly
     
     DO cgcnt = 2, CG_per_mon

        k = (i-1)*blob_per_ch + cgcnt
        bondtemp     = bondtemp + 1
        CALL CHECK_ARRAY_BOUNDS(i,k,bondtemp,arr_bound)
!!$        theta     = math_pi*RAN1(X)
!!$        phi       = 2*math_pi*RAN1(X)
        theta     = math_pi*RAND()
        phi       = 2*math_pi*RAND()
        rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
        rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
        rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
        aidvals(k,1) = k
        aidvals(k,2) = i
        aidvals(k,3) = cgcnt
        charge(k)    = (REAL(cgcnt)-0.5*(1+CG_per_mon))*2*charge_poly
        CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1),aidvals(k,1)&
             &,160)
        
     END DO

     j = 2
     DO WHILE (j .LE. M_poly+ideal_an_per_ch) !j runs over MONOMERS
        !and not BLOBS

        ! Second monomer onwards can be STSFI (anion)
        IF (j .NE. M_poly+ideal_an_per_ch .AND. RAND() .LE.&
             & an_poly_rat) THEN ! can be attached to
           ! blob-1 only
!!$                 theta     = math_pi*RAN1(X)
!!$                 phi       = 2*math_pi*RAN1(X)
           bondtemp       = bondtemp + 1
           ! Check if more blobs/bonds are created than the array
           ! bounds
           CALL CHECK_ARRAY_BOUNDS(i,k+1,bondtemp,arr_bound)

           IF (arr_bound .EQV. .TRUE.) THEN
              theta     = math_pi*RAND()
              phi       = 2*math_pi*RAND()          
              aidvals(k+1,1) = k + 1
              aidvals(k+1,2) = i
              aidvals(k+1,3) = CG_per_mon + 1
              charge(k+1)    = -1.0
           
              IF(aidvals(k,3) == CG_per_mon + 1) THEN !anion-anion
                 rxyz(k+1,1) = rxyz(k,1) + r0init*sin(theta)*cos(phi)
                 rxyz(k+1,2) = rxyz(k,2) + r0init*sin(theta)*sin(phi)
                 rxyz(k+1,3) = rxyz(k,3) + r0init*cos(theta)
                 CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k+1,1)&
                      &,aidvals(k,1),100)

              ELSEIF(aidvals(k,3) == 1) THEN ! Wrong choice
                 
                 PRINT *, "WRONG CHOICE" , k
                 STOP
                 
              ELSE !poly_mon - anion bond
                 rxyz(k+1,1) = rxyz(k+1-CG_per_mon,1) + r0init&
                      &*sin(theta)*cos(phi)
                 rxyz(k+1,2) = rxyz(k+1-CG_per_mon,2) + r0init&
                      &*sin(theta)*sin(phi)
                 rxyz(k+1,3) = rxyz(k+1-CG_per_mon,3) + r0init&
                      &*cos(theta)
                 CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k+1,1)&
                      &,aidvals(k+1-CG_per_mon,1),100)

              END IF
              
              an_per_ch = an_per_ch + 1 ! update an_per_ch
              k = k + 1 ! pseudo update if anion-anion bond is formed

           ELSE

              j = M_poly+ideal_an_per_ch

           END IF
              
        ELSE
           
           !CG Blobs
           cgcnt = 1
           DO WHILE (cgcnt .LE. CG_per_mon)
             
              k = (i-1)*blob_per_ch + CG_per_mon*poly_mon + cgcnt + an_per_ch
              bondtemp     = bondtemp + 1
              ! Check if more blobs/bonds are created than the array
              ! bounds
              CALL CHECK_ARRAY_BOUNDS(i,k,bondtemp,arr_bound)

              IF (arr_bound .EQV. .TRUE.) THEN
!!$           theta     = math_pi*RAN1(X)
!!$           phi       = 2*math_pi*RAN1(X)
                 theta     = math_pi*RAND()
                 phi       = 2*math_pi*RAND()
                 aidvals(k,1) = k
                 aidvals(k,2) = i
                 aidvals(k,3) = cgcnt
                 charge(k)    = (REAL(cgcnt)-0.5*(1+CG_per_mon))*2&
                      &*charge_poly
                 
                 IF (cgcnt == 1) THEN
                    IF(aidvals(k-1,3) == CG_per_mon + 1) THEN !poly_mon-anion
                       rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)&
                            &*cos(phi)
                       rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)&
                            &*sin(phi)
                       rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
                       CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1)&
                            &,aidvals(k,1),200)
                    ELSE !poly_mon_(k-1) - poly_mon_k
                       rxyz(k,1) = rxyz(k-CG_per_mon,1) + r0init&
                            &*sin(theta)*cos(phi)
                       rxyz(k,2) = rxyz(k-CG_per_mon,2) + r0init&
                            &*sin(theta)*sin(phi)
                       rxyz(k,3) = rxyz(k-CG_per_mon,3) + r0init&
                            &*cos(theta)
                       CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k&
                            &-CG_per_mon,1),aidvals(k,1),200)
                    END IF
                 ELSE !poly_mon_k - poly_mon_k (blob-blob connection)
                    rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)&
                         &*cos(phi)
                    rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)&
                         &*sin(phi)
                    rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
                    CALL ASSIGN_BOND_TOPO(bondtemp,aidvals(k-1,1)&
                         &,aidvals(k,1),200)
                 END IF
                 cgcnt = cgcnt + 1

              ELSE

                 j = M_poly+ideal_an_per_ch
                 cgcnt = CG_per_mon+1

              END IF
             
           END DO

           poly_mon = poly_mon + 1
           
        END IF

        j = j + 1
        
     END DO


     IF (an_per_ch == ideal_an_per_ch) THEN

        i = i + 1
        bondid_new = bondtemp

     END IF
     
  END DO

  PRINT *, "Ideal number of polymer blobs: ", N_poly*blob_per_ch
  PRINT *, "Total number of polymer blobs: ", k
  PRINT *, "NCations, NAnions ", N_cations, N_anions
! Create lithium cations

  DO i = k+1, k + N_cations
  
!!$     rxyz(i,1) = RAN1(X)*boxl_x
!!$     rxyz(i,2) = RAN1(X)*boxl_y
!!$     rxyz(i,3) = RAN1(X)*boxl_z
     rxyz(i,1) = RAND()*boxl_x
     rxyz(i,2) = RAND()*boxl_y
     rxyz(i,3) = RAND()*boxl_z
     aidvals(i,1) = i
     aidvals(i,2) = N_poly+1
     aidvals(i,3) = CG_per_mon + 2
     charge(i)    = 1.0
        
  END DO
     
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

SUBROUTINE ASSIGN_BOND_TOPO(bid,aid1,aid2,iderr)

  USE PARAMS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: bid, aid1, aid2,iderr
  INTEGER :: btype
  INTEGER :: atype1, atype2

  atype1 = aidvals(aid1,3)
  atype2 = aidvals(aid2,3)
  
  IF (atype1 == atype2) THEN

     IF (atype1 == 1) THEN
        btype = 1
     ELSEIF (atype1 == 3) THEN
        btype = 2
     ELSE
        PRINT *, "Error type 1", bid, aid1, aid2,atype1,atype2,iderr
     END IF
        
  ELSEIF (atype1 .NE. 2 .AND. atype2 .NE. 2) THEN
     btype = 3

  ELSEIF((atype1 == 1 .AND. atype2 == 2) .OR. (atype1==2 .AND. atype2&
       & == 1)) THEN
     btype = 4
  ELSE
     PRINT *, "Error type 2", bid, aid1, aid2,atype1,atype2,iderr
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

!!$  INTEGER,INTENT(IN)  :: narg
!!$  CHARACTER (LEN = 7) :: ch_arch
  CHARACTER (LEN = 6) :: nmon_char
  CHARACTER (LEN = 7) :: prefix = "VECdata"
 
  WRITE(nmon_char,"(I0)") M_poly

  datafile = prefix//trim(adjustl(nmon_char))//'.dat'
     
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
  ALLOCATE(topo_bond(1:nbonds,4), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate topo_bonds"


END SUBROUTINE ALLOCATE_ARRAYS

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
