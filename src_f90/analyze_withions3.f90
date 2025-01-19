!********************************************************************
!-------------------------Main Program File--------------------------
!--------------------------Ver:Oct-04-2016---------------------------
!-------For analyzing static properties of PS-PEO with ions system---
!--------------------Parameter file - ionparameters.f90--------------
!--------------------Analysis input file - outputfile.txt------------
!-----------Calculates RDF/Clusters/Complexed-Free O RDF-------------
!-----------Calculates Layer wise properties ------------------------
!-----------Definition of subroutines in the parameter file----------
!********************************************************************

PROGRAM MAIN

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  PRINT *, "Static analysis of PS-PEO system with ions .."
  PRINT *, "Starting OMP Threads .."

!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
  PRINT *, "Number of threads are: ", nproc
!$OMP END PARALLEL

  CALL READ_ANA_IP_FILE()
  CALL READFILE()
  CALL ANALYZE_PSPEO_WITHIONS_TRAJECTORYFILE()
  CALL FINDINTERFACES()
  CALL ALLOUTPUTS()
  IF(all_tana) CALL SPATIAL_ANALYSIS()
  IF(rdf_dom_ana) CALL OUTPUT_ALLCROSSRDF()
  CALL DEALLOCATE_ARRAYS()

  PRINT *, "All Calculations Completed Succesfully :)"

END PROGRAM MAIN

!--------------------------------------------------------------------

SUBROUTINE READ_ANA_IP_FILE()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag
  CHARACTER(256) :: dumchar

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"

  logflag = 0

  CALL GETARG(nargs,ana_fname)

  OPEN(unit = anaread,file=trim(ana_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(ana_fname), "not found"
     STOP

  END IF

  DO

     READ(anaread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     IF(dumchar == 'datafile') THEN
        
        READ(anaread,*,iostat=ierr) data_fname

     ELSEIF(dumchar == 'trajectory_file') THEN

        READ(anaread,*,iostat=ierr) traj_fname

     ELSEIF(dumchar == 'interface_file') THEN

        READ(anaread,*,iostat=ierr) inter_fname

     ELSEIF(dumchar == 'rdf_file') THEN

        READ(anaread,*,iostat=ierr) rdf_fname
        
     ELSEIF(dumchar == 'oxy_file') THEN

        READ(anaread,*,iostat=ierr) oxy_fname

     ELSEIF(dumchar == 'clust_file') THEN

        READ(anaread,*,iostat=ierr) clust_fname

     ELSEIF(dumchar == 'clust_dom_file') THEN

        READ(anaread,*,iostat=ierr) clust_fname

     ELSEIF(dumchar == 'dens_file') THEN

        READ(anaread,*,iostat=ierr) dens_fname

     ELSEIF(dumchar == 'lipclust_file') THEN

        READ(anaread,*,iostat=ierr) lipclust_fname

     ELSEIF(dumchar == 'rdf2d_file') THEN

        READ(anaread,*,iostat=ierr) rdf2d_fname
        
     ELSEIF(dumchar == 'rdf2d_dom_file') THEN

        READ(anaread,*,iostat=ierr) rdf2d_dom_fname

     ELSEIF(dumchar == 'crossrdf_file') THEN

        READ(anaread,*,iostat=ierr) crossrdf_fname

     ELSEIF(dumchar == 'oxyspatcoord_file') THEN

        READ(anaread,*,iostat=ierr) lioxypos_file

     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword", trim(dumchar)
        STOP

     END IF

  END DO

  IF(logflag == 0) STOP "no logfile name found"
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  PRINT *, "Analysis input file read finished .."

END SUBROUTINE READ_ANA_IP_FILE

!--------------------------------------------------------------------

SUBROUTINE READFILE()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,k,ierr,u,AllocateStatus
  INTEGER :: flag, cntr, nwords
  INTEGER :: aid,molid,atype,ix,iy,iz
  REAL    :: charge,rx,ry,rz
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(256) :: rline,dumchar


  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"

  WRITE(logout,*) "Datafile used is :", trim(adjustl(data_fname))

  OPEN(unit=ionwrite,file = 'iondistdata.txt',action ='write', status&
       &='replace',iostat=ierr) 

  
  ntotatoms = 0;ntotbonds=0;ntotangls=0;ntotdihds=0;ntotimprs=0
  atomflag =0;velflag = 0;bondflag=0;anglflag=0;dihdflag=0;imprflag=0

  READ(inpread,*)
  READ(inpread,*)

  DO i = 1,8 !Change here according to convenience
       
     READ(inpread,*) u, dumchar
!!$     CALL COUNTWORDS(adjustl(trim(rline)), nwords)
!!$
!!$     IF(nwords /= 0) THEN
!!$
!!$        ALLOCATE(chararr(256,nwords),stat = AllocateStatus)
!!$        IF(AllocateStatus/=0) STOP "did not allocate chararr"
!!$
!!$        CALL PARSECHAR(adjustl(trim(rline)), nwords)
!!$     
!!$        pause

        
        IF(dumchar == "atoms") THEN
           ntotatoms = u
        ELSEIF(dumchar == "bonds") THEN
           ntotbonds = u
        ELSEIF(dumchar == "angles") THEN
           ntotangls = u
        ELSEIF(dumchar == "dihedrals") THEN
           ntotdihds = u
        ELSEIF(dumchar == "atom" .OR. dumchar == "atomtypes") THEN
           ntotatomtypes = u
        ELSEIF(dumchar == "bond" .OR. dumchar == "bondtypes") THEN
           ntotbondtypes = u
        ELSEIF(dumchar == "angle" .OR. dumchar == "atomtypes") THEN
           ntotangltypes = u
        ELSEIF(dumchar == "dihedral" .OR. dumchar == "dihedraltypes") THEN
           ntotdihdtypes = u
        ELSEIF(dumchar == "improper" .OR. dumchar == "impropertypes") THEN
           ntotimprtypes = u
        ELSEIF(dumchar == "Masses") THEN
           
           ALLOCATE(masses(ntotatomtypes,1),stat = AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate masses"
           
           DO j = 1,ntotatomtypes
              
              READ(inpread,*) u, masses(u,1)
              
           END DO
           
        END IF
        
!!$        DEALLOCATE(chararr)
!!$
!!$     END IF

  END DO

  READ(inpread,*) 
  READ(inpread,*) xlo, xhi
  READ(inpread,*) ylo, yhi
  READ(inpread,*) zlo, zhi
  
  box_xl = xhi - xlo
  box_yl = yhi - ylo
  box_zl = zhi - zlo

  IF(major_axis == 1) boxval = box_xl
  IF(major_axis == 2) boxval = box_yl
  IF(major_axis == 3) boxval = box_zl

  PRINT *, "x-box  ", "y-box  ", "z-box  "
  PRINT *, box_xl, box_yl, box_zl

  PRINT *, "STATISTICS WITH IONS"
  PRINT *, "Number of atoms/atomtypes: " , ntotatoms,ntotatomtypes
  PRINT *, "Number of bonds/bondtypes: " , ntotbonds,ntotbondtypes
  PRINT *, "Number of angles/angletypes: " , ntotangls,ntotangltypes
  PRINT *, "Number of diheds/dihedtypes: " , ntotdihds,ntotdihdtypes
  flag = 0; cntr = 0

  CALL ALLOCATE_ARRAYS()

  DO 

     READ(inpread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     !READ DATA HERE FOR CHARGES AND MOLID
     !READ EVERYTHING AND OVERWRITE LATER
     IF(trim(dumchar) == "Atoms") THEN
             
        atomflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms

           READ(inpread,*) aid,molid,atype,charge,rx,ry,rz

           rx = rx - xlo
           ry = ry - ylo
           rz = rz - zlo

           aidvals(aid,1)     = aid
           aidvals(aid,2)     = molid
           aidvals(aid,3)     = atype
           charge_lmp(aid,1)  = charge
           rxyz_lmp(aid,1)    = rx
           rxyz_lmp(aid,2)    = ry
           rxyz_lmp(aid,3)    = rz
!!$           image_lmp(aid,1)   = ix
!!$           image_lmp(aid,2)   = iy
!!$           image_lmp(aid,3)   = iz

        END DO

        npolymol   = nchains
        PRINT *, "Number of polymer molecules : ", npolymol

     END IF

     IF(trim(dumchar) == "Masses") THEN

        ALLOCATE(masses(ntotatomtypes,1),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate masses"
         
        DO j = 1,ntotatomtypes

           READ(inpread,*) u, masses(u,1)

        END DO

     END IF

     IF(trim(dumchar) == "Velocities") THEN
             
        velflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms

           READ(inpread,*) vel_xyz(j,1),vel_xyz(j,2),vel_xyz(j,3)&
                &,vel_xyz(j,4)

        END DO


     END IF

     IF(trim(dumchar) == "Bonds") THEN
             
        bondflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotbonds

           READ(inpread,*) bond_lmp(j,1),bond_lmp(j,2),bond_lmp(j,3)&
                &,bond_lmp(j,4)

        END DO

     END IF

     IF(trim(dumchar) == "Angles") THEN
             
        anglflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotangls

           READ(inpread,*) angl_lmp(j,1),angl_lmp(j,2),angl_lmp(j,3)&
                &,angl_lmp(j,4),angl_lmp(j,5)

        END DO

     END IF

     IF(trim(dumchar) == "Dihedrals") THEN
             
        dihdflag = 1
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ntotdihds

           READ(inpread,*) dihd_lmp(j,1),dihd_lmp(j,2),dihd_lmp(j,3)&
                &,dihd_lmp(j,4),dihd_lmp(j,5),dihd_lmp(j,6)

        END DO

     END IF
  
     IF(trim(dumchar) == "Impropers") THEN
             
        imprflag = 1
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ntotimprs

           READ(inpread,*) impr_lmp(j,1),impr_lmp(j,2),impr_lmp(j,3)&
                &,impr_lmp(j,4),impr_lmp(j,5),impr_lmp(j,6)

        END DO

     END IF

  END DO
  
  PRINT *, "Fileread finish .."

  CALL COARSEGRAIN_PF6()

  IF(initdist) THEN
     
     CALL CALCDENSPROFILES()
     CALL OUTPUTDENS(1)

  END IF

END SUBROUTINE READFILE

!--------------------------------------------------------------------

!!$SUBROUTINE COUNTWORDS(linechar,nwords)
!!$
!!$  USE ALLCHAR
!!$  IMPLICIT NONE
!!$
!!$  CHARACTER(LEN=*), INTENT(IN) :: linechar
!!$  INTEGER, INTENT(OUT)  :: nwords
!!$
!!$  INTEGER :: pos, i
!!$
!!$  pos = 1
!!$  nwords = 0
!!$
!!$  DO 
!!$
!!$     i = VERIFY(linechar(pos:),' ')
!!$     IF(i == 0) EXIT
!!$     
!!$     nwords = nwords + 1
!!$     pos = pos + i - 1
!!$     
!!$     i = SCAN(linechar(pos:),' ')
!!$     IF(i == 0) EXIT
!!$     
!!$     pos = pos + i - 1
!!$     
!!$  END DO
!!$
!!$END SUBROUTINE COUNTWORDS

!--------------------------------------------------------------------

!!$SUBROUTINE PARSECHAR(linechar,nwords)
!!$  
!!$  USE ALLCHAR
!!$  IMPLICIT NONE
!!$
!!$  INTEGER, INTENT(IN) :: nwords
!!$  CHARACTER(LEN=*), INTENT(IN) :: linechar
!!$  INTEGER :: i, pos1, pos2, n,j
!!$  CHARACTER,DIMENSION(1:256) :: string
!!$
!!$  pos1 = 1; n = 0
!!$  
!!$  DO 
!!$
!!$     pos2 = INDEX(linechar(pos1:),' ')
!!$
!!$     IF(pos2 == 0) THEN
!!$        
!!$        n = n + 1
!!$
!!$        DO i = 1,len(linechar)
!!$
!!$           keywords(i,n) = linechar(i:i)
!!$
!!$        END DO
!!$
!!$        EXIT 
!!$
!!$     ELSE
!!$
!!$        n = n + 1
!!$
!!$        j = 0
!!$
!!$        DO i = pos1,pos1+pos2-2
!!$
!!$           j = j + 1
!!$           keywords(j,n) = linechar(i:i)
!!$
!!$        END DO
!!$
!!$        pos1 = pos2 + pos1
!!$
!!$     END IF
!!$
!!$  END DO
!!$
!!$     
!!$END SUBROUTINE PARSECHAR
!!$
!!$!--------------------------------------------------------------------

SUBROUTINE ANALYZE_PSPEO_WITHIONS_TRAJECTORYFILE()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,aid,ierr,atchk,atype
  REAL :: xlo,xhi,ylo,yhi,zlo,zhi

  OPEN(unit = 15,file =trim(traj_fname),action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "trajectory file not found"

  WRITE(logout,*) "trajectory file used is :",&
       & trim(adjustl(traj_fname))

  binavg = 0.0; volavg = 0.0; boxlzavg = 0.0
  adensavg = 0.0; bdensavg = 0.0; middensavg = 0.0; lithdensavg = 0.0
  pf6densavg = 0.0

  OPEN(unit = rgwrite,file ="rganalysis.lammpstrj",action="write"&
       &,status="replace")

  OPEN(unit = rgavgwrite,file ="rgavg.lammpstrj",action="write"&
       &,status="replace")

  DO i = 1,skipfr

     DO j = 1,ntotatoms+9

        READ(15,*) 

     END DO

     IF(mod(i,100) == 0) PRINT *, "Skipped ", i, "frames"

  END DO

  DO i = 1,nframes

     IF(i == 1)   WRITE(fwnum,'(I10)') int(i)
     IF(mod(i,50) == 0) PRINT *, "Processing ", i,"th frame"

     READ(15,*)
     READ(15,*) timestep

     READ(15,*) 
     READ(15,*) atchk
     IF(atchk /= ntotatoms) STOP "Number of atoms do not match"

     READ(15,*) 
     READ(15,*) xlo, xhi
     READ(15,*) ylo, yhi
     READ(15,*) zlo, zhi

     READ(15,*)

     box_xl = xhi - xlo
     box_yl = yhi - ylo
     box_zl = zhi - zlo
     
     IF(major_axis == 1) boxval = box_xl
     IF(major_axis == 2) boxval = box_yl
     IF(major_axis == 3) boxval = box_zl

     boxx_arr(i)  = box_xl
     boxy_arr(i)  = box_yl
     boxz_arr(i)  = box_zl
     boxvalarr(i) = boxval

     DO j = 1,ntotatoms

        READ(15,*) aid,atype,rxyz_lmp(aid,1),rxyz_lmp(aid,2)&
             &,rxyz_lmp(aid,3)

        IF(atype .NE. aidvals(aid,3)) THEN

           PRINT *, "Incorrect atom ids"
           PRINT *, i,j,aid,atype,aidvals(aid,3)
           STOP

        END IF

     END DO

!$OMP PARALLEL FIRSTPRIVATE(i) PRIVATE(j)
!$OMP DO
     DO j = 1,ntotatoms

        rxyz_lmp(j,1) = rxyz_lmp(j,1) - xlo
        rxyz_lmp(j,2) = rxyz_lmp(j,2) - ylo
        rxyz_lmp(j,3) = rxyz_lmp(j,3) - zlo
        
        IF(all_tana) THEN
           
           trx_lmp(j,i) = rxyz_lmp(j,1)
           try_lmp(j,i) = rxyz_lmp(j,2)
           trz_lmp(j,i) = rxyz_lmp(j,3)
           
        END IF

     END DO
!$OMP END DO        
!$OMP END PARALLEL

     IF(i == 1) CALL SORTALLARRAYS()
     CALL CALCDENSPROFILES()
     CALL STRUCTANALYSIS(i)


! Test

!!$     IF(i == 1) THEN
!!$
!!$        OPEN(unit = 91,file ="chkcoord.txt",action="write",status&
!!$             &="replace")
!!$        
!!$        DO j = 1,ntotatoms
!!$
!!$           WRITE(91,'(2(I0,1X),3(F16.9,1X))') j, aidvals(j),&
!!$                & rxyz_lmp(j,1), rxyz_lmp(j,2), rxyz_lmp(j,3)
!!$
!!$
!!$        END DO
!!$
!!$        CLOSE(91)
!!$
!!$     END IF

  END DO

  IF(lioxysort) CALL BINDING_LOCATION

  CALL OUTPUTDENS(0)

  CLOSE(15)

END SUBROUTINE ANALYZE_PSPEO_WITHIONS_TRAJECTORYFILE

!--------------------------------------------------------------------

SUBROUTINE CALCDENSPROFILES()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,ibin,acntr,bcntr,midcntr,iflag,interf,atype,dval
  INTEGER :: lithcntr,pf6cntr
  REAL,DIMENSION(0:maxbin-1) :: adens, bdens, middens,lithdens,pf6dens
  REAL :: binval,rval,aval,bval,volval,vnorm
  REAL :: minval,min2val,maxval,max2val
  INTEGER, DIMENSION(1:4) :: lidomcnt, pf6domcnt

  lidomcnt = 0; pf6domcnt = 0
  adens = 0.0
  bdens = 0.0
  middens = 0.0
  lithdens = 0.0
  pf6dens = 0.0

  binval = boxval/maxbin
  binavg = binavg + binval
  volval = box_xl*box_yl*box_zl
  volavg = volavg + volval
  boxlzavg = boxlzavg + boxval
  vnorm = boxval/(volval*binval)

  acntr = 0; bcntr = 0; midcntr = 0;lithcntr = 0;pf6cntr=0

!Find minimum position

!!$  IF(major_axis == 1) minval = rxyz_lmp(1,1)
!!$  IF(major_axis == 2) minval = rxyz_lmp(1,2)
!!$  IF(major_axis == 3) minval = rxyz_lmp(1,3)
!!$  
!!$  maxval = minval
!!$
!!$  DO i = 2,ntotatoms
!!$     
!!$     IF(major_axis == 1) min2val = rxyz_lmp(i,1)
!!$     IF(major_axis == 2) min2val = rxyz_lmp(i,2)
!!$     IF(major_axis == 3) min2val = rxyz_lmp(i,3)
!!$
!!$     max2val = min2val
!!$
!!$     IF(min2val .LT. minval) minval = min2val
!!$     IF(max2val .GT. maxval) maxval = max2val
!!$
!!$  END DO

!$OMP PARALLEL PRIVATE(i,rval,atype,ibin,dval) 
!$OMP DO REDUCTION(+:lithdens,adens,bdens,middens,pf6dens,lidomcnt&
!$OMP&  ,pf6domcnt,acntr,bcntr,midcntr,lithcntr,pf6cntr)

  DO i = 1,ntotatoms

     IF(major_axis == 1) rval = rxyz_lmp(i,1) 
     IF(major_axis == 2) rval = rxyz_lmp(i,2) 
     IF(major_axis == 3) rval = rxyz_lmp(i,3) 
     
     atype = aidvals(i,3)
     IF(atype == 1 .OR. atype == 2 .OR. atype == 3 .OR. atype == 4&
          & .OR. atype == 5) THEN
        
        !rval = rval - minval
        rval = rval - boxval*floor(rval/boxval)
        ibin = FLOOR(rval/binval)

        IF(ibin .LT. maxbin) THEN

           adens(ibin) = adens(ibin) + masses(atype,1)
           acntr  = acntr + 1

        END IF
        
     ELSEIF(atype == 6 .OR. atype == 7 .OR. atype == 8  .OR. atype ==&
          & 9) THEN
           

        !rval = rval - minval
        rval = rval - boxval*floor(rval/boxval)
        ibin = FLOOR(rval/binval)

        IF(ibin .LT. maxbin) THEN

           bdens(ibin) = bdens(ibin) +  masses(atype,1)
           bcntr  = bcntr + 1

        END IF

     ELSEIF(atype == 10) THEN

        !rval = rval - minval
        rval = rval - boxval*floor(rval/boxval)
        ibin = FLOOR(rval/binval)

        IF(ibin .LT. maxbin) THEN

           lithdens(ibin) = lithdens(ibin) + masses(atype,1)
           lithcntr  = lithcntr + 1

        END IF

        CALL FINDDOMAIN(rval,boxval,dval)
        lidomcnt(dval) = lidomcnt(dval) + 1

     ELSEIF(atype == 11 .OR. atype == 12 .OR. atype == 13  .OR. atype&
          & == 14) THEN

                !rval = rval - minval
        rval = rval - boxval*floor(rval/boxval)
        ibin = FLOOR(rval/binval)

        IF(ibin .LT. maxbin) THEN

           pf6dens(ibin) = pf6dens(ibin) + masses(atype,1)
           pf6cntr  = pf6cntr + 1

        END IF
        
        CALL FINDDOMAIN(rval,boxval,dval)
        pf6domcnt(dval) = pf6domcnt(dval) + 1
        

     END IF


     IF(atype == midattype) THEN

        IF(major_axis == 1) rval = rxyz_lmp(i,1) 
        IF(major_axis == 2) rval = rxyz_lmp(i,2) 
        IF(major_axis == 3) rval = rxyz_lmp(i,3) 

        !rval = rval - minval
        rval = rval - boxval*floor(rval/boxval)
        ibin = FLOOR(rval/binval)

        IF(ibin .LT. maxbin) THEN

           middens(ibin) = middens(ibin) + masses(atype,1)
           midcntr  = midcntr + 1

        END IF

     END IF

  END DO

!$OMP END DO

  IF(acntr + bcntr + lithcntr + pf6cntr .NE. ntotatoms) THEN
     
     PRINT *, "Unequal atom count"
     PRINT *, acntr + bcntr + lithcntr + pf6cntr, ntotatoms
     STOP

  END IF

 ! print *, acntr, bcntr,midcntr, lithcntr,pf6cntr

!$OMP DO PRIVATE(i) 
  DO i = 0,maxbin-1

     adensavg(i) = adensavg(i) + adens(i)*vnorm
     bdensavg(i) = bdensavg(i) + bdens(i)*vnorm
     middensavg(i) = middensavg(i) + middens(i)*vnorm
     lithdensavg(i) = lithdensavg(i) + lithdens(i)*vnorm
     pf6densavg(i) = pf6densavg(i) + pf6dens(i)*vnorm

  END DO
!$OMP END DO
!$OMP END PARALLEL

  WRITE(ionwrite,*) lidomcnt
  WRITE(ionwrite,*) pf6domcnt

END SUBROUTINE CALCDENSPROFILES

!--------------------------------------------------------------------

SUBROUTINE FINDINTERFACES()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

!!$ To calculate the interface positions
  REAL*8  :: dum1, dum2, dum3, volnorm, duma, dumb
  INTEGER :: i,j,AllocateStatus, maxmid, maxdens,icount,iflag,jinit
  REAL, DIMENSION(0:maxbin-1)  :: hist_mavg, hist_aavg, hist_bavg
  REAL, DIMENSION(maxbin*maxinter) :: dummid_interarr, dumdens_interarr
  !Number of elements is arbitrary 

  dummid_interarr  = -1
  dumdens_interarr = -1

  OPEN(unit = intmface,file = trim(inter_fname), status="replace",&
       & action ="write")

  DO i = 0,maxbin -1
     
     hist_mavg(i)  = nchains*middensavg(i)/(2*densmult*REAL(nframes))
     hist_aavg(i)  = adensavg(i)
     hist_bavg(i)  = bdensavg(i)
 
  END DO

!Using the middle monomer distribution. 
  icount = 0
  
!Left side of the box
  j = 0

  dum1 = hist_mavg(maxbin - 1)
  dum2 = hist_mavg(j)
  dum3 = hist_mavg(j+1)
  
  IF(dum2 > dum1 .and. dum2 > dum3 .and. dum2 > 0.1) THEN
     
!The last condition is arbitrary to avoid small bumps near the minima
!of the monomer distribution curve        
     icount = icount + 1
     dummid_interarr(icount) = (j+0.5)*binavg 
     
  END IF

!Inside the box
  DO j = 1, maxbin-2
     
     dum1 = hist_mavg(j-1)
     dum2 = hist_mavg(j)
     dum3 = hist_mavg(j+1)
      
     IF(dum2 > dum1 .and. dum2 > dum3 .and. dum2 > 0.1) THEN

        icount = icount + 1
        dummid_interarr(icount) = (j+0.5)*binavg
             
     END IF

  END DO


!Right side of the box
  j = maxbin - 1

  dum1 = hist_mavg(j - 1)
  dum2 = hist_mavg(j)
  dum3 = hist_mavg(0)
  
  IF(dum2 > dum1 .and. dum2 > dum3 .and. dum2 > 0.1) THEN

     icount = icount + 1
     dummid_interarr(icount) = (j+0.5)*binavg
     
  END IF

  maxmid   = icount

  PRINT *, "Number of interfaces from midmon distribution ", maxmid
  PRINT *, dummid_interarr(1:maxmid)

!Using the density profiles
  icount = 0
  j = 0
  jinit = 1
  duma = hist_aavg(j)
  dumb = hist_bavg(j)
!Left most end
  IF(duma > dumb) THEN
     iflag = 1
  ELSEIF(dumb > duma) THEN
     iflag = 2
  ELSE
     icount = icount + 1
     dumdens_interarr(icount) = (j+0.25)*binavg
     duma = hist_aavg(j+1)
     dumb = hist_bavg(j+1)
     jinit = 2
     IF(duma > dumb) iflag = 1
     IF(dumb > duma) iflag = 2
     IF(duma == dumb) STOP "Wrong density profiles"
  END IF
!Inside the box  
  DO j = jinit, maxbin -1
                
                
     IF(iflag == 1) THEN

        IF(hist_aavg(j) .LE. hist_bavg(j)) THEN !Change sign
           
           icount = icount + 1
           dumdens_interarr(icount) = (j+0.5)*binavg
           iflag  = 2

        END IF

     ELSEIF(iflag == 2) THEN

        IF(hist_bavg(j) .LE. hist_aavg(j)) THEN !Change sign
           
           icount = icount + 1
           dumdens_interarr(icount) = (j+0.5)*binavg
           iflag  = 1

        END IF

     END IF

  END DO
!Right most end
  IF((hist_aavg(maxbin-1) .GE. hist_bavg(maxbin-1)) .AND.&
       & (hist_aavg(0) .LE. hist_bavg(0))) THEN
     
     icount = icount + 1
     dumdens_interarr(icount) = (maxbin - 1 + 0.75)*binavg

  END IF

  maxdens = icount

  PRINT *, "Number of interfaces obtained from density: ", maxdens
  PRINT *, dumdens_interarr(1:maxdens)

! Need to determine which one needs to be used.
! Allocate to actual array. This can keep things clean

  maxinter = MIN(maxdens, maxmid)

  ALLOCATE (inter(maxinter), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "***Allocation inter not proper***"

  ALLOCATE (domtyp(maxinter), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "***Allocation domtyp not proper***"

  IF(maxmid > maxdens) THEN

     WRITE(logout,*) "Using density profiles to write interfaces. Check&
          &midmondistribution to understand discrepancies. "

     DO i = 1, maxinter

        IF(dumdens_interarr(i) == -1) THEN
           
           PRINT *, "Array updated incorrectly" 
           STOP
           
        END IF
     
        inter(i) = dumdens_interarr(i)
        
     END DO

  ELSEIF(maxmid == maxdens) THEN
     
     WRITE(logout,*) "Equal number of interfaces obtained using both ana&
          &lysis, Using interface distribution to update interfaces" 

     DO i = 1, maxinter

        IF(dumdens_interarr(i) == -1) THEN
           
           PRINT *, "Array updated incorrectly" 
           STOP
           
        END IF
     
        inter(i) = dumdens_interarr(i)
        
     END DO

  ELSE

     WRITE(logout,*) "Strange occurence, maxmid < maxdens. Check o/p&
          & profiles "
     STOP

  END IF

  WRITE(intmface,*) maxinter
  WRITE (logout,*) "Number of interfaces: ", maxinter
  WRITE (logout,*) inter

  DO i = 1, maxinter

     WRITE(intmface,*) i, inter(i)

  END DO



  CLOSE(intmface)

END SUBROUTINE FINDINTERFACES

!--------------------------------------------------------------------

SUBROUTINE INTPOSFILE()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i, ierror, val,AllocateStatus

  OPEN(unit = intmface,file = trim(inter_fname),status="old", action&
       &="read", iostat=ierror)
  IF(ierror /= 0) THEN

     print *, inter_fname, "not found"
     STOP
     
  ELSE

     READ(intmface,*) maxinter

     ALLOCATE (inter(maxinter), stat = AllocateStatus)
     IF(AllocateStatus /=0 ) STOP "***Allocation inter not proper***"

     ALLOCATE (domtyp(maxinter), stat = AllocateStatus)
     IF(AllocateStatus /=0 ) STOP "***Allocation domtyp not proper***"
     
     DO i = 1, maxinter

        READ(intmface,*) val, inter(i)

     END DO

     WRITE (logout,*) "Number of interfaces read are ", maxinter
     WRITE (logout,*)  inter

  END IF
 
END SUBROUTINE INTPOSFILE

!--------------------------------------------------------------------

SUBROUTINE RGVALS(iframe)

  USE ANALYZE_PSPEO_WITHIONS
 
  IMPLICIT NONE

  INTEGER :: i,j,molid,atype
  REAL, DIMENSION(1:nchains) :: rgxx, rgyy, rgzz, rgsq
  REAL, DIMENSION(1:nchains) :: rxcm, rycm, rzcm, totmass
  REAL :: rgsqavg, rgxxavg, rgyyavg, rgzzavg
  INTEGER, INTENT(IN) :: iframe

  rgxx = 0.0; rgyy =0.0; rgzz = 0.0; rgsq = 0.0
  totmass = 0.0
  rgsqavg = 0.0; rgxxavg = 0.0; rgyyavg = 0.0; rgzzavg = 0.0
  
  IF(iframe == 1) PRINT *, "Atoms/molecule: ", atperchain
  IF(iframe == 1) PRINT *, "Masses : ", masses

  DO i = 1,npolyatoms

     molid = aidvals(i,2)
     atype = aidvals(i,3)

     IF(molid .GT. nchains) STOP "Counting wrong number of chains"
 
     totmass(molid) = totmass(molid) + masses(atype,1)

     rxcm(molid) = rxcm(molid)+ rxyz_lmp(i,1)*masses(atype,1)
     rycm(molid) = rycm(molid)+ rxyz_lmp(i,2)*masses(atype,1)
     rzcm(molid) = rzcm(molid)+ rxyz_lmp(i,3)*masses(atype,1)

  END DO

  DO i = 1,nchains

     rxcm(i) = rxcm(i)/totmass(i)
     rycm(i) = rycm(i)/totmass(i)
     rzcm(i) = rzcm(i)/totmass(i)

  END DO


  IF(iframe == 1) THEN

     OPEN(unit = 98,file ="totmasschk.txt",action="write",status="repl&
          &ace")

     DO i = 1,nchains

        WRITE(98,'(I0,1X,4(F14.9,1X))') i, totmass(i),rxcm(i),&
             & rycm(i), rzcm(i)

     END DO

     CLOSE(98)

     OPEN(unit = 98,file ="molidchk.txt",action="write",status="repl&
          &ace")

     DO i = 1,npolyatoms

        WRITE(98,'(I0,1X,I0)') i, aidvals(i,2)

     END DO

     CLOSE(98)

  END IF
  

  DO i = 1,npolyatoms

     molid = aidvals(i,2)
     atype = aidvals(i,3)

     rgxx(molid) = rgxx(molid)+masses(atype,1)*((rxyz_lmp(i,1)-rxcm(molid))**2)
     rgyy(molid) = rgyy(molid)+masses(atype,1)*((rxyz_lmp(i,2)-rycm(molid))**2)
     rgzz(molid) = rgzz(molid)+masses(atype,1)*((rxyz_lmp(i,3)-rzcm(molid))**2)

     rgsq(molid) = rgsq(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
          &-rxcm(molid))**2 + (rxyz_lmp(i,2)-rycm(molid))**2 +&
          & (rxyz_lmp(i,3)-rzcm(molid))**2)

  END DO

  DO i = 1,nchains

     rgsq(i) = rgsq(i)/totmass(i)
     rgxx(i) = rgxx(i)/totmass(i)
     rgyy(i) = rgyy(i)/totmass(i)
     rgzz(i) = rgzz(i)/totmass(i)

  END DO

  DO i = 1,nchains

     rgsqavg = rgsqavg + rgsq(i)
     rgxxavg = rgxxavg + rgxx(i)
     rgyyavg = rgyyavg + rgyy(i)
     rgzzavg = rgzzavg + rgzz(i)
     
  END DO
  
  rgsqavg = rgsqavg/REAL(nchains)
  rgxxavg = rgxxavg/REAL(nchains)
  rgyyavg = rgyyavg/REAL(nchains)
  rgzzavg = rgzzavg/REAL(nchains)

!  print *, timestep, sqrt(rgsqavg),rgxxavg,rgyyavg,rgzzavg

  WRITE(rgavgwrite,'(I8,1X,4(F14.6,1X))') timestep, sqrt(rgsqavg),&
       & sqrt(rgxxavg), sqrt(rgyyavg), sqrt(rgzzavg)

  WRITE(rgwrite,'(2(I0,1X))') timestep, nchains

  DO i = 1,nchains
     
     WRITE(rgwrite,'(I4,1X,4(F14.6,1X))') i,rgxx(i),rgyy(i)&
          &,rgzz(i),sqrt(rgsq(i))

  END DO
     

END SUBROUTINE RGVALS

!--------------------------------------------------------------------
!Ref: Borodin and Smith Macromolecules 31, 23, 1998
!To compute g(r) between different ions
!Also computes g(r) within a distance cutoff
!To see whether I can club both?

SUBROUTINE SORTALLARRAYS()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,a1type,cnt,AllocateStatus
  INTEGER, DIMENSION(1:ntotatoms,2) :: dumsortarr

  dumsortarr = -1
  cnt = 0;litype = 0;ptype = 0; otype = 0
  otypearr = -1;litypearr = -1;ptypearr = -1

  DO i = 1,ntotatoms

     a1type = aidvals(i,3)

     IF(a1type == 9 .OR. a1type == 10 .OR. a1type == 11) THEN
        
        cnt = cnt + 1
        dumsortarr(cnt,1) = i
        dumsortarr(cnt,2) = a1type

        IF(a1type == 9)  THEN
           otype  = otype  + 1
           otypearr(otype) = aidvals(i,1)
        ELSEIF(a1type == 10) THEN
           litype = litype + 1
           litypearr(litype) = aidvals(i,1)
           allionids(litype+ptype) = aidvals(i,1)
        ELSEIF(a1type == 11) THEN
           ptype  = ptype  + 1
           ptypearr(ptype) = aidvals(i,1)
           allionids(litype+ptype) = aidvals(i,1)
        END IF

     END IF

  END DO

  PRINT *, "Sorting began.. Number of types.. "
  PRINT *, cnt, otype, litype, ptype
  
  lippair = litype*ptype
  liopair = litype*otype
  popair  = ptype*otype

  IF(cnt .NE. noxyplusions) THEN
     PRINT *, "unequal cnt while sorting"
     PRINT *, cnt, noxyplusions
     STOP
  END IF
  ALLOCATE(sortedarray(cnt,2),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate sortedarray"
  
  i = 0

  DO WHILE(dumsortarr(i+1,1) .NE. -1) 

     i = i + 1
     sortedarray(i,1) = dumsortarr(i,1)
     sortedarray(i,2) = dumsortarr(i,2)

  END DO

  IF(i .NE. cnt) STOP "Wrong total count in sortedarray"

  DO i = 1,cnt

     IF(sortedarray(i,1) == -1 .OR. sortedarray(i,2) == -1) THEN
     
        PRINT *, i,sortedarray(i,1), sortedarray(i,2)
        PRINT *, "Something wrong in assigning sortedarray"
        STOP

     END IF

     IF(sortedarray(i,2) .LT. 9 .OR. sortedarray(i,2) .GT. 11) THEN
     
        PRINT *, i,sortedarray(i,1), sortedarray(i,2)
        PRINT *, "Something wrong in sortedarray type"
        STOP

     END IF

  END DO
        
  rdfpaircnt = cnt

  PRINT *, "Sorting done. Number of pairs.. ",rdfpaircnt
     

END SUBROUTINE SORTALLARRAYS

!--------------------------------------------------------------------

SUBROUTINE ALL_RDF()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id
  REAL :: rxval,ryval,rzval,rval,rboxval
  INTEGER :: LiPcntr,LiOcntr,POcntr
  INTEGER,DIMENSION(0:rmaxbin-1,3) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  rvolavg = rvolavg + rvolval 

  dumrdfarray = 0; LiPcntr = 0; LiOcntr = 0; POcntr = 0

!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval&
!$OMP& ,ibin) REDUCTION(+:dumrdfarray,LiPcntr,LiOcntr,POcntr)
  DO i = 1,rdfpaircnt

     a1id   = sortedarray(i,1)     
     a1type = sortedarray(i,2)

     IF(a1type .LT. 9 .OR. a1type .GT. 11) THEN
        PRINT *, a1type,a1id,i
        STOP "No i type found"
     END IF

     DO j = 1,rdfpaircnt

        a2id   = sortedarray(j,1)        
        a2type = sortedarray(j,2)

        IF(a2type .LT. 9 .OR. a2type .GT. 11) THEN
           PRINT *, a2type,a2id,j
           STOP "No j type found"
        END IF
              
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/rbinval)

        IF(a1type == 10 .AND. a2type == 11) THEN        
                  
           IF(ibin .LT. rmaxbin) THEN
              
              dumrdfarray(ibin,1) = dumrdfarray(ibin,1) + 1
              LiPcntr   = LiPcntr + 1
              
           END IF

        ELSEIF(a1type == 10 .AND. a2type == 9) THEN

           IF(ibin .LT. rmaxbin) THEN
              
              dumrdfarray(ibin,2) = dumrdfarray(ibin,2) + 1
              LiOcntr   = LiOcntr + 1
              
           END IF

        ELSEIF(a1type == 11 .AND. a2type == 9) THEN
           
           IF(ibin .LT. rmaxbin) THEN

              dumrdfarray(ibin,3) = dumrdfarray(ibin,3) + 1
              POcntr   = POcntr + 1
              
           END IF

        END IF
           
     END DO

  END DO
!$OMP END DO

!$OMP DO
  DO i = 0,rmaxbin-1

     rdfarray(i,1) = rdfarray(i,1) + REAL(dumrdfarray(i,1))&
          &*rvolval/(REAL(lippair))
     rdfarray(i,2) = rdfarray(i,2) + REAL(dumrdfarray(i,2))&
          &*rvolval/(REAL(liopair))
     rdfarray(i,3) = rdfarray(i,3) + REAL(dumrdfarray(i,3))&
          &*rvolval/(REAL(popair))
  END DO
!$OMP END DO

!$OMP END PARALLEL


!  PRINT *, "cntr", LiPcntr, LiOcntr, POcntr

END SUBROUTINE ALL_RDF

!--------------------------------------------------------------------

SUBROUTINE SORTOXYFREECOMPLEX()

  USE OMP_LIB
  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,boundflag,ocnt,licnt,tid,cnt
  INTEGER :: ofree,obound,oboundtot,ofreetot,AllocateStatus
  REAL :: rxval,ryval,rzval,rval
  INTEGER, DIMENSION(1:n_o_mons,0:nproc-1) :: oxybounddum,oxyfreedum

 ofree = 0; obound = 0; oboundtot=0; ofreetot=0;licnt = 1

!!$!$OMP PARALLEL
!!$!$OMP DO PRIVATE(i,j)
  DO i = 1,n_o_mons
     
     DO j = 0,nproc-1
        
        oxybounddum(i,j) = -1
        oxyfreedum(i,j)  = -1

     END DO

  END DO
!!$!$OMP END DO

!!$!$OMP DO PRIVATE(ocnt,a1id,boundflag,a2id,rxval,ryval,rzval,rval,tid)& 
!!$!$OMP& REDUCTION(+:oboundtot,ofreetot) FIRSTPRIVATE(ofree,obound,licnt)
  DO ocnt = 1,n_o_mons

     tid  = OMP_GET_THREAD_NUM()
     a1id = otypearr(ocnt)
     licnt = 1
     boundflag = 0

     DO WHILE(licnt .LE. nliions)
     
        a2id  = litypearr(licnt)

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        IF(rval .LT. rclus_cut) THEN

           obound = obound + 1
           oboundtot = oboundtot + 1
           oxybounddum(obound,tid) = a1id
           licnt = nliions + 1
           boundflag = 1

        ELSE
           
           licnt = licnt + 1

        END IF

     END DO
           
     IF(boundflag == 0) THEN

        ofree = ofree + 1
        ofreetot = ofreetot + 1
        oxyfreedum(ofree,tid) = a1id

     END IF

  END DO
!!$!$OMP END DO
!!$!$OMP END PARALLEL

  IF(oboundtot == 0 .OR. ofreetot == 0) THEN

     PRINT *, "Zero oboundtot/ofreetot"
     PRINT *, oboundtot, ofreetot
     STOP

  END IF

  ALLOCATE(oxyboundarr(1:oboundtot),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate oxyboundarr"

  ALLOCATE(oxyfreearr(1:ofreetot),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate oxyfreearr"
  
  ofree = 0; obound = 0
 
  DO i = 0,nproc-1
     
     cnt = 0

     DO WHILE(oxyfreedum(cnt+1,i) .NE. -1)

        cnt = cnt + 1
        ofree = ofree + 1
        oxyfreearr(ofree) = oxyfreedum(cnt,i)

     END DO
 
     cnt = 0

     DO WHILE(oxybounddum(cnt+1,i) .NE. -1)
        
        cnt = cnt + 1
        obound = obound + 1
        oxyboundarr(obound) = oxybounddum(cnt,i)

     END DO

  END DO

  IF(ofree .NE. ofreetot .OR. obound .NE. oboundtot) THEN

     PRINT *, "Unequal assignment in free and bound oxygens"
     PRINT *, ofree,ofreetot,obound,oboundtot
     STOP

  END IF

  CALL FREE_BOUND_OXYRDF(ofreetot,oboundtot)

  DEALLOCATE(oxyfreearr)
  DEALLOCATE(oxyboundarr)

END SUBROUTINE SORTOXYFREECOMPLEX

!--------------------------------------------------------------------

SUBROUTINE FREE_BOUND_OXYRDF(ofreetot,oboundtot)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,otot,a1id,a2id,ibin
  INTEGER :: obound_free, obound_bound, ofree_free
  REAL    :: rval,rxval,ryval,rzval,rboxval
  INTEGER,INTENT(IN) :: ofreetot,oboundtot
  INTEGER,DIMENSION(0:rmaxbin-1) ::dumrdf_o_ff,dumrdf_o_bb,dumrdf_o_fb

  rvolval = box_xl*box_yl*box_zl
  IF(rdfcalc == .false.) THEN !Already in RDF. So if it is true it is
     !already accounted in rdf computation
  
     rvolavg = rvolavg + rvolval
  
  END IF

  otot = ofreetot + oboundtot

  obound_free = 0; obound_bound = 0; ofree_free = 0

!$OMP PARALLEL 

!$OMP DO PRIVATE(i)
  DO i = 0,rmaxbin-1
     
     dumrdf_o_fb(i) = 0
     dumrdf_o_ff(i) = 0
     dumrdf_o_bb(i) = 0
     
  END DO
!$OMP END DO

!$OMP DO REDUCTION(+:dumrdf_o_fb,obound_free,ofree_free,obound_bound)&
!$OMP& PRIVATE(i,j,a1id,a2id,rval,rxval,ryval,rzval,ibin)
  DO i = 1,ofreetot

     a1id = oxyfreearr(i)
     
     IF(aidvals(a1id,3) .NE. 9) STOP "Wrong Oxygen atom for a1id"

     DO j = i,oboundtot

        a2id   = oxyboundarr(j)

        IF(aidvals(a2id,3) .NE. 9) STOP "Wrong Oxygen atom for a2id"

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/roxybinval)

        IF(ibin .LT. rmaxbin) THEN

           dumrdf_o_fb(ibin) = dumrdf_o_fb(ibin) + 2
           obound_free = obound_free  + 1

        END IF

     END DO

  END DO
!$OMP END DO

!$OMP DO REDUCTION(+:dumrdf_o_ff,obound_free,ofree_free,obound_bound)&
!$OMP& PRIVATE(i,j,a1id,a2id,rval,rxval,ryval,rzval,ibin)
  DO i = 1,ofreetot

     a1id = oxyfreearr(i)
     
     IF(aidvals(a1id,3) .NE. 9) STOP "Wrong Oxygen atom for a1id"

     DO j = i+1,ofreetot

        a2id   = oxyfreearr(j)

        IF(aidvals(a2id,3) .NE. 9) STOP "Wrong Oxygen atom for a2id"

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/roxybinval)
        
        IF(ibin .LT. rmaxbin) THEN
           
           dumrdf_o_ff(ibin) = dumrdf_o_ff(ibin) + 2
           ofree_free = ofree_free + 1
           
        END IF

     END DO

  END DO
!$OMP END DO

!$OMP DO REDUCTION(+:dumrdf_o_bb,obound_free,ofree_free,obound_bound) &
!$OMP& PRIVATE(i,j,a1id,a2id,rval,rxval,ryval,rzval,ibin)
  DO i = 1,oboundtot

     a1id = oxyboundarr(i)
     
     IF(aidvals(a1id,3) .NE. 9) STOP "Wrong Oxygen atom for a1id"

     DO j = i+1,oboundtot

        a2id   = oxyboundarr(j)

        IF(aidvals(a2id,3) .NE. 9) STOP "Wrong Oxygen atom for a2id"

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/roxybinval)
        
        IF(ibin .LT. rmaxbin) THEN
           
           dumrdf_o_bb(ibin) = dumrdf_o_bb(ibin) + 2
           obound_bound = obound_bound + 1
           
        END IF

     END DO

  END DO
!$OMP END DO


!$OMP DO

  DO i = 0,rmaxbin-1

     rdf_o_fb(i) = rdf_o_fb(i) + REAL(dumrdf_o_fb(i))*rvolval&
          &/(REAL(oboundtot)*REAL(ofreetot))
     rdf_o_ff(i) = rdf_o_ff(i) + REAL(dumrdf_o_ff(i))*rvolval&
          &/(REAL(ofreetot)*REAL(ofreetot))
     rdf_o_bb(i) = rdf_o_bb(i) + REAL(dumrdf_o_bb(i))*rvolval&
          &/(REAL(oboundtot)*REAL(oboundtot))

  END DO
!$OMP END DO
!$OMP END PARALLEL



END SUBROUTINE FREE_BOUND_OXYRDF

!--------------------------------------------------------------------

SUBROUTINE COARSEGRAIN_PF6()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,molid,molinit, molcnt
  INTEGER, DIMENSION(nPF6ions,3) :: flagp

  flagp = -1

  DO i = 1,ntotatoms
     
     IF(aidvals(i,3) == 11) THEN !First PF6 ion

        molinit = aidvals(i,2)
        EXIT

     END IF
       
  END DO
  
  PRINT *, "First PF6 molid : ", molinit
  WRITE(logout,*) "First PF6 molid : ", molinit

  DO i = 1,ntotatoms

     IF(aidvals(i,3) .GE. 11) THEN

        molid  = aidvals(i,2)
        molcnt = molid - molinit + 1

        IF(aidvals(i,3) == 11) THEN

           PF6aids(molcnt,1) = aidvals(i,1)

        ELSEIF(aidvals(i,3) == 12) THEN
!Sorting F atoms into rows corresponding to each molecule
           IF(flagp(molcnt,1) == -1) THEN

              PF6aids(molcnt,2) = aidvals(i,1) !First F
              flagp(molcnt,1) = 1

           ELSEIF(flagp(molcnt,1) == 1) THEN
              
              PF6aids(molcnt,5) = aidvals(i,1) !Second F
              
           END IF

        ELSEIF(aidvals(i,3) == 13) THEN

           IF(flagp(molcnt,2) == -1) THEN

              PF6aids(molcnt,3) = aidvals(i,1)
              flagp(molcnt,2) = 1

           ELSEIF(flagp(molcnt,2) == 1) THEN
              
              PF6aids(molcnt,6) = aidvals(i,1)
              
           END IF
              

        ELSEIF(aidvals(i,3) == 14) THEN

           IF(flagp(molcnt,3) == -1) THEN

              PF6aids(molcnt,4) = aidvals(i,1)
              flagp(molcnt,3) = 1

           ELSEIF(flagp(molcnt,3) == 1) THEN
              
              PF6aids(molcnt,7) = aidvals(i,1)
              
           END IF
              
        ELSE 

           PRINT *, "Unknown atype :", aidvals(i,3),molid
           STOP

        END IF

     END IF

  END DO

  OPEN(unit = 90,file="pf6neigh.txt",action="write",status="replace")

  DO i = 1,nPF6ions

     WRITE(90,'(8(I0,1X))') i,PF6aids(i,:)

  END DO

  CLOSE(90)

END SUBROUTINE COARSEGRAIN_PF6

!--------------------------------------------------------------------

SUBROUTINE LiP_ION_CLUSTERS()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,clust_cnt,tid
  INTEGER,DIMENSION(1:maxclustsize,0:nproc-1) :: lipclust_inst&
       &,p_liclust_inst 
  REAL :: rxval, ryval, rzval, rval

  lipclust_inst = 0; p_liclust_inst = 0

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,clust_cnt,tid)
!$OMP DO
  DO i = 1,nLiions
     
     clust_cnt = 0
     a1id = litypearr(i)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,nPF6ions

        a2id = ptypearr(j)
        
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rclus_cut) THEN
           
           clust_cnt = clust_cnt + 1
           
        END IF

     END DO

     IF(clust_cnt + 1 .GT. maxclustsize) THEN

        PRINT *, "Cluster count exceeded maxclustsize"
        PRINT *, clust_cnt, maxclustsize
        STOP

     END IF

     lipclust_inst(clust_cnt+1,tid) = lipclust_inst(clust_cnt+1,tid) &
          &+ 1

  END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,clust_cnt,tid)
!$OMP DO
  DO i = 1,nPF6ions
     
     clust_cnt = 0
     a1id = ptypearr(i)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,nLiions

        a2id = litypearr(j)
        
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rclus_cut) THEN
           
           clust_cnt = clust_cnt + 1
           
        END IF

     END DO

     IF(clust_cnt + 1 .GT. maxclustsize) THEN

        PRINT *, "Cluster count exceeded maxclustsize"
        PRINT *, clust_cnt, maxclustsize
        STOP

     END IF

     p_liclust_inst(clust_cnt+1,tid) = p_liclust_inst(clust_cnt+1&
          &,tid) + 1

  END DO
!$OMP END DO


!$OMP DO 
  DO  i = 1,maxclustsize
     DO j = 0,nproc-1
        lipclustavg(i) = lipclustavg(i) + lipclust_inst(i,j)
        p_liclustavg(i) = p_liclustavg(i) + p_liclust_inst(i,j)
     END DO
  END DO
!$OMP END DO


!$OMP END PARALLEL


END SUBROUTINE LIP_ION_CLUSTERS

!--------------------------------------------------------------------

SUBROUTINE LiP_ION_CLUSTERS_PEO_DOMAIN(tval)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,clust_cnt,tid
  INTEGER,DIMENSION(1:maxclustsize,0:nproc-1) :: lipclust_inst&
       &,p_liclust_inst 
  REAL :: rxval, ryval, rzval, rval
  INTEGER,INTENT(IN) :: tval

  lipclust_inst = 0; p_liclust_inst = 0


!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,clust_cnt,tid)
!$OMP DO
  DO i = 1,nLiions
     
     clust_cnt = 0
     a1id = litypearr(i)
     tid = OMP_GET_THREAD_NUM()

     IF(seg_dtype(a1id) == 2) THEN !For PEO based

        DO j = 1,nPF6ions

           a2id = ptypearr(j)
           
           rxval = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
           ryval = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
           rzval = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LT. rclus_cut) THEN
              
              clust_cnt = clust_cnt + 1
              
           END IF
           
        END DO
        
        IF(clust_cnt + 1 .GT. maxclustsize) THEN
           
           PRINT *, "Cluster count exceeded maxclustsize"
           PRINT *, clust_cnt, maxclustsize
           STOP
           
        END IF
        
        lipclust_inst(clust_cnt+1,tid) = lipclust_inst(clust_cnt+1&
             &,tid) + 1

     END IF

  END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,clust_cnt)
!$OMP DO
  DO i = 1,nPF6ions
     
     clust_cnt = 0
     a1id = ptypearr(i)
     tid = OMP_GET_THREAD_NUM()

     IF(seg_dtype(a1id) == 2) THEN !For PEO

        DO j = 1,nLiions

           a2id = litypearr(j)
        
           rxval = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
           ryval = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
           rzval = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LT. rclus_cut) THEN
              
              clust_cnt = clust_cnt + 1
              
           END IF
           
        END DO

        IF(clust_cnt + 1 .GT. maxclustsize) THEN
        
           PRINT *, "Cluster count exceeded maxclustsize"
           PRINT *, clust_cnt, maxclustsize
           STOP
           
        END IF

        p_liclust_inst(clust_cnt+1,tid) = p_liclust_inst(clust_cnt+1&
             &,tid) + 1

     END IF
     
  END DO
!$OMP END DO
  
!$OMP DO 

  DO  i = 1,maxclustsize
     DO j = 0,nproc-1
        lipclustavg_peo(i) = lipclustavg_peo(i) + lipclust_inst(i,j)
        p_liclustavg_peo(i) = p_liclustavg_peo(i) + p_liclust_inst(i&
             &,j)
     END DO
  END DO
!$OMP END DO


!$OMP END PARALLEL

END SUBROUTINE LIP_ION_CLUSTERS_PEO_DOMAIN

!--------------------------------------------------------------------

SUBROUTINE CLUSTER_ANALYSIS(frnum)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

!Ref Sevick et.al ., J Chem Phys 88 (2)

  INTEGER :: i,j,k,a2ptr,a1id,a2id,itype,jtype,jptr,idum,jflag,jcnt&
       &,iflag,jtot,jind,jprev
  INTEGER, DIMENSION(ntotion_centers,ntotion_centers) :: all_direct&
       &,LiP_direct,all_neigh,LiP_neigh
  INTEGER, DIMENSION(1:ntotion_centers) :: union_all,flag_LiP,scnt&
       &,all_linked
  REAL :: rxval, ryval, rzval, rval
  INTEGER, INTENT(IN) :: frnum

!$OMP PARALLEL SHARED(LiP_direct,LiP_neigh)

!$OMP DO PRIVATE(i,j)
  DO i = 1,ntotion_centers

     scnt(i) = 0; all_linked(i)  = 0
     union_all(i) = -1; flag_LiP(i) = -1

     DO j = 1,ntotion_centers

        IF(i .NE. j) THEN
           all_direct(i,j) = 0
           LiP_direct(i,j) = 0
        END IF

        IF(i == j) THEN
           all_direct(i,j) = 1
           LiP_direct(i,j) = 0
        END IF

        all_neigh(i,j) = 0
!        LiP_neigh(i,j) = 0

     END DO

  END DO
!$OMP END DO

!Create Direct connectivity matrix
!allneigh - does not distinguish between Li and P neigh
!LiP_neigh - neighbors with sequence Li-P-Li-P.. or P-Li-P-Li...
!$OMP DO PRIVATE(i,j,a1id,a2ptr,a2id,rxval,ryval,rzval,rval,itype&
!$OMP& ,jptr,jtype)  
  DO i = 1,ntotion_centers
     
     a1id = allionids(i)
     a2ptr = 1
     itype = aidvals(a1id,3)
     jptr  = 1
     all_neigh(i,i) = a1id
     
     DO j = 1,ntotion_centers
        
        a2id = allionids(j)
        
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)

        IF(rval .LT. rclus_rdf_cut .AND. a1id .NE. a2id) THEN
           
           all_direct(i,j) = 1
           all_neigh(i,j)  = a2id
!           all_neigh(i,a2ptr) = a2id
!           a2ptr = a2ptr + 1
           
           jtype = aidvals(a2id,3)
           
           IF(itype .NE. jtype) THEN
              
              LiP_direct(i,j) = 1
!              LiP_neigh(i,j)  = a2id
!              LiP_neigh(i,jptr+1) = a2id
              itype = jtype
!              jptr  = jptr + 1
              
           END IF
           
        END IF

     END DO
     
  END DO
!$OMP END DO  

!Check for symmetry
  IF(frnum == 1) THEN
!$OMP DO
     DO i = 1,ntotion_centers

        DO j = 1,ntotion_centers

           IF(all_direct(i,j) .NE. all_direct(j,i)) STOP "Unsymmetric&
                & all_direct"

          IF(all_neigh(i,j) .NE. 0) THEN

              IF(all_neigh(i,j) .NE. all_neigh(j,j) .OR. all_neigh(j&
                   &,i) .NE. all_neigh(i,i)) THEN

                 PRINT *, i,j,all_direct(i,j),all_direct(j,i)&
                      &,all_neigh(j,i),all_neigh(i,i)
                 STOP "Unsymmetric neighbor list"
                 
              END IF
              
           END IF

        END DO

     END DO
!$OMP END DO
  END IF

!$OMP END PARALLEL        


!Intersection

  DO i = 1,ntotion_centers-1 !Ref row counter

     iflag = 0
     idum  = i

     DO WHILE(iflag == 0 .AND. union_all(i) == -1)
     
        jflag = 0
        k    = 1 !Column counter
        j    = idum+1 !Other row counter
        
        DO WHILE(jflag == 0 .AND. k .LE. ntotion_centers) 
           
           IF((all_direct(i,k) == all_direct(j,k)).AND. all_direct(i&
                &,k)== 1) THEN
              
              jflag = 1
!!$              jprev = 0

              DO jcnt = 1,ntotion_centers
               
                 
!!$                 IF(all_direct(j,jcnt) == 1) jprev = 1

                 !Replace highest row by union of two rows

                 all_direct(j,jcnt) = all_direct(i,jcnt) .OR.&
                      & all_direct(j,jcnt)

!!$                 IF((all_direct(j,jcnt) == 1 .AND. all_direct(i,jcnt)&
!!$                      &==1) .AND. jprev == 0) THEN
!!$                    
!!$                    all_neigh(j,jcnt) = all_neigh(i,jcnt)
!!$                    jprev = 0 !Other condition is already
!!$                    ! incorporated before
!!$                 END IF
!!$                 
              END DO
              
              union_all(i) = 1 !One match implies the low ranked row
              ! is present in high ranked row
              
           ELSE
              
              k = k + 1
              
           END IF
           
        END DO
        
        IF(union_all(i) == 1) THEN
           
           iflag = 1
           
        ELSE
           
           idum  = idum + 1
           
        END IF
        
        IF(idum == ntotion_centers) iflag = 1

     END DO

  END DO
  
!Count
  jtot = 0
!$OMP PARALLEL PRIVATE(i,j,jind) 
!$OMP DO
  DO i = 1,ntotion_centers

     IF(union_all(i) == -1) THEN
        
        jind = 0

        DO j = 1,ntotion_centers

           IF(all_direct(i,j) == 1) jind = jind + 1

        END DO

        scnt(jind) = scnt(jind) + 1
        all_linked(i) = jind

     END IF
     
  END DO
!$OMP END DO

!$OMP DO

  DO i = 1,ntotion_centers

     clust_avg(i) = clust_avg(i) + scnt(i)

  END DO
!$OMP END DO

!$OMP END PARALLEL

  IF(frnum == 1) THEN
     OPEN(unit =90,file ="scnt.txt",action="write",status="replace")  
  END IF

  jtot = 0
  
  DO i = 1,ntotion_centers
     
     IF(frnum == 1) WRITE(90,*) i,scnt(i)
     jtot = jtot + all_linked(i)
     
  END DO
  
  IF(jtot .NE. ntotion_centers) THEN
     
     PRINT *, "Sum of ions not equal to total ion centers"
     PRINT *, jtot, ntotion_centers
     STOP
     
  END IF

  IF(frnum == 1) CLOSE(90)
  
  IF(frnum == 1) THEN

     OPEN(unit =90,file ="all_neigh.txt",action="write",status="replace")  
    
     DO i = 1,ntotion_centers
        
        IF(union_all(i) == -1) THEN
           
           WRITE(90,*) i,all_linked(i)
        
           DO j = 1,ntotion_centers
           
              IF(all_direct(i,j) == 1) WRITE(90,*) j,allionids(j),&
                   & all_direct(i,j)
           
           END DO

        END IF
        
     END DO

     CLOSE(90)

  END IF
  
END SUBROUTINE CLUSTER_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE COMPARTMENTALIZE_PARTICLES(tval)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,k,p,allocatestatus,c1,c2,c3,flagz,flagz2,flagbin,a1id
  REAL    :: par_pos, segeps, domcheck, eps2,eps1
  REAL    :: interin, interout, interin2, interout2
  INTEGER, INTENT(IN) :: tval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: dum_seg, dum_typ

  domcheck = 0.0; epsinit = 0
  DO i = 1, maxinter-1
     domcheck = domcheck + (inter(i+1)-inter(i))
  END DO
  !One can use any of the two following definitions. One will be a
  !lesser average for domavg
  !Method 1
  domcheck = domcheck + boxval - inter(maxinter) + inter(1)
  domavg   = domcheck/REAL(maxinter)
  !Method 2
!!$  domavg   = domcheck/REAL(maxinter-1)

  segeps   = segper*domavg !Divide half wavelength

  IF(tval == 1) THEN
     WRITE(logout,*) "Domain width is ", domavg
     WRITE(logout,*) "Basic segmental width is ", segeps
  END IF
   
  ALLOCATE (dum_seg(ntotatoms), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "*** Allocation dum_seg not proper ***"
  ALLOCATE (dum_typ(ntotatoms), stat = AllocateStatus)
  IF(AllocateStatus /=0 ) STOP "*** Allocation dum_typ not proper ***"

  IF(int_possort == .true.) THEN
     ALLOCATE (seg_dtype(ntotatoms), stat = AllocateStatus)
     IF(AllocateStatus /=0 ) STOP "*** Allocation seg_dtype not proper&
          & ***"
     CALL DOMVAL(tval)
     IF(sort_ana) CALL SORT_WRT_DOM(tval)
  ELSE
     ALLOCATE (seg_dtype(ntotatoms), stat = AllocateStatus)
     DEALLOCATE(seg_dtype)
  END IF
  
  IF(tval == 1) THEN
     WRITE(logout,*) "Increment in divisions: ", epsinc*segeps
     WRITE(logout,*) "Width for each segment: ", epspre*epsinc*segeps
  END IF

  DO p = 1,pmax

     dum_seg = -1 !initialize everything to -1
     dum_typ = -1 !initialize everything to -1
     
     eps1 = epsinit !takes monomers between eps1*segeps and eps2*segeps
     eps2 = eps1 + epspre*REAL(epsinc) !prefactor is a fraction of epsinc
 
     IF(eps2*segeps > 0.45*domavg) THEN !difference of 0.1 to account
        !for imperfections in the domain width
        
        WRITE(logout,*) "Domain between ", eps1*segeps, "and", eps2&
             &*segeps, "is greater than half domain size"
        WRITE(logout,*) "Domain more than half domain size"
        PRINT*,  "Domain width more than half domain size at", tval
        EXIT

     END IF

     IF(tval == 1) THEN
        WRITE(epsnum,'(I3)') int(epsinit)
        dum_fname  = 'segcoords.'//trim(adjustl(fwnum))&
             &//"_"//trim(adjustl(epsnum))//".txt"

        OPEN(unit = 14,file = dum_fname, status="replace", action &
             &="write")

     END IF

     c1 = 0;c2 = 0;c3 = 0
     k = 0
     IF(tval == 1) print *, boxval
     DO i = 1, ntotion_centers+n_o_mons

        a1id = sortedarray(i,1)
!        IF(i .LE. ntotion_centers) a1id = allionids(i)
!        IF(i .GT. ntotion_centers) a1id = otypearr(i-ntotion_centers)
        
        IF(major_axis == 1) THEN
           par_pos = trx_lmp(a1id,tval) - boxval*floor(trx_lmp(a1id&
                &,tval)/boxval)
        ELSE IF(major_axis == 2) THEN 
           par_pos = try_lmp(a1id,tval) - boxval*floor(try_lmp(a1id&
                &,tval)/boxval)
        ELSE
           par_pos = trz_lmp(a1id,tval) - boxval*floor(trz_lmp(a1id&
                &,tval)/boxval)
        END IF

        j = 1
        flagbin = 0
!Both interfaces are considered that is to the left and the right.
!Can use ANINT conditions like for PBC?? Need to check.
!Changing order of DO loops can result in efficiency, but memory to
!store arrays for each interface has to be created 
       DO WHILE(j <= maxinter)
         
          !To right of interface
          IF((inter(j) + eps1*segeps) > boxval) THEN
             !if inner interface is outside box, outer interface will
             !also be. Flag will not change as this is equivalent to
             !all the particles and interfaces inside a "new" box
             !with the interface location at the other end of the box.
             interin  = inter(j) + eps1*segeps - boxval 
             interout = inter(j) + eps2*segeps - boxval
             flagz    = 0 !Keep old flag
             
          ELSE IF((inter(j) + eps2*segeps) > boxval) THEN
             !If inner interface is inside box and outer interface is
             !outside box. Only in this case, one needs to
             !redefine the conditions for selecting a particle.
             interin  = inter(j) + eps1*segeps
             interout = inter(j) + eps2*segeps - boxval
             flagz    = 1 !Update flag
             
          ELSE
             !Other normal conditions
             interin  = inter(j) + eps1*segeps
             interout = inter(j) + eps2*segeps
             flagz    = 0 !Keep old flag
             
          END IF
          
          !To left of interface
          
          IF((inter(j) - eps1*segeps) < 0.0) THEN
             !if inner interface is outside box, outer interface will
             !also be. Flag will not change as this is equivalent to
             !all the particles and interfaces inside a "new" box
             !with the interface location at the other end of the box.
             interin2  = inter(j) - eps1*segeps + boxval 
             interout2 = inter(j) - eps2*segeps + boxval
             flagz2    = 0 !Keep old flag
             
          ELSE IF((inter(j) - eps2*segeps) < 0.0) THEN
             !If inner interface is inside box and outer interface is
             !outside box. Only in this case, one needs to
             !redefine the conditions for selecting a particle.
             interin2  = inter(j) - eps1*segeps
             interout2 = inter(j) - eps2*segeps + boxval
             flagz2    = 1 !Update flag
             
          ELSE
             !Other normal conditions
             interin2  = inter(j) - eps1*segeps
             interout2 = inter(j) - eps2*segeps
             flagz2    = 0 !Keep old flag
              
          END IF
          
          !Right of interface      
          IF(flagz == 0) THEN !Conditions 1 and 3 above
             
             IF(par_pos .GE. interin .AND. par_pos .LT. interout)&
                  & THEN
              
                IF(aidvals(a1id,3) == 10) THEN 

                   c1 = c1 + 1
                   
                ELSEIF(aidvals(a1id,3) == 11) THEN

                   c2 = c2 + 1

                ELSEIF(aidvals(a1id,3) == 9) THEN
                   
                   c3 = c3 + 1
                   
                END IF
                
                flagbin = flagbin + 1
                k = k + 1
                dum_seg(k) = a1id
                dum_typ(k) = aidvals(a1id,3)
                j = maxinter
                IF(tval == 1) WRITE(14,'(2(I0,1X),3(F14.8,1X))') i,&
                     & aidvals(a1id,3),trx_lmp(a1id,tval)&
                     &,try_lmp(a1id,tval),trz_lmp(a1id,tval)
                
             END IF

          ELSEIF(flagz == 1) THEN!Condition 2 in the previous IF-ELSE LOOP
             !Change conditions for inner interface inside and outer
             !interface outside the box. AND changes to OR
             
             !Right of interface
             IF(par_pos .GE. interin .OR. par_pos .LT. interout) THEN
                
                IF(aidvals(a1id,3) == 10) THEN 
                   
                   c1 = c1 + 1
                   
                ELSEIF(aidvals(a1id,3) == 11) THEN

                   c2 = c2 + 1

                ELSEIF(aidvals(a1id,3) == 9) THEN

                   c3 = c3 + 1
                   
                END IF
                
                k = k + 1
                dum_seg(k) = a1id
                dum_typ(k) = aidvals(a1id,3)
                j = maxinter
                flagbin = flagbin + 1
                IF(tval == 1) WRITE(14,'(2(I0,1X),3(F14.8,1X))') i,&
                     & aidvals(a1id,3),trx_lmp(a1id,tval)&
                     &,try_lmp(a1id,tval),trz_lmp(a1id,tval)
                
             END IF
             
          END IF
          
          !Left of interface. Interfacial conditions remain identical
          IF(flagz2 == 0) THEN
             
             IF(par_pos .GE. interout2 .AND. par_pos .LT. interin2)&
                  & THEN
                
                IF(aidvals(a1id,3) == 10) THEN 
                   c1 = c1 + 1
                    
                ELSEIF(aidvals(a1id,3) == 11) THEN 

                   c2 = c2 + 1

                ELSEIF(aidvals(a1id,3) == 9) THEN

                   c3 = c3 + 1

                END IF
                
                k = k + 1
                dum_seg(k) = a1id
                dum_typ(k) = aidvals(a1id,3)
                j = maxinter
                flagbin = flagbin + 1

                IF(tval == 1) WRITE(14,'(2(I0,1X),3(F14.8,1X))') i,&
                     & aidvals(a1id,3),trx_lmp(a1id,tval)&
                     &,try_lmp(a1id,tval),trz_lmp(a1id,tval)

                
             END IF
             
          ELSEIF(flagz2 == 1) THEN
             
             IF(par_pos .GE. interout2 .OR. par_pos .LT. interin2)&
                  & THEN
                
                IF(aidvals(a1id,3) == 10) THEN 
                   
                   c1 = c1 + 1

                ELSEIF(aidvals(a1id,3) == 11) THEN 

                   c2 = c2 + 1

                ELSEIF(aidvals(a1id,3) == 9) THEN

                   c3 = c3 + 1

                END IF
                
                k = k + 1
                dum_seg(k) = a1id
                dum_typ(k) = aidvals(a1id,3)
                j = maxinter
                flagbin = flagbin + 1
                IF(tval == 1) WRITE(14,'(2(I0,1X),3(F14.8,1X))') i,&
                     & aidvals(a1id,3),trx_lmp(a1id,tval)&
                     &,try_lmp(a1id,tval),trz_lmp(a1id,tval)

             END IF
             
          END IF
          
          j = j + 1

          IF(flagbin > 1) THEN
             
             WRITE(logout,*) "Particle binned twice"
             WRITE(logout,*) interin, interout, interin2, interout2
             WRITE(logout,'(I0,1X,I0,1X,3F14.8)') i, aidvals(a1id,3),&
                  & trx_lmp(a1id,tval),try_lmp(a1id,tval)&
                  &,trz_lmp(a1id,tval)
             STOP
             
          END IF
          
       END DO
       
    END DO
    
    IF(tval == 1)WRITE(logout,*) "Width considered is ", (eps2-eps1)&
         &*segeps
    
    IF(tval == 1) PRINT *,"The number of mon for ",eps1*segeps,"< z &
         &<",eps2*segeps,"is :", k, "with ", c1, "type Li", c2, "type &
         &P and ", c3, "type O at tval: ",tval
    IF(tval == 1) THEN
       WRITE(logout,*),"The number of mon for ",eps1*segeps, "< z <",&
            & eps2*segeps,"is:",k,"with ", c1, "type Li", c2, "type P&
            & and ", c3, "type O at tval: ",tval
    END IF
    
    
    num_mons = k
    
    ALLOCATE (seg_act(num_mons), stat = AllocateStatus)
    IF(AllocateStatus /=0 ) STOP "*** Allocation seg_act not proper ***"
    
    ALLOCATE (seg_typ(num_mons), stat = AllocateStatus)
    IF(AllocateStatus /=0 ) STOP "*** Allocation seg_typ not proper ***"
    
    DO i = 1,num_mons
       
       IF(dum_seg(i) == -1 .OR. dum_typ(i) == -1) THEN
          
          PRINT *, "dummy arrays for seg_mons updated incorrectly"
          
          STOP
          
       END IF
       
       seg_act(i) = dum_seg(i)
       seg_typ(i) = dum_typ(i)
        
    END DO

    CALL LAMPLANE_ANALYSIS(tval,c1,c2,c3,p,(eps2-eps1)*segeps)
    
    DEALLOCATE(seg_act)
    DEALLOCATE(seg_typ)
    
    epsinit = epsinit + epsinc
    
 END DO
 
 
 pfin = p - 1
 IF(tval == 1) PRINT *, "Total number of divisions: ", pfin
 IF(pmax .LT. pfin) PRINT *, "Warning: Not all divisions are accounted&
      & for: pmax/pfin", pmax, pfin
 DEALLOCATE(dum_seg)
 DEALLOCATE(dum_typ)
 DEALLOCATE(seg_dtype)
 
 WRITE(logout,*) "Segmental diffusion analysis complete .. "
 WRITE(logout,*) "Segmental mode relaxation analysis complete .. "
 
END SUBROUTINE COMPARTMENTALIZE_PARTICLES

!--------------------------------------------------------------------

SUBROUTINE LAMPLANE_ANALYSIS(tval,c1,c2,c3,ipos,widbin)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  REAL, INTENT(IN) :: widbin
  INTEGER, INTENT(IN) :: tval,c1,c2,c3,ipos
  INTEGER :: t1,t2,clock_rate,clock_max

  IF(rdf_2d .AND. tval == 1) THEN
     IF(ipos == 1) THEN
        rdf2darray = 0.0; r2dvolavg = 0.0
     END IF
     CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
     CALL RDF2D_LAMELLAR_PLANE(tval,c1,c2,c3,ipos,widbin)
     CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
     PRINT *, 'Elapsed real time for 2D RDF= ',REAL(t2-t1)/&
          & REAL(clock_rate)           
  ELSEIF(rdf_2d .AND. mod(tval,rdffreq) == 0) THEN
     CALL RDF2D_LAMELLAR_PLANE(tval,c1,c2,c3,ipos,widbin)
  END IF

  IF(rdf_2ddom .AND. tval == 1) THEN
     IF(ipos == 1) THEN
        rdf2deo = 0.0; rdf2dps = 0.0
     END IF
     CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
     CALL RDF2D_LAMELLAR_PLANE_WRTDOM(tval,c1,c2,c3,ipos,widbin)
     CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
     PRINT *, 'Elapsed real time for 2D RDF= ',REAL(t2-t1)/&
          & REAL(clock_rate)           
  ELSEIF(rdf_2ddom .AND. mod(tval,rdffreq) == 0) THEN
     CALL RDF2D_LAMELLAR_PLANE_WRTDOM(tval,c1,c2,c3,ipos,widbin)
  END IF
     
!    IF(int_sgdif) CALL SEGDIF(c1, c2)
!    CALL MIDSEGDIF()

!    IF(int_vanhv) CALL VANHOVE(c1, c2)
!    IF(int_qrelx) CALL SEGMODRELAX(c1,c2)
!    IF(int_segsc) CALL SEGSELFCONC(c1,c2)


END SUBROUTINE LAMPLANE_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE DOMVAL(tval)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  REAL :: par_pos
  INTEGER :: flagdom, domnum,i,j,csolvA, csolvB
  INTEGER, INTENT(IN) :: tval

  seg_dtype = 0
  IF(tval == 1) THEN
     dum_fname  = 'domdist.'//trim(adjustl(fwnum))//".txt"
     OPEN(unit = 18,file = dum_fname, status="replace", action = "writ&
          &e")
  END IF

!$OMP PARALLEL DO PRIVATE(i,par_pos,flagdom,domnum)
  DO i = 1, ntotatoms
     
     IF(major_axis == 1) THEN
        par_pos = trx_lmp(i,tval) - boxval*floor(trx_lmp(i,tval)&
             &/boxval)
     ELSE IF(major_axis == 2) THEN 
        par_pos = try_lmp(i,tval) - boxval*floor(try_lmp(i,tval)&
             &/boxval)
     ELSE
        par_pos = trz_lmp(i,tval) - boxval*floor(trz_lmp(i,tval)&
             &/boxval)
     END IF

     flagdom = 0
     domnum  = 0

     IF(par_pos .GE. 0.0 .AND. par_pos .LT. inter(1)) THEN!Changed her
        !eon 2/8/16
        
        flagdom = 1
        domnum  = domtyp(1)
        
     END IF

     IF(flagdom == 0) THEN
        
        DO j = 1, maxinter-1
        
           IF(par_pos .GE. inter(j) .AND. par_pos .LT. inter(j+1)) THEN
              
              flagdom = 1
              domnum  = domtyp(j+1)
              
           END IF
           
        END DO
        
     END IF
     
     IF(flagdom == 0) THEN
        
        j = maxinter
        
        IF(par_pos .GE. inter(j) .AND. par_pos .LT. boxval) THEN
           
           flagdom = 1
           domnum  = domtyp(1)
           
        END IF
        
     END IF


     IF(flagdom == 0 .OR. domnum == 0 .OR. domnum > maxinter) THEN
        
        PRINT *, "Particle Coordinate not assigned domain properly"
        PRINT *, "i,Coordinate,domnum,tval", i,par_pos,domnum,tval
        
        WRITE(logout,*),"Particle Coordinate not assinged domain properl&
             &y"
        WRITE(logout,*), "i,Coordinate,domnum,time", i,par_pos,domnum&
             &,tval
     
        !STOP
        
     ELSE
        
        seg_dtype(i) = domnum
        
     END IF
     
  END DO
!$OMP END PARALLEL DO

  IF(tval == 1) THEN
     DO i = 1,ntotatoms
        WRITE(18,'(2(I0,1X),3(F14.8,1X))') i, seg_dtype(i)&
             &,trx_lmp(i,tval), try_lmp(i,tval), trz_lmp(i,tval)
        
     END DO

     CLOSE(18)
  END IF

  IF(lipclust .AND. tval == 1) THEN
     CALL LiP_ION_CLUSTERS_PEO_DOMAIN(tval)
  ELSEIF(lipclust .AND. mod(tval,clusfreq) == 0) THEN
     CALL LiP_ION_CLUSTERS_PEO_DOMAIN(tval)
  END IF


END SUBROUTINE DOMVAL

!--------------------------------------------------------------------

SUBROUTINE DOMCLASSIFY(tval)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j, flagdom, AllocateStatus, flagnew
  REAL :: par_pos
  INTEGER, DIMENSION(1:maxinter) :: domAdum, domBdum
  INTEGER, INTENT(IN) :: tval

  DO i = 1, maxinter
        
     domAdum(i) = 0
     domBdum(i) = 0
     
  END DO

  DO i = 1,npolyatoms
     
     IF(major_axis == 1) THEN
        par_pos = trx_lmp(i,tval) - boxval*floor(trx_lmp(i,tval)&
             &/boxval)
     ELSE IF(major_axis == 2) THEN 
        par_pos = try_lmp(i,tval) - boxval*floor(try_lmp(i,tval)&
             &/boxval)
     ELSE
        par_pos = trz_lmp(i,tval) - boxval*floor(trz_lmp(i,tval)&
             &/boxval)
     END IF

     flagdom = 0

     IF(par_pos .GE. 0.0 .AND. par_pos .LT. inter(1)) THEN
        
        flagdom = 1
        IF(aidvals(i,3) .LE. 5) THEN

           domAdum(1) = domAdum(1) + 1

        ELSEIF(aidvals(i,3) .LE. 9) THEN

           domBdum(1) = domBdum(1) + 1

        END IF

     END IF
     
     flagnew = 0

     IF(flagdom == 0) THEN

        DO j = 1, maxinter-1
     
           IF(par_pos > inter(j) .AND. par_pos .LT. inter(j+1) ) THEN
              
              IF(flagnew == 1 .AND. flagdom == 0) STOP "particle binne&
                   &d twice"

              flagdom = 1
              flagnew = 1

              IF(aidvals(i,3) .LE. 5) THEN
                 
                 domAdum(j+1) = domAdum(j+1) + 1
                 
              ELSEIF(aidvals(i,3) .LE. 9) THEN
                 
                 domBdum(j+1) = domBdum(j+1) + 1
                 
              END IF

           END IF
           
        END DO
        
     END IF

     IF(flagdom == 0) THEN

        j = maxinter

        IF(par_pos > inter(j) .AND. par_pos .LT. boxval ) THEN

           flagdom = 1
           
           IF(aidvals(i,3) .LE. 5) THEN
              
              domAdum(1) = domAdum(1) + 1
              
           ELSEIF(aidvals(i,3) .LE. 9) THEN
              
              domBdum(1) = domBdum(1) + 1
              
           END IF

           dompos(i) = 1

        END IF

     END IF

     IF(par_pos == 0.0000000) par_pos = 10**(-8)

     IF(flagdom == 0) THEN
        
        PRINT *, "Flag for domain bin ", flagdom
        PRINT *, "Particle Coordinate not binned properly"
        PRINT *, "i,Coordinate", i,par_pos

        WRITE(logout,*), "Flag for domain bin ", flagdom
        WRITE(logout,*), "Particle Coordinate not binned properly"
        WRITE(logout,*), "i,Coordinate", i,par_pos
       
        !STOP

     END IF

  END DO

  DO i = 1, maxinter
     
     IF(domAdum(i) > domBdum(i)) THEN
        
        domtyp(i) = 1
        
     ELSEIF(domBdum(i) > domAdum(i)) THEN
        
        domtyp(i) = 2
        
     ELSE
        
        PRINT *,"Something wrong in distribution"
        PRINT *,"Number of particles in",i,"dom", domAdum(i),&
             & domBdum(i)
        
        WRITE(logout,*),"Something wrong in distribution"
        WRITE(logout,*),"Number of particles in",i,"dom", domAdum(i),&
             & domBdum(i)
        
        !STOP
        
     END IF
     
  END DO
  
  WRITE(logout,*) "The domain types are "
  
  DO i = 1,maxinter
     
     IF(domtyp(i) == 1) THEN
        
        WRITE(logout,*) i,domtyp(i), "PS"
        
     ELSE
        
        WRITE(logout,*) i,domtyp(i), "EO"
        
     END IF
     
  END DO
     
  ALLOCATE(widdoms(maxinter),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate widdoms"
  
  i = 1
  widdoms(i) = inter(i) + (boxval-inter(maxinter))
  
  DO i = 1,maxinter-1
     
     widdoms(i+1) = inter(i+1)-inter(i)
     
  END DO
  
  WRITE(logout,*) "Widths of each domain"
  DO i= 1,maxinter
     
     WRITE(logout,*) i,widdoms(i)
     
  END DO
  
END SUBROUTINE DOMCLASSIFY
           
!--------------------------------------------------------------------

SUBROUTINE SORT_WRT_DOM(tval)

  USE ANALYZE_PSPEO_WITHIONS
  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,mon,typ
  INTEGER, INTENT(IN) :: tval
  INTEGER :: ctot
  INTEGER :: AllocateStatus,cptot,cltot,cotot

  IF(tval == 1) THEN

     WRITE(logout,*) "Cross Segmental analysis for", epsinit
     PRINT *, "Cross Segmental analysis for", epsinit
     
     WRITE(epsnum,'(I3)') int(epsinit)

  END IF

  ALLOCATE(crossflagarr(ntotion_centers+n_o_mons),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate crossflagarr"
  crossflagarr = -1

!Classify A and B monomers into same and cross segments
!IMPORTANT: seg_dtype stores information of particle. Hence size is
!totpart. To see the type of a particular particle, access using
!sortedarray(i,1), whereas sortedarray only consists of ions+oxy
!particles and to access the type of the specified partilce use the
!second index

!Flags: OinPS-1,OinEO-2,LinPS-3,LinEO-4,PinPS-5,PinEO-6
  tL_EO =0; tL_PS =0; tP_EO=0; tP_PS=0; tO_EO=0; tO_PS =0
!$OMP PARALLEL DO PRIVATE(i,mon) REDUCTION(+:tL_EO,tL_PS,tP_EO&
!$OMP& ,tP_PS,tO_EO,tO_PS)  
  DO i = 1,noxyplusions

     mon = sortedarray(i,1)
     !Sort O
     IF(sortedarray(i,2) == 9) THEN
        
        IF(seg_dtype(mon) == 1) THEN

           crossflagarr(i) = 1
           tO_PS = tO_PS + 1
           
        ELSEIF(seg_dtype(mon) == 2) THEN
           
           crossflagarr(i) = 2
           tO_EO = tO_EO + 1
           
        END IF
           
     !Sort Li
     ELSEIF(sortedarray(i,2) == 10) THEN

        IF(seg_dtype(mon) == 1) THEN

           crossflagarr(i) = 3
           tL_PS = tL_PS + 1

        ELSEIF(seg_dtype(mon) == 2) THEN
           
           crossflagarr(i) = 4
           tL_EO = tL_EO + 1

        END IF
        
     !Sort P   
     ELSEIF(sortedarray(i,2) == 11) THEN

        IF(seg_dtype(mon) == 1) THEN

           crossflagarr(i) = 5
           tP_PS = tP_PS + 1

        ELSEIF(seg_dtype(mon) == 2) THEN
           
           crossflagarr(i) = 6
           tP_EO = tP_EO + 1
           
        END IF

     ELSE

        PRINT *, "Wrong assignment in seg_dtype or seg_typ"
        PRINT *, i, mon, sortedarray(i,2), seg_dtype(mon)
        STOP

     END IF

  END DO
!$OMP END PARALLEL DO

  ctot = tO_PS + tO_EO + tL_PS + tL_EO + tP_PS + tP_EO
  cotot = tO_EO + tO_PS; cptot = tP_PS+tP_EO; cltot = tL_PS+tL_EO
  IF(ctot .NE. noxyplusions) THEN

     PRINT *, "Some problem in assigning particles to domains"
     PRINT *, tO_PS,tO_EO,tL_PS,tL_EO,tP_PS,tP_EO,ctot,noxyplusions
     WRITE(logout,*) "Some problem in assigning particles to domains"
     WRITE(logout,*) tO_PS,tO_EO,tL_PS,tL_EO,tP_PS,tP_EO,ctot&
          &,noxyplusions

     STOP

  ELSEIF(cotot .NE. n_o_mons) THEN

     PRINT *, "Some problem in assigning oxygens"
     PRINT *, tO_PS,tO_EO,cotot,n_o_mons
     WRITE(logout,*) "Some problem in assigning particles oxygens"
     WRITE(logout,*) tO_PS,tO_EO,cotot,n_o_mons

     STOP

  ELSEIF(cltot .NE. nLiions) THEN

     PRINT *, "Some problem in assigning Lithiums"
     PRINT *, tL_PS,tL_EO,cltot,nLiions
     WRITE(logout,*) "Some problem in assigning particles oxygens"
     WRITE(logout,*) tL_PS,tL_EO,cotot,nLiions

     STOP

  ELSEIF(cptot .NE. nPF6ions) THEN

     PRINT *, "Some problem in assigning oxygens"
     PRINT *, tP_PS,tP_EO,cPtot,nPF6ions
     WRITE(logout,*) "Some problem in assigning particles oxygens"
     WRITE(logout,*) tP_PS,tP_EO,cPtot,nPF6ions

     STOP

  ELSE
     
     IF(tval == 1) THEN

        PRINT *, "Number of particless in A and B domains are", &
        & tO_PS,tO_EO,tL_PS,tL_EO,tP_PS,tP_EO
        WRITE(logout,*) "Number of particles in A and B domai&
             &ns are",tO_PS,tO_EO,tL_PS,tL_EO,tP_PS,tP_EO

     END IF

  END IF

  IF(rdf_dom_ana) THEN
     
     IF(tval==1) THEN
        rdflocal = 0.0
!!$        CALL RDF3D_DOM_ANALYSIS(tO_PS,tO_EO,tL_PS,tL_EO,tP_PS&
!!$             &,tP_EO,tval)
        CALL RDF3D_DOM_ANALYSIS(tval)
     ELSEIF(mod(tval,rdffreq) == 0) THEN
!!$        CALL RDF3D_DOM_ANALYSIS(tO_PS,tO_EO,tL_PS,tL_EO,tP_PS&
!!$             &,tP_EO,tval)
        CALL RDF3D_DOM_ANALYSIS(tval)
     END IF

  END IF
  DEALLOCATE(crossflagarr)

END SUBROUTINE SORT_WRT_DOM

!--------------------------------------------------------------------

SUBROUTINE ION_DOM_SUBDIVIDE(cL_EO,cL_PS,cO_EO,cO_PS,cP_EO,cP_PS,tval)

  USE ANALYZE_PSPEO_WITHIONS
  IMPLICIT NONE

  INTEGER :: i,a1id,a1typ,AllocateStatus
  REAL :: par_pos
  INTEGER, INTENT(OUT) :: cL_EO,cL_PS,cO_EO,cO_PS,cP_EO,cP_PS
  INTEGER, INTENT(IN)  :: tval

  cL_EO = 0; cL_PS =0; cO_EO =0; cO_PS =0; cP_EO=0; cP_PS=0

!$OMP PARALLEL DO PRIVATE(i,a1id,a1typ,par_pos) REDUCTION(+:cL_EO&
!$OMP& ,cL_PS,cO_EO,cO_PS,cP_EO,cP_PS)
  DO i = 1,noxyplusions

     a1id   = sortedarray(i,1)     
     a1typ  = sortedarray(i,2)     

     IF(major_axis == 1) THEN
        par_pos = trx_lmp(i,tval) - boxval*floor(trx_lmp(i,tval)&
             &/boxval)
     ELSE IF(major_axis == 2) THEN 
        par_pos = try_lmp(i,tval) - boxval*floor(try_lmp(i,tval)&
             &/boxval)
     ELSE
        par_pos = trz_lmp(i,tval) - boxval*floor(trz_lmp(i,tval)&
             &/boxval)
     END IF
     
     IF(par_pos .LT. inter(2) .OR. par_pos .GE. inter(maxinter)) THEN
        
        IF(a1typ == 9 .AND. seg_dtype(a1id) == 1) THEN

           cO_PS = cO_PS + 1
           
        ELSEIF(a1typ == 9 .AND. seg_dtype(a1id) == 2) THEN

           cO_EO = cO_EO + 1

        ELSEIF(a1typ == 10 .AND. seg_dtype(a1id) == 1) THEN

           cL_PS = cL_PS + 1
           
        ELSEIF(a1typ == 10 .AND. seg_dtype(a1id) == 2) THEN

           cL_EO = cL_EO + 1

        ELSEIF(a1typ == 11 .AND. seg_dtype(a1id) == 1) THEN

           cP_PS = cP_PS + 1
           
        ELSEIF(a1typ == 11 .AND. seg_dtype(a1id) == 2) THEN

           cP_EO = cP_EO + 1

        ELSE

           print *, a1id, a1typ,seg_dtype(a1id)

           STOP "Unknown ion in further binning"

        END IF

      

     END IF

  END DO
!$OMP END PARALLEL DO


END SUBROUTINE ION_DOM_SUBDIVIDE

!--------------------------------------------------------------------

SUBROUTINE RDF3D_DOM_ANALYSIS(tval)

  USE ANALYZE_PSPEO_WITHIONS  
  IMPLICIT NONE
  

  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,intval
  REAL :: rxval,ryval,rzval,rval,rboxval,lEOdom,lPSdom
  REAL :: rlower,rupper,pos,ldist,rdist,rbase,drval
  REAL :: videal,ideal_cont
  INTEGER :: domid1,domid2,cL_EO,cL_PS,cO_EO,cO_PS,cP_EO,cP_PS
  INTEGER, INTENT(IN) :: tval
!  INTEGER,INTENT(IN) :: tO_PS,tO_EO,tL_PS,tL_EO,tP_PS,tP_EO,tval
  REAL,DIMENSION(0:rmaxbin-1,6) :: dumrdfarray
  REAL,DIMENSION(6) :: normfacarr
  REAL,DIMENSION(6) :: volarr
  INTEGER :: l1cntr,l2cntr,p1cntr,p2cntr,flag

!Add as necessary. Change dumrdfarray. Also rdflocal in ALLOCATE
!Num  - Type1  - Type2  (a1type-a2type-crossflaga1-crossflaga2)
!Col1 - Li(EO) - O(EO)  (10-9-2-2)
!Col2 - P(EO)  - O(EO)  (11-9-2-2)
!Col3 - Li(EO) - P(EO)  (10-11-2-2)
!Col4 - Li(PS) - 0(PS)  (10-9-1-1)
!Col5 - P(PS)  - O(PS)  (11-9-1-1)
!Col6 - Li(EO) - P(PS)  (10-11-1-1)

  lEOdom = 0.0; lPSdom = 0.0

  DO i = 1,maxinter

     IF(domtyp(i) == 1) lPSdom = lPSdom + widdoms(i)
     
     IF(domtyp(i) == 2) lEOdom = lEOdom + widdoms(i)

  END DO

  IF(tval == 1) PRINT *, "lPSdom/lEOdom", lPSdom, lEOdom

  rvolval = box_xl*box_yl*box_zl/boxval  


  CALL ION_DOM_SUBDIVIDE(cL_EO,cL_PS,cO_EO,cO_PS,cP_EO,cP_PS,tval)

  print *, cL_EO,cL_PS,cO_EO,cO_PS,cP_EO,cP_PS
  print *, tL_EO,tL_PS,tO_EO,tO_PS,tP_EO,tP_PS
 
  volarr(1) = rvolval*lEOdom
  volarr(2) = rvolval*lEOdom
  volarr(3) = rvolval*lEOdom
  volarr(4) = rvolval*lPSdom
  volarr(5) = rvolval*lPSdom
  volarr(6) = rvolval*lPSdom
  
!!$  print *, volarr

  normfacarr(1) = REAL(tL_EO*tO_EO)/volarr(1)
  normfacarr(2) = REAL(tP_EO*tO_EO)/volarr(2)
  normfacarr(3) = REAL(tL_EO*tP_EO)/volarr(3)
  normfacarr(4) = REAL(tL_PS*tO_PS)/volarr(4)
  normfacarr(5) = REAL(tP_PS*tO_PS)/volarr(5)
  normfacarr(6) = REAL(tL_PS*tP_PS)/volarr(6)

!!$  print *, normfacarr
!!$  pause;

!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j)
  DO i = 0,rmaxbin-1

     DO j = 1,6

        dumrdfarray(i,j) = 0.0
        
     END DO

  END DO

!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval&
!$OMP& ,ibin,domid1,domid2,rlower,rupper,pos,ldist,rdist,rbase&
!$OMP& ,drval,videal,ideal_cont) REDUCTION(+:dumrdfarray)
  DO i = 1,noxyplusions

     a1id   = sortedarray(i,1)     
     a1type = sortedarray(i,2)
     domid1 = seg_dtype(a1id)

     IF(a1type .LT. 9 .OR. a1type .GT. 11) THEN
        PRINT *, a1type,a1id,i
        STOP "No i type found"
     END IF

     IF(major_axis == 1) THEN
        pos = trx_lmp(a1id,tval) - box_xl*FLOOR(trx_lmp(a1id,tval)&
             &/box_xl)

     ELSEIF(major_axis == 2) THEN 
        pos = try_lmp(a1id,tval) - box_yl*FLOOR(try_lmp(a1id,tval)&
             &/box_yl)

     ELSEIF(major_axis == 3) THEN 
        pos = trz_lmp(a1id,tval) - box_zl*FLOOR(trz_lmp(a1id,tval)&
             &/box_zl)
        
     END IF

     flag = 0

     IF(pos .LT. inter(2) .OR. pos .GE.inter(maxinter)) flag=1

     IF(pos .LT. inter(1)) THEN
        ldist = pos + boxval - inter(maxinter)
        rdist = inter(1) - pos
        intval = 1
     ELSEIF(pos .GE. inter(maxinter)) THEN
        ldist = pos - inter(maxinter)
        rdist = inter(1) - (pos - boxval)
        intval = 1
     ELSE
        DO j = 1,maxinter-1
           IF(pos .GE. inter(j) .AND. pos .LT. inter(j+1)) THEN
              ldist = pos - inter(j)
              rdist = inter(j+1) - pos
              intval = j+1
           END IF
        END DO
     END IF

     IF(ldist .LT. 0.0 .OR. rdist .LT. 0.0) THEN
        PRINT *, "distances less than zero. Recheck"
        PRINT *, ldist, rdist, pos, intval
     END IF

     DO j = 1,noxyplusions

        a2id   = sortedarray(j,1)        
        a2type = sortedarray(j,2)
        domid2 = seg_dtype(a2id)

        IF(a2type .LT. 9 .OR. a2type .GT. 11) THEN
           PRINT *, a2type,a2id,j
           STOP "No j type found"
        END IF

        rxval = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
        ryval = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
        rzval = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/rdombinval)

        IF((ibin .LT. rmaxbin) .AND. (a1id .NE. a2id) .AND. (rval&
             & .LE. rdomcut)) THEN
     
           rlower = real(ibin)*rdombinval
           rupper = rlower + rdombinval
           videal = (rupper**3 - rlower**3)
           rbase  = 0.5*(rlower+rupper)
           drval  = rupper-rlower

           IF(rval .GT. ldist .AND. rval .GT. rdist) THEN
           
              PRINT *, "Warning: Assumes a shape which is not accounte&
                   &d for. Check rdomcut"
              PRINT *, a1id, a2id, pos, rval, ldist, rdist, rdomcut
              PRINT *, domid1, domid2
              PRINT *, trx_lmp(a1id,tval), trx_lmp(a2id,tval)
              PRINT *, trz_lmp(a1id,tval), trz_lmp(a2id,tval)
              PRINT *, rxval, ryval, rzval, tval
              STOP "Terminating abnormally RDF3D_DOM_ANALYSIS"
              
           END IF

           IF(rval .LT. ldist .AND. rval .LT. rdist) THEN 
              !Full Sphere
              
              ideal_cont = 3.0/(4.0*pival*videal)

              !Spherical Cap              
           ELSEIF(rval .GT. ldist) THEN

              ideal_cont = 1.0/(pi2val*rval*(rval+ldist)*drval)

           ELSEIF(rval .GT. rdist) THEN

              ideal_cont = 1.0/(pi2val*rval*(rval+rdist)*drval)

           ELSE

              PRINT *, "WARNING: Unknown volume contribution"
              PRINT *, a1id, a2id, pos, rval, ldist, rdist,rdomcut
              PRINT *, rxval, ryval, rzval,rval,tval
              STOP "Terminating abnormally"

           END IF

           !Col1 - Li(EO) - O(EO)  (10-9-2-2)
           IF(a1type == 10 .AND. domid1 == 2 .AND. a2type == 9 .AND. &
                & domid2 ==2) THEN        
              
              dumrdfarray(ibin,1) = dumrdfarray(ibin,1) + ideal_cont
              
              !Col2 - P(EO)  - O(EO)  (11-9-2-2)
           ELSEIF(a1type == 11 .AND. domid1 == 2 .AND. a2type == 9&
                & .AND. domid2 == 2) THEN

              dumrdfarray(ibin,2) = dumrdfarray(ibin,2) + ideal_cont
              
              !Col3 - Li(EO) - P(EO)  (10-11-2-2)
           ELSEIF(a1type == 10 .AND. domid1 == 2 .AND. a2type ==&
                & 11.AND. domid2 == 2) THEN

              dumrdfarray(ibin,3) = dumrdfarray(ibin,3) + ideal_cont
              
              !Col4 - Li(PS) - 0(PS)  (10-9-1-1)
           ELSEIF(a1type == 10 .AND. domid1 == 1 .AND.a2type == 9&
                & .AND. domid2 == 1) THEN

              dumrdfarray(ibin,4) = dumrdfarray(ibin,4) + ideal_cont
              
              !Col5 - P(PS)  - O(PS)  (11-9-1-1)
           ELSEIF(a1type == 11 .AND. domid1 == 1 .AND. a2type == 9&
                & .AND. domid2 == 1) THEN

              dumrdfarray(ibin,5) = dumrdfarray(ibin,5) + ideal_cont
              
              !Col6 - Li(PS) - P(PS)  (10-11-1-1)
           ELSEIF(a1type == 10 .AND. domid1 == 1 .AND. a2type ==&
                & 11.AND. domid2 == 1) THEN

              dumrdfarray(ibin,6) = dumrdfarray(ibin,6) + ideal_cont
              
           END IF

        END IF
        
     END DO

  END DO
!$OMP END DO


!$OMP DO PRIVATE(i,j)
  DO i = 0,rmaxbin-1

     DO j = 1,6

        rdflocal(i,j) = rdflocal(i,j) + REAL(dumrdfarray(i,j))/REAL(normfacarr(j))

     END DO

  END DO
!$OMP END DO

!$OMP END PARALLEL




!  PRINT *, "cntr", LiPcntr, LiOcntr, POcntr

END SUBROUTINE RDF3D_DOM_ANALYSIS

!--------------------------------------------------------------------
                                                                   
SUBROUTINE RDF2D_LAMELLAR_PLANE(tval,cLi,cP,cO,ipos,widbin)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,a1type,a2type,a1id,a2id,ibin
  REAL :: rxval,ryval,rzval,rval,rboxval
  INTEGER :: LiPcntr,LiOcntr,POcntr,lic,pc
  REAL, INTENT(IN) :: widbin
  INTEGER, INTENT(IN) :: tval,cLi,cP,cO,ipos
  INTEGER,DIMENSION(0:r2dmaxbin-1,3) :: dumrdfarray

!!$  r2dbinval = boxval/REAL(2.0*rmaxbin)
!!$
!!$  IF(ipos == 1) THEN
!!$
!!$     r2dbinavg = r2dbinavg + r2dbinval
!!$     rareavg   = rareavg + (box_xl*box_yl*box_zl)/(2.0*boxval) !Half the
!!$     ! box is considered
!!$
!!$  END IF

  rvolval   = box_xl*box_yl*box_zl
  
  IF(ipos == 1 .AND. rdfcalc .EQ. .false. .AND. oxylocal .EQ.&
       & .false.) THEN


     r2dvolavg = r2dvolavg + rvolval

  END IF

  dumrdfarray = 0; LiPcntr = 0; LiOcntr = 0; POcntr = 0
  
!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval,&
!$OMP& ibin) REDUCTION(+:LiPcntr,LiOcntr,POcntr,dumrdfarray)
  DO i = 1,num_mons

     a1id   = seg_act(i)     
     a1type = seg_typ(i)

     IF(a1type .NE. aidvals(a1id,3)) THEN
        PRINT *, "1", a1type, a1id, aidvals(a1id,3)
        STOP "Invalid j type found"
     END IF
     IF(a1type .LT. 9 .OR. a1type .GT. 11) THEN
        PRINT *, a1type,a1id,i
        STOP "Invalid i type found"
     END IF
     
     DO j = 1,n_o_mons+ntotion_centers
        
        a2id   = sortedarray(j,1)
        a2type = sortedarray(j,2)

        IF(a2type .LT. 9 .OR. a2type .GT. 11) THEN
           PRINT *, a2type,a2id,j
           STOP "Invalid j type found"
        END IF

        IF(a2type .NE. aidvals(a2id,3)) THEN
           PRINT *, "2",a2type, a2id, aidvals(a2id,3)
           STOP "Invalid j type found"
        END IF

        
!!$        IF(major_axis == 1) THEN
!!$           
!!$           r1val = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
!!$           r2val = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
!!$           
!!$           r1val = r1val - box_yl*ANINT(r1val/box_yl)
!!$           r2val = r2val - box_zl*ANINT(r2val/box_zl)
!!$           
!!$        ELSEIF(major_axis == 2) THEN
!!$           
!!$           r1val = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
!!$           r2val = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
!!$
!!$           r1val = r1val - box_xl*ANINT(r1val/box_xl)
!!$           r2val = r2val - box_zl*ANINT(r2val/box_zl)
!!$           
!!$        ELSEIF(major_axis == 3) THEN
!!$           
!!$           r1val = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
!!$           r2val = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
!!$           
!!$           r1val = r1val - box_xl*ANINT(r1val/box_xl)
!!$           r2val = r2val - box_yl*ANINT(r2val/box_yl)
!!$           
!!$        END IF
!!$        rval = sqrt(r1val**2 + r2val**2)
!!$        ibin = FLOOR(rval/r2dbinval)

        rxval = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
        ryval = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
        rzval = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        ibin = FLOOR(rval/r2dbinval)
        
        IF(a1type == 10 .AND. a2type == 11) THEN        
           
           IF(ibin .LT. r2dmaxbin) THEN
              
              dumrdfarray(ibin,1) = dumrdfarray(ibin,1) + 1
              LiPcntr   = LiPcntr + 1
              
           END IF
           
        ELSEIF(a1type == 10 .AND. a2type == 9) THEN
           
           IF(ibin .LT. r2dmaxbin) THEN

              dumrdfarray(ibin,2) = dumrdfarray(ibin,2) + 1
              LiOcntr   = LiOcntr + 1
              
           END IF
           
        ELSEIF(a1type == 11 .AND. a2type == 9) THEN
           
           IF(ibin .LT. r2dmaxbin) THEN
              
              dumrdfarray(ibin,3) = dumrdfarray(ibin,3) + 1
              POcntr   = POcntr + 1
              
           END IF
           
        END IF
           
     END DO

  END DO
!$OMP END DO


!$OMP DO PRIVATE(i)
  DO i = 0,r2dmaxbin-1

     rdf2darray(i,1,ipos) = rdf2darray(i,1,ipos) + REAL(dumrdfarray(i&
          &,1))*rvolval/(REAL(cli)*REAL(ptype))
     rdf2darray(i,2,ipos) = rdf2darray(i,2,ipos) + REAL(dumrdfarray(i&
          &,2))*rvolval/(REAL(cli)*REAL(otype))
     rdf2darray(i,3,ipos) = rdf2darray(i,3,ipos) + REAL(dumrdfarray(i&
          &,3))*rvolval/(REAL(cp)*REAL(otype))
  END DO
!$OMP END DO

!$OMP END PARALLEL

!!$  print *, rdf2darray(:,3,ipos), LiPcntr;pause;
!  PRINT *, tval,ipos,LiPcntr, LiOcntr, POcntr

END SUBROUTINE RDF2D_LAMELLAR_PLANE
  
!--------------------------------------------------------------------

SUBROUTINE RDF2D_LAMELLAR_PLANE_WRTDOM(tval,cLi,cP,cO,ipos,widbin)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,a1type,a2type,a1id,a2id,ibin
  REAL :: rxval,ryval,rzval,rparval,rboxval,par_pos
  REAL :: rlower,rupper,pos,ldist,rdist,rbase,drval
  REAL :: r1val,r2val,r3val,videal,ideal_cont,pos1,pos2,rvalfull
  REAL :: lEOdom, lPSdom, volEO, volPS
  INTEGER :: LiPcntr,LiOcntr,POcntr,lic,pc,domid1,domid2,intval
  REAL, INTENT(IN) :: widbin
  INTEGER, INTENT(IN) :: tval,cLi,cP,cO,ipos
  INTEGER :: flagmon, flagnew, mon
  INTEGER :: LiinEO, LiinPS, PinEO, PinPS
  REAL,DIMENSION(0:rmaxbin-1,6) :: dumrdfarray

  dumrdfarray = 0; LiPcntr = 0; LiOcntr = 0; POcntr = 0
  LiinEO = 0;LiinPS = 0; PinEO = 0; PinPS = 0

!  r2dbinval = boxval/REAL(2.0*rmaxbin)

!!$  IF(ipos == 1) THEN
!!$
!!$     r2dbinavg = r2dbinavg + r2dbinval
!!$     rareavg   = rareavg + (box_xl*box_yl*box_zl)/(2.0*boxval) !Half the
!!$     ! box is considered
!!$
!!$  END IF

  !1 - Li in PS - P in PS
  !2 - Li in PS - Li in PS
  !3 - P in PS  - P in PS
  !4 - Li in EO - P in EO
  !5 - Li in EO - Li in EO
  !6 - P in EO  - P in EO

  lEOdom = 0.0; lPSdom = 0.0
  DO i = 1,maxinter

     IF(domtyp(i) == 1) lPSdom = lPSdom + widdoms(i)
     
     IF(domtyp(i) == 2) lEOdom = lEOdom + widdoms(i)

  END DO

  IF(tval == 1) PRINT *, "lPSdom/lEOdom", lPSdom, lEOdom

  rvolval = box_xl*box_yl*box_zl/boxval  
  volEO   = rvolval*lEOdom
  volPS   = rvolval*lPSdom

  

!$OMP PARALLEL DO PRIVATE(i,a1id,a1type,domid1) REDUCTION(+:LiinEO&
!$OMP& ,LiinPS,PinEO,PinPS)

  DO i = 1,num_mons

     a1id   = seg_act(i)     
     a1type = seg_typ(i)
     domid1 = seg_dtype(a1id)
     IF(a1type == 10 .AND. domid1 == 1) LiinPS = LiinPS+1
     IF(a1type == 10 .AND. domid1 == 2) LiinEO = LiinEO+1
     IF(a1type == 11 .AND. domid1 == 1) PinPS  = PinPS+1
     IF(a1type == 11 .AND. domid1 == 2) PinEO  = PinEO+1
     
  END DO
!$OMP END PARALLEL DO

  IF(cLi .NE. LiinPS + LiinEO) STOP "Unequal Li count"
  IF(cP  .NE. PinPS  + PinEO) STOP "Unequal PF6 count"


!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rparval,rxval,ryval,rzval,intval&
!$OMP &,rbase,drval,ibin,r1val,r2val,r3val,pos,ldist,rdist,rlower,rvalfull&
!$OMP &,rupper,videal,ideal_cont,pos1,pos2,domid1,domid2)& 
!$OMP & REDUCTION(+:LiPcntr,LiOcntr,POcntr,dumrdfarray,lic,pc)
  DO i = 1,num_mons

     a1id   = seg_act(i)     
     a1type = seg_typ(i)
     domid1 = seg_dtype(a1id)

     IF(a1type .NE. aidvals(a1id,3)) THEN
        PRINT *, "1", a1type, a1id, aidvals(a1id,3)
        STOP "Invalid j type found"
     END IF
     IF(a1type .LT. 9 .OR. a1type .GT. 11) THEN
        PRINT *, a1type,a1id,i
        STOP "Invalid i type found"
     END IF

     IF(major_axis == 1) THEN
        pos = trx_lmp(a1id,tval) - box_xl*FLOOR(trx_lmp(a1id,tval)&
             &/box_xl)

     ELSEIF(major_axis == 2) THEN 
        pos = try_lmp(a1id,tval) - box_yl*FLOOR(try_lmp(a1id,tval)&
             &/box_yl)

     ELSEIF(major_axis == 3) THEN 
        pos = trz_lmp(a1id,tval) - box_zl*FLOOR(trz_lmp(a1id,tval)&
             &/box_zl)

     END IF

     IF(pos .LT. inter(1)) THEN
        ldist = pos + boxval - inter(maxinter)
        rdist = inter(1) - pos
        intval = 1
     ELSEIF(pos .GE. inter(maxinter)) THEN
        ldist = pos - inter(maxinter)
        rdist = inter(1) - (pos - boxval)
        intval = 1
     ELSE
        DO j = 1,maxinter-1
           IF(pos .GE. inter(j) .AND. pos .LT. inter(j+1)) THEN
              ldist = pos - inter(j)
              rdist = inter(j+1) - pos
              intval = j+1
           END IF
        END DO
     END IF

     IF(ldist .LT. 0.0 .OR. rdist .LT. 0.0) THEN
        PRINT *, "distances less than zero. Recheck"
        PRINT *, ldist, rdist, pos, intval
     END IF

     DO j = 1,noxyplusions
        
        a2id   = sortedarray(j,1)
        a2type = aidvals(a2id,3)
        domid2 = seg_dtype(a2id)

        IF(a2type .LT. 9 .OR. a2type .GT. 11) THEN
           PRINT *, a2type,a2id,j
           STOP "Invalid j type found"
        END IF

             
!!$        IF(major_axis == 1) THEN
!!$           
!!$           r1val = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
!!$           r2val = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
!!$           r3val = tr
!!$
!!$           r1val = r1val - box_yl*ANINT(r1val/box_yl)
!!$           r2val = r2val - box_zl*ANINT(r2val/box_zl)
!!$
!!$           pos1  = trx_lmp(a1id,tval) - box_xl*FLOOR(trx_lmp(a1id&
!!$                &,tval)/box_xl)
!!$           pos2  = trx_lmp(a2id,tval) - box_xl*FLOOR(trx_lmp(a2id&
!!$                &,tval)/box_xl)
!!$           r3val = pos1 - pos2
!!$           
!!$        ELSEIF(major_axis == 2) THEN
!!$           
!!$           r1val = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
!!$           r2val = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
!!$
!!$
!!$           r1val = r1val - box_xl*ANINT(r1val/box_xl)
!!$           r2val = r2val - box_zl*ANINT(r2val/box_zl)
!!$
!!$           pos1  = try_lmp(a1id,tval) - box_yl*FLOOR(try_lmp(a1id&
!!$                &,tval)/box_yl)
!!$           pos2  = try_lmp(a2id,tval) - box_yl*FLOOR(try_lmp(a2id&
!!$                &,tval)/box_yl)
!!$           r3val = pos1 - pos2
!!$           
!!$        ELSEIF(major_axis == 3) THEN
!!$           
!!$           r1val = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
!!$           r2val = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
!!$           
!!$           r1val = r1val - box_xl*ANINT(r1val/box_xl)
!!$           r2val = r2val - box_yl*ANINT(r2val/box_yl)
!!$
!!$           pos1  = trz_lmp(a1id,tval) - box_zl*FLOOR(trz_lmp(a1id&
!!$                &,tval)/box_zl)
!!$           pos2  = trz_lmp(a2id,tval) - box_zl*FLOOR(trz_lmp(a2id&
!!$                &,tval)/box_zl)
!!$           r3val = pos1 - pos2
!!$           
!!$        END IF

        r1val = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
        r2val = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
        r3val = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
        
        r1val = r1val - box_xl*ANINT(r1val/box_xl)
        r2val = r2val - box_yl*ANINT(r2val/box_yl)
        r3val = r3val - box_zl*ANINT(r3val/box_zl)
        
!        rparval = sqrt(r1val**2 + r2val**2)
        rvalfull = sqrt(r1val**2 + r2val**2 + r3val**2)
        ibin = FLOOR(rvalfull/rdombinval)

        IF(ibin .LT. rmaxbin .AND. a1id .NE. a2id .AND. domid1 ==&
             & domid2 .AND. rvalfull .LE. rdomcut) THEN

           rlower = real(ibin)*rdombinval
           rupper = rlower + rdombinval
           videal = (rupper**3 - rlower**3)
           rbase  = 0.5*(rlower+rupper)
           drval  = rupper-rlower

           IF(rvalfull .GT. ldist .AND. rvalfull .GT. rdist) THEN
           
              PRINT *, "Warning: Assumes a shape which is not accounte&
                   &d for. Check rdomcut"
              PRINT *, a1id, a2id, pos, rvalfull, ldist, rdist, rdomcut
              PRINT *, pos1, pos2, domid1, domid2
              PRINT *, trx_lmp(a1id,tval), trx_lmp(a2id,tval)
              PRINT *, trz_lmp(a1id,tval), trz_lmp(a2id,tval)
              PRINT *, r1val, r2val, r3val,rparval,tval
              STOP "Terminating abnormally for RDF2D_LAMELLAR_PLANE_WR&
                   &T_DOM"
              
           END IF
           
           
           IF(rvalfull .LT. ldist .AND. rvalfull .LT. rdist) THEN 
              !Full Sphere
              
              ideal_cont = 3.0/(4.0*pival*videal)

              !Spherical Cap              
           ELSEIF(rvalfull .GT. ldist) THEN

              ideal_cont = 1.0/(pi2val*rvalfull*(rvalfull+ldist)&
                   &*drval)

           ELSEIF(rvalfull .GT. rdist) THEN

              ideal_cont = 1.0/(pi2val*rvalfull*(rvalfull+rdist)&
                   &*drval)

           ELSE

              PRINT *, "WARNING: Unknown volume contribution"
              PRINT *, a1id, a2id, pos, rvalfull, ldist, rdist,rdomcut
              PRINT *, r1val, r2val, r3val,tval
              STOP "Terminating abnormally"

           END IF

!        rxval = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
!        ryval = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
!        rzval = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
        
!        rxval = rxval - box_xl*ANINT(rxval/box_xl)
!        ryval = ryval - box_yl*ANINT(ryval/box_yl)
!        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
!        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
!!$!        ibin = FLOOR(rval/r2dbinval)
        
           IF(a1type == 10 .AND. a2type == 11 .AND. domid1 == 1 .AND.&
                & domid2 == 1) THEN        
              !Li-P in PS
              dumrdfarray(ibin,1)= dumrdfarray(ibin,1) + ideal_cont
              
           ELSEIF(a1type == 10 .AND. a2type == 9 .AND. domid1 == 1 &
                & .AND. domid2 == 1) THEN
              !Li-O in PS
              dumrdfarray(ibin,2) = dumrdfarray(ibin,2) + ideal_cont
                      
           ELSEIF(a1type == 11 .AND. a2type == 9 .AND. domid1 == 1 &
                & .AND. domid2 == 1) THEN
              !P-O in PS
              
              dumrdfarray(ibin,3) = dumrdfarray(ibin,3) + ideal_cont
                 
           ELSEIF(a1type == 10 .AND. a2type == 11 .AND. domid1 == 2 &
                & .AND. domid2 == 2) THEN
              !Li-P in EO
              
              dumrdfarray(ibin,4) = dumrdfarray(ibin,4) + ideal_cont
              
           ELSEIF(a1type == 10 .AND. a2type == 9 .AND. domid1 == 2 &
                & .AND. domid2 == 2) THEN
              !Li-O in EO
              dumrdfarray(ibin,5) = dumrdfarray(ibin,5) + ideal_cont
                 
           ELSEIF(a1type == 11 .AND. a2type == 9 .AND. domid1 == 2 &
                & .AND. domid2 == 2) THEN
              !P - O in EO
              dumrdfarray(ibin,6) = dumrdfarray(ibin,6) + ideal_cont
           END IF

        END IF
              
     END DO
           
  END DO
!$OMP END DO
  

!$OMP DO PRIVATE(i)
  DO i = 0,rmaxbin-1

     rdf2dps(i,1,ipos) = rdf2dps(i,1,ipos) + REAL(dumrdfarray(i,1))&
          &*volPS/(REAL(LiinPS)*REAL(tP_PS))
     rdf2dps(i,2,ipos) = rdf2dps(i,2,ipos) + REAL(dumrdfarray(i,2))&
          &*volPS/(REAL(LiinPS)*REAL(tO_PS))
     rdf2dps(i,3,ipos) = rdf2dps(i,3,ipos) + REAL(dumrdfarray(i,3))&
          &*volPS/(REAL(PinPS)*REAL(tO_PS))
     rdf2deo(i,1,ipos) = rdf2deo(i,1,ipos) + REAL(dumrdfarray(i,4))&
          &*volEO/(REAL(LiinEO)*REAL(tP_EO))
     rdf2deo(i,2,ipos) = rdf2deo(i,2,ipos) + REAL(dumrdfarray(i,5))&
          &*volEO/(REAL(LiinEO)*REAL(tO_EO))
     rdf2deo(i,3,ipos) = rdf2deo(i,3,ipos) + REAL(dumrdfarray(i,6))&
          &*volEO/(REAL(PinEO)*REAL(tO_EO))

  END DO
!$OMP END DO


!$OMP END PARALLEL

!!$  IF(ipos == 3) PRINT *, rdf2deo(:,1,ipos)

!  PRINT *, tval,ipos,LiPcntr, LiOcntr, POcntr

END SUBROUTINE RDF2D_LAMELLAR_PLANE_WRTDOM

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLRDF()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,ierr,p
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdffrnorm
  CHARACTER(LEN =3) :: intnum

  IF(rdffreq == 1) rdffrnorm = nframes
  IF(rdffreq .NE. 1) rdffrnorm = INT(nframes/rdffreq) + 1

!!$  IF(rdf_2d) THEN
!!$     rareavg = rareavg/REAL(rdffrnorm)
!!$     r2dbinavg = r2dbinavg/REAL(rdffrnorm)
!!$  END IF
  
  IF(rdfcalc .OR. oxylocal) PRINT *, "rvolavg: ", rvolavg

  IF(rdf_2d) PRINT *, "r2dvolavg: ", r2dvolavg

  IF(rdfcalc == .true.) THEN
     OPEN(unit = 90,file =trim(rdf_fname),action="write",status&
          &="replace",iostat=ierr)
     IF(ierr /= 0) PRINT *, "Unknown rdf_filename"
  END IF

  IF(oxylocal == .true.) THEN
     OPEN(unit = 95,file =trim(oxy_fname),action="write",status&
          &="replace",iostat=ierr)
     IF(ierr /= 0) PRINT *, "Unknown oxyrdf_filename"
  END IF

  IF(rdf_2d == .true.) THEN

     DO p = 1,pfin
   
        WRITE(intnum,'(I0)') p
        dum_fname  = trim(rdf2d_fname)//"_"//trim(intnum)//".txt"

        OPEN(unit = 92,file = trim(dum_fname), status="replace",&
             & action = "write")

        DO i = 0,r2dmaxbin-1

           rlower = real(i)*r2dbinval
           rupper = rlower + r2dbinval
           nideal = vconst*(rupper**3 - rlower**3)
   
           WRITE(92,'(4(F16.9,1X))') 0.5*r2dbinval*(REAL(2*i+1))&
                &,rdf2darray(i,1,p)/(rdffrnorm*nideal),rdf2darray(i,2&
                &,p)/(rdffrnorm*nideal),rdf2darray(i,3,p)/(rdffrnorm*nideal)

!!$           rlower = real(i)*r2dbinavg
!!$           rupper = rlower + r2dbinavg
!!$           nideal = pi2val*(rupper**2 - rlower**2)
!!$           WRITE(92,'(4(F16.9,1X))') 0.5*r2dbinavg*(REAL(2*i+1))&
!!$                &,rareavg*rdf2darray(i,1,p)/(rdffrnorm*densmult&
!!$                &*nideal),rareavg*rdf2darray(i,2,p)/(rdffrnorm&
!!$                &*densmult*nideal),rareavg*rdf2darray(i,3,p)&
!!$                &/(rdffrnorm*densmult*nideal)
           
        END DO
        
        IF(rdf_2d) CLOSE(92)
        
     END DO

  END IF

  IF(rdf_2ddom == .true.) THEN

     DO p = 1,pfin

        WRITE(intnum,'(I0)') p
        dum_fname  = trim(rdf2d_dom_fname)//"_"//trim(intnum)//".txt"
        OPEN(unit = 92,file = trim(dum_fname), status="replace",&
             & action = "write")
        DO i = 0,rmaxbin-1

   
           WRITE(92,'(7(F16.9,1X))') 0.5*rdombinval*(REAL(2*i+1))&
                &,rdf2deo(i,1,p)/REAL(rdffrnorm),rdf2deo(i,2,p)&
                &/REAL(rdffrnorm),rdf2deo(i,3,p)/REAL(rdffrnorm)&
                &,rdf2dps(i,1,p)/REAL(rdffrnorm),rdf2dps(i,2,p)&
                &/REAL(rdffrnorm),rdf2dps(i,3,p)/REAL(rdffrnorm)
        END DO
           
        IF(rdf_2ddom) CLOSE(92)
        
     END DO

  END IF
  

  DO i = 0,rmaxbin-1

     rlower = real(i)*rbinval
     rupper = rlower + rbinval
     nideal = vconst*(rupper**3 - rlower**3)

     IF(rdfcalc == .true.) THEN

        WRITE(90,'(4(F16.9,1X))') 0.5*rbinval*(REAL(2*i+1)),&
             &rdfarray(i,1)/(rdffrnorm*nideal),rdfarray(i,2)&
             &/(rdffrnorm*nideal),rdfarray(i,3)/(rdffrnorm*nideal)

     END IF

     IF(oxylocal == .true.) THEN

        rlower = real(i)*roxybinval
        rupper = rlower + roxybinval
        nideal = vconst*(rupper**3 - rlower**3)

        WRITE(95,'(4(F16.9,1X))') 0.5*roxybinval*(REAL(2*i+1))&
             &,rdf_o_fb(i)/(rdffrnorm*nideal),rdf_o_ff(i)/(rdffrnorm&
             &*nideal),rdf_o_bb(i)/(rdffrnorm*nideal)
     END IF

  END DO

  
  IF(rdfcalc) CLOSE(90)
  IF(oxylocal) CLOSE(95)
  
END SUBROUTINE OUTPUT_ALLRDF

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLCROSSRDF()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,ierr,rdffrnorm

  IF(rdffreq == 1) rdffrnorm = nframes
  IF(rdffreq .NE. 1) rdffrnorm = INT(nframes/rdffreq) + 1
  
  OPEN(unit = 92,file = crossrdf_fname, status="replace", action = "wr&
       &ite",iostat=ierr)
  
  IF(ierr .ne. 0) PRINT *, "crossrdffile not found"

  DO i = 0,rmaxbin-1
     
     WRITE(92,'(7(F16.9,1X))') 0.5*rdombinval*(REAL(2*i+1))&
          &,rdflocal(i,1)/REAL(rdffrnorm),rdflocal(i,2)&
          &/REAL(rdffrnorm),rdflocal(i,3)/REAL(rdffrnorm),rdflocal(i&
          &,4)/REAL(rdffrnorm),rdflocal(i,5)/REAL(rdffrnorm)&
          &,rdflocal(i,6)/(rdffrnorm)     
     
  END DO
  
  CLOSE(92)

END SUBROUTINE OUTPUT_ALLCROSSRDF

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLCLUSTERS()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,frnorm,ierr
  REAL :: totlipclust,totp_liclust

  IF(clusfreq == 1) frnorm = nframes
  IF(clusfreq .NE. 1) frnorm = nframes/clusfreq + 1

  totlipclust = 0.0; totp_liclust = 0.0

  IF(lipclust .AND. all_tana) THEN
     OPEN(unit = 90,file =trim(lipclust_fname),action="write",status&
          &="replace",iostat=ierr)
     IF(ierr /= 0) PRINT *, "Unknown lipclust_filename"
!$OMP PARALLEL DO REDUCTION(+:totlipclust,totp_liclust) PRIVATE(i)
  
     DO i = 1,maxclustsize

        totlipclust  = totlipclust + REAL(lipclustavg(i))
        totp_liclust = totp_liclust + REAL(p_liclustavg(i))

     END DO

!$OMP END PARALLEL DO

     
     DO i = 1,maxclustsize     

        WRITE(90,'(I0,1X,6(F14.8,1X))') i, REAL(lipclustavg(i))&
             &/REAL(frnorm),100.0*REAL(lipclustavg(i))/totlipclust&
             &,REAL(p_liclustavg(i))/REAL(frnorm),100.0&
             &*REAL(p_liclustavg(i))/totp_liclust&
             &,REAL(lipclustavg_peo(i))/REAL(frnorm)&
             &,REAL(p_liclustavg_peo(i))/REAL(frnorm)

     END DO

  END IF

  IF(lipclust .AND. all_tana == .false.) THEN
     OPEN(unit = 90,file =trim(lipclust_fname),action="write",status&
          &="replace",iostat=ierr)
     IF(ierr /= 0) PRINT *, "Unknown lipclust_filename"
!$OMP PARALLEL DO REDUCTION(+:totlipclust,totp_liclust) PRIVATE(i)
  
     DO i = 1,maxclustsize

        totlipclust  = totlipclust + REAL(lipclustavg(i))
        totp_liclust = totp_liclust + REAL(p_liclustavg(i))

     END DO

!$OMP END PARALLEL DO

     
     DO i = 1,maxclustsize     

        WRITE(90,'(I0,1X,4(F14.8,1X))') i, REAL(lipclustavg(i))&
             &/REAL(frnorm),100.0*REAL(lipclustavg(i))/totlipclust&
             &,REAL(p_liclustavg(i))/REAL(frnorm),100.0&
             &*REAL(p_liclustavg(i))/totp_liclust

     END DO

     CLOSE(90)
  END IF

  IF(clus_ana) THEN

     OPEN(unit = 90,file =trim(clust_fname),action="write",status&
          &="replace",iostat=ierr)
     IF(ierr /= 0) PRINT *, "Unknown clust_filename"
     DO i = 1,ntotion_centers

        WRITE(90,'(I0,1X,F14.8,1X)') i, REAL(clust_avg(i))&
             &/REAL(nframes)

     END DO
     CLOSE(90)

  END IF
  

END SUBROUTINE OUTPUT_ALLCLUSTERS
     
!--------------------------------------------------------------------

SUBROUTINE OUTPUTDENS(distinit)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,ierr
  INTEGER, INTENT(IN) :: distinit

  IF(distinit ==  1) THEN
     OPEN(unit = 90,file ="density_avg_init.txt",action="write"&
          &,status="replace")
     binavg = binavg
     volavg = volavg
     boxlzavg = boxlzavg
     print *, "binavg_init: ", binavg
     print *, "volavg_init: ", volavg
     print *, "LZ-avg_init: ", boxlzavg
  ELSE
     OPEN(unit = 90,file =dens_fname,action="write",status="replace"&
          &,iostat=ierr)
     IF(ierr .NE. 0) PRINT *, "Unknown dens_fname"
     binavg = binavg/nframes
     volavg = volavg/nframes
     boxlzavg = boxlzavg/nframes
     print *, "binavg_fin: ", binavg
     print *, "volavg_fin: ", volavg
     print *, "LZ-avg_fin: ", boxlzavg
  END IF


  adensavg    = adensavg/densmult
  bdensavg    = bdensavg/densmult
  middensavg  = middensavg/densmult
  lithdensavg = lithdensavg/densmult
  pf6densavg  = pf6densavg/densmult

  IF(distinit == 1) THEN
     
     DO i = 0,maxbin-1
        
        WRITE(90,'(6(F16.9,1X))') 0.5*binavg*(REAL(2*i+1)),&
             & adensavg(i),bdensavg(i),middensavg(i)&
             &,lithdensavg(i),pf6densavg(i)
        
     END DO

  ELSE
     
     DO i = 0,maxbin-1
        
        WRITE(90,'(6(F16.9,1X))') 0.5*binavg*(REAL(2*i+1)),&
             & adensavg(i)/REAL(nframes),bdensavg(i)/REAL(nframes)&
             &,middensavg(i)/REAL(nframes),lithdensavg(i)&
             &/REAL(nframes),pf6densavg(i)/REAL(nframes)
        
     END DO
     
  END IF

  CLOSE(90)
 
END SUBROUTINE OUTPUTDENS

!--------------------------------------------------------------------

SUBROUTINE BINDING_LOCATION()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: i,j,tval,a1id,a2id,ierr,mapval
  REAL :: rxval,ryval,rzval,rval
!This algorith will work assuming the oxygens are numbered exactly the
!same way in each chain. ie., ID of any oxygen in SECOND chain =
!ID of SAME oxygen in FIRST chain + num_of_atomsperchain and so on
!If not one need to create mapping as was done for Santosh's PEO work
  INTEGER, DIMENSION(1:oxyperchain) :: oxy_index
  

  oxy_index = 0

!$OMP PARALLEL DO PRIVATE(tval,i,j,a1id,a2id,rxval,ryval,rzval,rval&
!$ &,mapval) REDUCTION(+:oxy_index)
  DO i = 1,n_o_mons

     a1id = otypearr(i)
     mapval = mod(i,oxyperchain)

     IF(mapval == 0) mapval = oxyperchain

     DO tval = 1,nframes

        DO j = 1,nLiions

           a2id = litypearr(j)

           rxval = trx_lmp(a1id,tval) - trx_lmp(a2id,tval) 
           ryval = try_lmp(a1id,tval) - try_lmp(a2id,tval) 
           rzval = trz_lmp(a1id,tval) - trz_lmp(a2id,tval) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LT. r_oli_cut) THEN

              oxy_index(mapval) = oxy_index(mapval) + 1

           END IF

        END DO

     END DO

  END DO

!$OMP END PARALLEL DO


  OPEN(unit = 42,file=lioxypos_file,action="write",status="replace"&
       &,iostat=ierr)

  IF(ierr /= 0) THEN

     PRINT *, lioxypos_file, "not found"
     STOP

  END IF

  DO i = 1,oxyperchain
           
     WRITE(42,"(I0,1X,F14.8)") i, REAL(oxy_index(i))/REAL(nframes)

  END DO


END SUBROUTINE BINDING_LOCATION

!--------------------------------------------------------------------

SUBROUTINE SPATIAL_ANALYSIS()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE
  
  INTEGER :: i,AllocateStatus
  INTEGER :: rdffrnorm

  IF(rdffreq == 1) rdffrnorm = nframes
  IF(rdffreq .NE. 1) rdffrnorm = INT(nframes/rdffreq) + 1

  box_xl = boxx_arr(1)
  box_yl = boxy_arr(1)
  box_zl = boxz_arr(1)
  boxval = boxvalarr(1)
  rvolval = box_xl*box_yl*box_zl

  CALL DOMCLASSIFY(1)

  DO i = 1,nframes

     box_xl = boxx_arr(i)
     box_yl = boxy_arr(i)
     box_zl = boxz_arr(i)
     boxval = boxvalarr(i)

     CALL COMPARTMENTALIZE_PARTICLES(i)

  END DO

  IF(rdf_2d) THEN
     r2dvolavg = r2dvolavg/REAL(rdffrnorm)
     CALL OUTPUT_ALLRDF()
  END IF

  IF(lipclust) CALL OUTPUT_ALLCLUSTERS()

END SUBROUTINE SPATIAL_ANALYSIS

!--------------------------------------------------------------------

SUBROUTINE STRUCTANALYSIS(frnum)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: frnum
  INTEGER :: t1,t2,clock_rate,clock_max

  IF(rdfcalc) THEN
     
     IF(frnum == 1) THEN           
        rdfarray = 0.0
        rvolavg = 0.0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL ALL_RDF()
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for RDF= ',REAL(t2-t1)/&
             & REAL(clock_rate)           
     END IF
     
     IF(mod(frnum,rdffreq) == 0) CALL ALL_RDF()
     
  END IF
  
  IF(oxylocal) THEN
     
     IF(frnum == 1) THEN
        
        rdf_o_fb = 0.0; rdf_o_ff=0.0; rdf_o_bb = 0.0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL SORTOXYFREECOMPLEX()
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for bound/free RDF= ',REAL(t2&
             &-t1)/REAL(clock_rate)           
     END IF
     
     IF(mod(frnum,rdffreq) == 0) CALL SORTOXYFREECOMPLEX()
     
  END IF
  
  IF(rgcalc)  CALL RGVALS(frnum)

  IF(lipclust) THEN
     
     IF(frnum == 1) THEN

        lipclustavg = 0.0; p_liclustavg=0.0
        lipclustavg_peo = 0.0; p_liclustavg_peo=0.0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL LiP_ION_CLUSTERS()
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for LiP-Clust Analysis= ',REAL(t2&
             &-t1)/REAL(clock_rate)           
     END IF
     
     IF(mod(frnum,clusfreq) == 0) CALL LiP_ION_CLUSTERS()
     
  END IF

  IF(clus_ana) THEN

     IF(frnum == 1) THEN

        clust_avg = 0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL CLUSTER_ANALYSIS(frnum)
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for Cluster Analysis= ',REAL(t2&
             &-t1)/REAL(clock_rate)           
     ELSE
        
        CALL CLUSTER_ANALYSIS(frnum)
     
     END IF

  END IF
        

END SUBROUTINE STRUCTANALYSIS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARRAYS()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(charge_lmp(ntotatoms,1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate charge_lmp"
!!$  ALLOCATE(image_lmp(ntotatoms,3),stat = AllocateStatus)
!!$  IF(AllocateStatus/=0) STOP "did not allocate image_lmp"
  ALLOCATE(vel_xyz(ntotatoms,4),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate vel_xyz"

  IF(ntotbonds /= 0) THEN
     ALLOCATE(bond_lmp(ntotbonds,4),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate bond_lmp"
  ELSE
     PRINT *, "Warning: No bonds - Not correct for bonded systems"
     ALLOCATE(bond_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(bond_lmp)
  END IF
  
  IF(ntotangls /= 0) THEN
     ALLOCATE(angl_lmp(ntotangls,5),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate angl_lmp"
  ELSE
     ALLOCATE(angl_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(angl_lmp)
  END IF
     
  IF(ntotdihds /= 0) THEN
     ALLOCATE(dihd_lmp(ntotdihds,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate dihd_lmp"
  ELSE
     ALLOCATE(dihd_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(dihd_lmp)
  END IF
  
  IF(ntotimprs /= 0) THEN
     ALLOCATE(impr_lmp(ntotimprs,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate zlmp"
  ELSE
     ALLOCATE(impr_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(impr_lmp)
  END IF

  ALLOCATE(dompos(ntotatoms),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate dompos"

  IF(all_tana) THEN
     ALLOCATE(trx_lmp(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate zlmp"
     ALLOCATE(try_lmp(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate zlmp"
     ALLOCATE(trz_lmp(ntotatoms,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate zlmp"
  ELSE
     ALLOCATE(trx_lmp(ntotatoms,nframes),stat = AllocateStatus)
     ALLOCATE(try_lmp(ntotatoms,nframes),stat = AllocateStatus)
     ALLOCATE(trz_lmp(ntotatoms,nframes),stat = AllocateStatus)
     DEALLOCATE(trx_lmp)
     DEALLOCATE(try_lmp)
     DEALLOCATE(trz_lmp)
  END IF

  IF(rdf_dom_ana == .true.) THEN

     ALLOCATE(rdflocal(0:rmaxbin-1,8),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdflocal"

  ELSE

     ALLOCATE(rdflocal(1,1),stat=AllocateStatus)
     DEALLOCATE(rdflocal)

  END IF

  PRINT *, "Successfully allocated memory"

END SUBROUTINE ALLOCATE_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE FINDDOMAIN(zval,boxw,dval)

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  REAL, INTENT(IN) :: zval,boxw
  INTEGER, INTENT(OUT) :: dval

  IF(zval .GE. 0.0 .AND. zval .LE. 18.50) THEN
     
     dval = 1

  ELSEIF(zval .GT. 18.5 .AND. zval .LE. 48.5) THEN

     dval = 2
     
  ELSEIF(zval .GT. 48.5 .AND. zval .LE. 72.86) THEN

     dval = 3

  ELSEIF(zval .GT. 72.86 .AND. zval .LE. 98.98) THEN

     dval = 4

  ELSEIF(zval .GT. 98.98 .AND. zval .LE. boxw) THEN

     dval = 1
     
  ELSE

     print *, "Out of bounds"

     STOP

  END IF

END SUBROUTINE FINDDOMAIN

!--------------------------------------------------------------------

SUBROUTINE ALLOUTPUTS()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  INTEGER :: rdffrnorm

  IF(rdffreq == 1) rdffrnorm = nframes
  IF(rdffreq .NE. 1) rdffrnorm = INT(nframes/rdffreq) + 1

  IF(rdfcalc == .true.) THEN
     rvolavg = rvolavg/REAL(rdffrnorm)
  END IF

  IF(rdfcalc) CALL OUTPUT_ALLRDF()

  IF(rdfcalc == .false. .AND. oxylocal == .true.) THEN
     rvolavg = rvolavg/REAL(rdffrnorm)
  END IF

  IF(oxylocal) CALL OUTPUT_ALLRDF()
  IF(lipclust .AND. all_tana == .false.) CALL OUTPUT_ALLCLUSTERS()
  IF(clus_ana) CALL OUTPUT_ALLCLUSTERS()

END SUBROUTINE ALLOUTPUTS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE ANALYZE_PSPEO_WITHIONS

  IMPLICIT NONE

  DEALLOCATE(aidvals)
  DEALLOCATE(rxyz_lmp)
  DEALLOCATE(vel_xyz)
!  DEALLOCATE(image_lmp)

  IF(ntotbonds /= 0) DEALLOCATE(bond_lmp)
  IF(ntotangls /= 0) DEALLOCATE(angl_lmp)
  IF(ntotdihds /= 0) DEALLOCATE(dihd_lmp)
  IF(ntotimprs /= 0) DEALLOCATE(impr_lmp)

END SUBROUTINE DEALLOCATE_ARRAYS

!--------------------------------------------------------------------
