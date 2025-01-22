!---------------To analyze properties of polymer-ion systems---------
!---------------Version 3: June-12-2024------------------------------
!---------------Parameter File: ana_params.f90-----------------------
!********************************************************************

PROGRAM PEMAIN

  USE ANALYZE_PARAMS
  IMPLICIT NONE

! Print headers

  PRINT *, "Static analysis of CG_VEC_MTFSI system .."
  PRINT *, "Starting OMP Threads .."
!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
  PRINT *, "Number of threads: ", nproc

! Call functions
  CALL READ_ANA_IP_FILE()
  CALL READ_DATAFILE()
  CALL SORTALLARRAYS()
  CALL ALLOCATE_ANALYSIS_ARRAYS()
  CALL ANALYZE_TRAJECTORYFILE()
  CALL ALLOUTPUTS()
  CALL DEALLOCATE_ARRAYS()
  
! Print completion
  PRINT *, "All Calculations Completed Succesfully :)"

END PROGRAM PEMAIN

!--------------------------------------------------------------------

SUBROUTINE READ_ANA_IP_FILE()

  USE ANALYZE_PARAMS

  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i,j
  CHARACTER(256) :: dumchar

  CALL DEFAULTVALUES()

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

     ELSEIF(dumchar == 'nframes') THEN

        READ(anaread,*,iostat=ierr) nframes

     ELSEIF(dumchar == 'skipfr') THEN

        READ(anaread,*,iostat=ierr) skipfr

     ELSEIF(dumchar == 'freqfr') THEN

        READ(anaread,*,iostat=ierr) freqfr
        
     ELSEIF(dumchar == 'nchains') THEN

        READ(anaread,*,iostat=ierr) nchains

     ELSEIF(dumchar == 'poly_types') THEN ! For Rg

        READ(anaread,*,iostat=ierr) npoly_types

        ALLOCATE(polytyp_arr(npoly_types),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate polytyp_arr"

        READ(anaread,*,iostat=ierr) (polytyp_arr(i),i=1,npoly_types)
        polyflag = 1

     ELSEIF(dumchar == 'ion_type') THEN
        
        READ(anaread,*,iostat=ierr) iontype
        
     ELSEIF(dumchar == 'cion_type') THEN

        READ(anaread,*,iostat=ierr) c_iontype

     ELSEIF(dumchar == 'iondiff') THEN

        READ(anaread,*,iostat=ierr) ion_diff
        ion_dynflag = 1

     ELSEIF(dumchar == 'ciondiff') THEN

        READ(anaread,*,iostat=ierr) cion_diff
        cion_dynflag = 1

     ELSEIF(dumchar == 'piondiff') THEN !stored in polyionarray

        READ(anaread,*,iostat=ierr) p_iontype
        pion_dynflag = 1

     ELSEIF(dumchar == 'compute_rdf') THEN

        rdfcalc = 1
        READ(anaread,*,iostat=ierr) rdffreq, rmaxbin, rdomcut,npairs
        
        ALLOCATE(pairs_rdf(npairs,3),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate pairs_rdf"
      
        DO i = 1,npairs

           READ(anaread,*,iostat=ierr) pairs_rdf(i,1), pairs_rdf(i,2)

        END DO

     ELSEIF(dumchar == 'compute_rg') THEN

        rgcalc = 1
        READ(anaread,*,iostat=ierr) rgfreq,rgall,rgavg

     ELSEIF(dumchar == 'catanneigh') THEN
        
        catan_neighcalc = 1
        READ(anaread,*,iostat=ierr) neighfreq,maxneighsize,rneigh_cut

     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword: ", trim(dumchar)
        STOP

     END IF

  END DO

  IF(logflag == 0) log_fname = "log."//trim(adjustl(traj_fname))
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  PRINT *, "Analysis input file read finished .."

  CALL SANITY_CHECK_IONTYPES()

END SUBROUTINE READ_ANA_IP_FILE

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  ! Frame, molecules and processor details
  nframes = 0; skipfr = 0; freqfr = 0; nfrcntr = 0
  nchains = 0; atperchain = 0

  ! Initialize flags
  polyflag = 0
  rgall = 0; rgcalc = 0; rdfcalc = 0
  ion_dynflag = 0; cion_dynflag = 0; pion_dynflag = 0
  ion_diff = 0; cion_diff = 0; pion_diff = 0

  ! Initialize iontypes
  c_iontype = -1; p_iontype = -1; iontype = -1

  !Initialize system quantities
  npoly_types = 0; ioncnt = 0; c_ioncnt = 0; p_ioncnt= 0

  ! Initialize distributions and frequencies
  rdffreq = 0; rgfreq = 0

  ! Initialzie structural quantities
  rdomcut = 0;  rmaxbin = 0

  ! Initialize structural averages
  rvolavg = 0; rgavg = 0

END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------

SUBROUTINE SANITY_CHECK_IONTYPES()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  IF(ion_dynflag .OR. catan_neighcalc) THEN

     IF(iontype == -1) THEN
        
        PRINT *, "ion type undefined for neigh or diff calculation"
        STOP 

     END IF

  END IF

  IF(cion_dynflag .OR. catan_neighcalc) THEN

     IF(c_iontype == -1) THEN
        
        PRINT *, "counter-ion type undefined for neigh or diff calculation"
        STOP 

     END IF

  END IF

  IF(pion_dynflag) THEN

     IF(p_iontype == -1) THEN
        
        PRINT *, "polymer-ion type undefined for diff calculation"
        STOP 

     END IF

  END IF


END SUBROUTINE SANITY_CHECK_IONTYPES

!--------------------------------------------------------------------

SUBROUTINE READ_DATAFILE()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,k,ierr,u,AllocateStatus,imax
  INTEGER :: flag, cntr, nwords
  INTEGER :: aid,molid,atype,ix,iy,iz
  REAL    :: charge,rx,ry,rz
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(256) :: rline,dumchar

  CALL COMPUTE_INIT_NLINES(imax)

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"

  WRITE(logout,*) "Datafile used: ", trim(adjustl(data_fname))

  ntotatoms = 0;ntotbonds=0;ntotangls=0;ntotdihds=0;ntotimprs=0
  atomflag =0;velflag = 0;bondflag=0;anglflag=0;dihdflag=0;imprflag=0

  READ(inpread,*)
  READ(inpread,*) 

  DO i = 1,imax-2 !Change here according to convenience
      
     READ(inpread,*) u, dumchar
     
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
     
  END DO

  READ(inpread,*)
  READ(inpread,*) xlo, xhi
  READ(inpread,*) ylo, yhi
  READ(inpread,*) zlo, zhi
  
  box_xl = xhi - xlo
  box_yl = yhi - ylo
  box_zl = zhi - zlo

  PRINT *, "x-box  ", "y-box  ", "z-box  "
  PRINT *, box_xl, box_yl, box_zl

  PRINT *, "STATISTICS"
  PRINT *, "Number of atoms/atomtypes: " , ntotatoms,ntotatomtypes
  PRINT *, "Number of bonds/bondtypes: " , ntotbonds,ntotbondtypes
  PRINT *, "Number of angles/angletypes: " , ntotangls,ntotangltypes
  PRINT *, "Number of diheds/dihedtypes: " , ntotdihds,ntotdihdtypes
  flag = 0; cntr = 0

  CALL ALLOCATE_TOPO_ARRAYS()

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

        END DO

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
           
           READ(inpread,*) aid,vel_xyz(aid,2),vel_xyz(aid,3)&
                &,vel_xyz(aid,4)
           vel_xyz(aid,1) = aid

        END DO

        CALL CHECK_MOMENTUM(0)

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
  
  PRINT *, "Fileread finished..."


END SUBROUTINE READ_DATAFILE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_INIT_NLINES(imax)

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: imax
  INTEGER :: init, pos, ipos,u,nwords,lcnt,ierr
  CHARACTER(LEN=120) :: charline

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"
  
  lcnt = 0

  READ(inpread,*)

  DO 

     READ(inpread,'(A)',iostat=ierr) charline     

     lcnt = lcnt + 1
     pos = 1
     nwords = 0

     DO

        ipos = VERIFY(charline(pos:),' ')
        IF(ipos == 0) EXIT
        nwords = nwords + 1
        pos = pos + ipos - 1
        ipos = SCAN(charline(pos:),' ')
        IF(ipos == 0) EXIT
        pos = pos + ipos - 1
        
     END DO

     IF(nwords .GE. 4) THEN

        imax = lcnt - 1
        EXIT
        
     END IF

  END DO

  CLOSE(inpread)

END SUBROUTINE COMPUTE_INIT_NLINES

!--------------------------------------------------------------------

SUBROUTINE ANALYZE_TRAJECTORYFILE()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,aid,ierr,atchk,atype,jumpfr,jout
  REAL :: xlo,xhi,ylo,yhi,zlo,zhi

  OPEN(unit = 15,file =trim(traj_fname),action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "trajectory file not found"

  PRINT *, "Trajectory file used: ",trim(adjustl(traj_fname))
  WRITE(logout,*) "Trajectory file used: "&
       &,trim(adjustl(traj_fname))

  
  PRINT *, "Analyzing trajectory file..."

  CALL STRUCT_INIT()
  CALL OPEN_STRUCT_OUTPUT_FILES()

  DO i = 1,skipfr

     DO j = 1,ntotatoms+9

        READ(15,*) 

     END DO

     IF(mod(i,100) == 0) PRINT *, "Skipped ", i, "frames"

  END DO

  DO i = 1,nframes

     nfrcntr = nfrcntr + 1
     IF(mod(i,100) == 0) PRINT *, "Processing ", i+1,"th frame"

     READ(15,*)
     READ(15,*) timestep

     READ(15,*) 
     READ(15,*) atchk

     READ(15,*) 
     READ(15,*) xlo, xhi
     READ(15,*) ylo, yhi
     READ(15,*) zlo, zhi

     READ(15,*)

     box_xl = xhi - xlo
     box_yl = yhi - ylo
     box_zl = zhi - zlo
     
     boxx_arr(i)  = box_xl
     boxy_arr(i)  = box_yl
     boxz_arr(i)  = box_zl

     DO j = 1,atchk

        READ(15,*) aid,atype,rxyz_lmp(aid,1),rxyz_lmp(aid,2)&
             &,rxyz_lmp(aid,3)

        IF(atype .NE. aidvals(aid,3)) THEN

           PRINT *, "Incorrect atom ids"
           PRINT *, i,j,aid,atype,aidvals(aid,3)
           STOP

        END IF

        ! Store to time-dependent array before subtracting xlo/ylo/zlo
        ! or else it will induce an error
        IF(ion_dynflag == 1 .AND. atype == iontype) THEN
                 
           CALL MAP_REFTYPE(aid,atype,jout)
           itrx_lmp(jout,i) = rxyz_lmp(aid,1)
           itry_lmp(jout,i) = rxyz_lmp(aid,2)
           itrz_lmp(jout,i) = rxyz_lmp(aid,3)

        ELSEIF(cion_dynflag == 1 .AND. atype == c_iontype) THEN

           CALL MAP_REFTYPE(aid,atype,jout)
           ctrx_lmp(jout,i) = rxyz_lmp(aid,1)
           ctry_lmp(jout,i) = rxyz_lmp(aid,2)
           ctrz_lmp(jout,i) = rxyz_lmp(aid,3)

        ELSEIF(pion_dynflag == 1 .AND. atype == p_iontype) THEN
           
           CALL MAP_REFTYPE(aid,atype,jout)
           ptrx_lmp(jout,i) = rxyz_lmp(aid,1)
           ptry_lmp(jout,i) = rxyz_lmp(aid,2)
           ptrz_lmp(jout,i) = rxyz_lmp(aid,3)

        END IF       

        rxyz_lmp(j,1) = rxyz_lmp(j,1) - xlo
        rxyz_lmp(j,2) = rxyz_lmp(j,2) - ylo
        rxyz_lmp(j,3) = rxyz_lmp(j,3) - zlo

     END DO
     
     IF(i == 1) PRINT *, "Beginning statics analysis..."
     CALL STRUCT_MAIN(nfrcntr)
     
     DO jumpfr = 1,freqfr

        READ(15,*)
        READ(15,*)        
        READ(15,*)
 
        READ(15,*) atchk

        DO j = 1,atchk+5

           READ(15,*) 

        END DO

     END DO

  END DO

  CLOSE(15)

  PRINT *, "Beginning dynamical analysis..."
  CALL DYNAMICS_MAIN()

END SUBROUTINE ANALYZE_TRAJECTORYFILE

!--------------------------------------------------------------------

SUBROUTINE OPEN_STRUCT_OUTPUT_FILES()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  IF(rgcalc) THEN
     
     IF(rgavg) THEN
        dum_fname = "rgavg_"//trim(adjustl(traj_fname))
        OPEN(unit = rgavgwrite,file =trim(dum_fname),action="write"&
             &,status="replace")
        WRITE(rgavgwrite,'(A5,1X,4(A3,1X))') "tstep", "rgt","rgx","&
             &rgy","rgz"

     END IF

     IF(rgall) THEN
        dum_fname = "rgall_"//trim(adjustl(traj_fname))
        OPEN(unit = rgwrite,file =trim(dum_fname),action="write"&
             &,status="replace")
     END IF
 
  END IF

   
END SUBROUTINE OPEN_STRUCT_OUTPUT_FILES

!--------------------------------------------------------------------

SUBROUTINE STRUCT_INIT()

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,t1,t2,norm,acnt,fcnt,a1id,molid,flagch,flagpr,jmax
  INTEGER :: AllocateStatus

  IF(rdfcalc) THEN

     rdfarray = 0.0
     rbinval = rdomcut/REAL(rmaxbin)

     DO i = 1, npairs

        t1 = 0; t2 = 0

        DO j = 1,ntotatoms

           IF(aidvals(j,3) == pairs_rdf(i,1)) t1 = t1+1
           IF(aidvals(j,3) == pairs_rdf(i,2)) t2 = t2+1
           
        END DO

        IF(pairs_rdf(i,1) == pairs_rdf(i,2)) THEN
           pairs_rdf(i,3) = t1*(t1-1) !g_AA(r)
        ELSE
           pairs_rdf(i,3) = t1*t2 !g_AB(r)
        END IF

     END DO

  END IF



END SUBROUTINE STRUCT_INIT

!--------------------------------------------------------------------

SUBROUTINE STRUCT_MAIN(tval)

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN):: tval
  INTEGER :: t1, t2
  INTEGER :: clock_rate, clock_max

  IF(rgcalc .AND. mod(tval-1,rgfreq)==0) CALL COMPUTE_RADGYR(tval)
  IF(rdfcalc .AND. mod(tval-1,rdffreq)==0) CALL COMPUTE_RDF(tval)

  IF(catan_neighcalc) THEN
     
     IF(tval == 1) THEN

        cat_an_neighavg = 0.0; an_cat_neighavg=0.0
        CALL SYSTEM_CLOCK(t1,clock_rate,clock_max)
        CALL CAT_AN_NEIGHS()
        CALL SYSTEM_CLOCK(t2,clock_rate,clock_max)
        PRINT *, 'Elapsed real time for Neighbor Analysis= ',REAL(t2&
             &-t1)/REAL(clock_rate)

     END IF
     
     IF(mod(tval,neighfreq) == 0) CALL CAT_AN_NEIGHS()
     
  END IF


END SUBROUTINE STRUCT_MAIN

!--------------------------------------------------------------------

SUBROUTINE CAT_AN_NEIGHS()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,a1id,a2id,neigh_cnt,tid
  INTEGER,DIMENSION(1:maxneighsize,0:nproc-1) :: cat_an_neigh_inst&
       &,an_cat_neigh_inst 
  REAL :: rxval, ryval, rzval, rval

  cat_an_neigh_inst = 0; an_cat_neigh_inst = 0

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,ioncnt
     
     neigh_cnt = 0
     a1id = ionarray(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,c_ioncnt

        a2id = counterarray(j,1)
        
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rneigh_cut) THEN
           
           neigh_cnt = neigh_cnt + 1
           
        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     cat_an_neigh_inst(neigh_cnt+1,tid) = cat_an_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(i,j,a1id,a2id,rxval,ryval,rzval,rval,neigh_cnt,tid)
!$OMP DO
  DO i = 1,c_ioncnt
     
     neigh_cnt = 0
     a1id = counterarray(i,1)
     tid = OMP_GET_THREAD_NUM()

     DO j = 1,ioncnt

        a2id = ionarray(j,1)
        
        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval - box_zl*ANINT(rzval/box_zl)
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LT. rneigh_cut) THEN
           
           neigh_cnt = neigh_cnt + 1
           
        END IF

     END DO

     IF(neigh_cnt + 1 .GT. maxneighsize) THEN

        PRINT *, "Neighbor count exceeded max size"
        PRINT *, neigh_cnt, maxneighsize
        STOP

     END IF

     an_cat_neigh_inst(neigh_cnt+1,tid) = an_cat_neigh_inst(neigh_cnt&
          &+1,tid) + 1

  END DO
!$OMP END DO


!$OMP DO 
  DO  i = 1,maxneighsize
     DO j = 0,nproc-1
        cat_an_neighavg(i) = cat_an_neighavg(i) + cat_an_neigh_inst(i&
             &,j)
        an_cat_neighavg(i) = an_cat_neighavg(i) + an_cat_neigh_inst(i&
             &,j)
     END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE CAT_AN_NEIGHS

!--------------------------------------------------------------------

SUBROUTINE DYNAMICS_MAIN()

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  IF(ion_diff) CALL DIFF_IONS()
  IF(cion_diff) CALL DIFF_COUNTERIONS()

END SUBROUTINE DYNAMICS_MAIN

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RDF(iframe)

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,paircnt,AllocateStatus
  REAL :: rxval,ryval,rzval,rval
  INTEGER :: a1ref,a2ref
  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  rvolavg = rvolavg + rvolval 

  ALLOCATE(dumrdfarray(0:rmaxbin-1,npairs),stat=AllocateStatus)
  IF(AllocateStatus/=0) STOP "dumrdfarray not allocated"
  dumrdfarray = 0


!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval&
!$OMP& ,ibin,paircnt,a1ref,a2ref) REDUCTION(+:dumrdfarray)
  DO paircnt = 1,npairs

     a1ref = pairs_rdf(paircnt,1); a2ref = pairs_rdf(paircnt,2)
     
     DO i = 1,ntotatoms
        
        a1id   = aidvals(i,1)     
        a1type = aidvals(i,3)
        
        DO j = 1,ntotatoms

           a2id   = aidvals(j,1)        
           a2type = aidvals(j,3)

           ! Remove identical IDs when computing g_AA(r)
           IF(a1id == a2id .AND. a1ref == a2ref) CYCLE


           IF(a1type == a1ref .AND. a2type == a2ref) THEN

              rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
              ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
              rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
              
              rxval = rxval - box_xl*ANINT(rxval/box_xl)
              ryval = ryval - box_yl*ANINT(ryval/box_yl)
              rzval = rzval - box_yl*ANINT(rzval/box_zl)
           
              rval = sqrt(rxval**2 + ryval**2 + rzval**2)
              ibin = FLOOR(rval/rbinval)
           
              IF(ibin .LT. rmaxbin) THEN
                 
                 dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                      &,paircnt) + 1
              
              END IF

           END IF
           
        END DO

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j)
  DO j = 1,npairs
     
     DO i = 0,rmaxbin-1

        rdfarray(i,j) = rdfarray(i,j) + REAL(dumrdfarray(i,j))&
             &*rvolval/(REAL(pairs_rdf(j,3)))
        
     END DO
     
  END DO
!$OMP END DO

!$OMP END PARALLEL

  DEALLOCATE(dumrdfarray)

END SUBROUTINE COMPUTE_RDF

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RADGYR(iframe)

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,molid,atype
  REAL    :: rgsqavg, rgxxavg, rgyyavg, rgzzavg
  REAL, DIMENSION(1:nchains) :: rgxx, rgyy, rgzz, rgsq
  REAL, DIMENSION(1:nchains) :: rxcm, rycm, rzcm, totmass
  INTEGER, INTENT(IN) :: iframe

  rgxx = 0.0; rgyy =0.0; rgzz = 0.0; rgsq = 0.0
  rgsqavg = 0.0; rgxxavg = 0.0; rgyyavg = 0.0; rgzzavg = 0.0
  rxcm = 0.0; rycm = 0.0; rzcm = 0.0
  totmass = 0.0


  IF(iframe == 1) THEN
     IF(polyflag == 0 .OR. npoly_types == 0) STOP "ERROR: No polymer t&
          &ypes defined!"
     PRINT *, masses
  END IF

!$OMP PARALLEL 
!$OMP DO PRIVATE(i,molid,atype) REDUCTION(+:totmass,rxcm,rycm,rzcm)

  DO i = 1,ntotatoms

     molid = aidvals(i,2)
     atype = aidvals(i,3)

     IF(ANY(polytyp_arr == atype)) THEN
       
        totmass(molid) = totmass(molid) + masses(atype,1)
        rxcm(molid) = rxcm(molid)+ rxyz_lmp(i,1)*masses(atype,1)
        rycm(molid) = rycm(molid)+ rxyz_lmp(i,2)*masses(atype,1)
        rzcm(molid) = rzcm(molid)+ rxyz_lmp(i,3)*masses(atype,1)

     END IF

  END DO
!$OMP END DO

!$OMP DO PRIVATE(i)
  DO i = 1, nchains
     rxcm(i) = rxcm(i)/totmass(i)
     rycm(i) = rycm(i)/totmass(i)
     rzcm(i) = rzcm(i)/totmass(i)
  END DO
!$OMP END DO


!$OMP DO PRIVATE(i,molid,atype) REDUCTION(+:rgxx, rgyy, rgzz, rgsq)
  DO i = 1,ntotatoms

     molid = aidvals(i,2)
     atype = aidvals(i,3)

     IF(ANY(polytyp_arr == atype)) THEN

        rgxx(molid) = rgxx(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
             &-rxcm(molid))**2)
        rgyy(molid) = rgyy(molid) + masses(atype,1)*((rxyz_lmp(i,2)&
             &-rycm(molid))**2)
        rgzz(molid) = rgzz(molid) + masses(atype,1)*((rxyz_lmp(i,3)&
             &-rzcm(molid))**2)

        rgsq(molid) = rgsq(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
             &-rxcm(molid))**2 + (rxyz_lmp(i,2)-rycm(molid))**2 +&
             & (rxyz_lmp(i,3)-rzcm(molid))**2)


     END IF

  END DO
!$OMP END DO

!$OMP DO PRIVATE(i)
  DO i = 1, nchains
     rgxx(i) = rgxx(i)/totmass(i)
     rgyy(i) = rgyy(i)/totmass(i)
     rgzz(i) = rgzz(i)/totmass(i)
     rgsq(i) = rgsq(i)/totmass(i)
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF(iframe == 1) THEN

     OPEN(unit = 98,file ="totmasschk.txt",action="write",status="repl&
          &ace")

     DO i = 1,nchains

        WRITE(98,'(I0,1X,8(F14.5,1X))') i, totmass(i),rxcm(i),&
             & rycm(i), rzcm(i), rgxx(i), rgyy(i), rgzz(i), rgsq(i)

     END DO

     CLOSE(98)

     OPEN(unit = 98,file ="molidchk.txt",action="write",status="repl&
          &ace")

     DO i = 1,ntotatoms

        WRITE(98,'(I0,1X,I0)') i, aidvals(i,2)

     END DO

     CLOSE(98)

  END IF


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
  

  IF(rgavg) THEN

     WRITE(rgavgwrite,'(I0,1X,4(F14.6,1X))') timestep, sqrt(rgsqavg),&
          & sqrt(rgxxavg), sqrt(rgyyavg), sqrt(rgzzavg)

  END IF

  IF(rgall) THEN
     
     WRITE(rgwrite,'(2(I0,1X),4(A3,1X))') timestep, nchains,"rgt","rg&
          &x","rgy","rgz"

     DO i = 1,nchains
     
        WRITE(rgwrite,'(I0,1X,4(F14.5,1X))') i,sqrt(rgsq(i))&
             &,sqrt(rgxx(i)),sqrt(rgyy(i)),sqrt(rgzz(i))

     END DO

  END IF
     
END SUBROUTINE COMPUTE_RADGYR

!--------------------------------------------------------------------

SUBROUTINE ALLOUTPUTS()

  USE ANALYZE_PARAMS
  IMPLICIT NONE
  INTEGER :: i

  PRINT *, "Number of frames from start to end: ", nframes/(freqfr+1)
  PRINT *, "Frequency of Frames: ", freqfr + 1
  PRINT *, "Total number of Frames analyzed: ", nfrcntr

  WRITE(logout,*) "Number of frames from start to end: ", nframes&
       &/(freqfr+1)
  WRITE(logout,*) "Frequency of Frames: ", freqfr+1
  WRITE(logout,*) "Total number of Frames analyzed: ", nfrcntr

  IF(rdfcalc) THEN
     PRINT *, "Writing RDFs .."
     CALL OUTPUT_ALLRDF()
  END IF

  IF(catan_neighcalc) THEN
     PRINT *, "Writing neighbors .."
     CALL OUTPUT_ALLNEIGHBORS()
  END IF


END SUBROUTINE ALLOUTPUTS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLRDF()

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdffrnorm,acrnorm
  
  IF(rdfcalc) THEN

     rdffrnorm = INT(nfrcntr/rdffreq)
     rvolavg = rvolavg/REAL(rdffrnorm)
     PRINT *, "Average volume of box", rvolavg
     
     dum_fname = "rdf_"//trim(adjustl(traj_fname))
     OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
          &,status="replace",iostat=ierr)
     
     IF(ierr /= 0) THEN
        PRINT *, "Could not open", trim(dum_fname)
     END IF
     
     WRITE(dumwrite,'(A,2X)',advance="no") "r"
     
     DO j = 1,npairs
        
        WRITE(dumwrite,'(2(I0,1X))',advance="no") pairs_rdf(j,1)&
             &,pairs_rdf(j,2)
        
     END DO
     
     WRITE(dumwrite,*)
     
     DO i = 0,rmaxbin-1
        
        rlower = real(i)*rbinval
        rupper = rlower + rbinval
        nideal = vconst*(rupper**3 - rlower**3)
        
        WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*rbinval*(REAL(2*i&
             &+1))
        
        DO j = 1,npairs
           
           WRITE(dumwrite,'(F16.9,1X)',advance="no")rdfarray(i,j)&
                &/(rdffrnorm*nideal)
           
        END DO
        
        WRITE(dumwrite,*)
        
     END DO
     
     CLOSE(dumwrite)

  END IF

END SUBROUTINE OUTPUT_ALLRDF

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLNEIGHBORS()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER :: i,frnorm,ierr
  REAL :: totcat_an_neigh,totan_cat_neigh

  IF(neighfreq == 1) frnorm = nframes
  IF(neighfreq .NE. 1) frnorm = nframes/neighfreq + 1

  totcat_an_neigh = 0.0; totan_cat_neigh = 0.0

  IF(catan_neighcalc) THEN
!$OMP PARALLEL DO REDUCTION(+:totcat_an_neigh,totan_cat_neigh) PRIVATE(i)
  
     DO i = 1,maxneighsize

        totcat_an_neigh = totcat_an_neigh + REAL(cat_an_neighavg(i))
        totan_cat_neigh = totan_cat_neigh + REAL(an_cat_neighavg(i))

     END DO

!$OMP END PARALLEL DO

     IF(catan_neighcalc) THEN
        
        dum_fname = "catanneigh_"//trim(adjustl(traj_fname))
        OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
             &,status="replace")

        
     END IF

     
     DO i = 1,maxneighsize     

        WRITE(dumwrite,'(I0,1X,4(F14.8,1X))') i,&
             & REAL(cat_an_neighavg(i))/REAL(frnorm),100.0&
             &*REAL(cat_an_neighavg(i))/totcat_an_neigh&
             &,REAL(an_cat_neighavg(i))/REAL(frnorm),100.0&
             &*REAL(an_cat_neighavg(i))/totan_cat_neigh

     END DO

     CLOSE(dumwrite)

  END IF

END SUBROUTINE OUTPUT_ALLNEIGHBORS
     
!--------------------------------------------------------------------
! pion_type is set to -1 in the beginning and is read if and only if
! pion_dynflag is activated
SUBROUTINE SORTALLARRAYS()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,a1type,cnt,AllocateStatus
  INTEGER, DIMENSION(1:ntotatoms,2) :: dumsortarr,dumcionarr&
       &,dumpionarr

  dumsortarr = -1; dumcionarr = -1; dumpionarr = -1
  cnt = 0

  DO i = 1,ntotatoms

     a1type = aidvals(i,3)

     IF(a1type == iontype) THEN
        ioncnt = ioncnt + 1
        dumsortarr(ioncnt,1) = i
        dumsortarr(ioncnt,2) = a1type

     ELSEIF(a1type == c_iontype) THEN
        c_ioncnt = c_ioncnt + 1
        dumcionarr(c_ioncnt,1) = i
        dumcionarr(c_ioncnt,2) = a1type

     ELSEIF(a1type == p_iontype) THEN
        p_ioncnt = p_ioncnt + 1
        dumpionarr(p_ioncnt,1) = i
        dumpionarr(p_ioncnt,2) = a1type

     END IF

  END DO

  IF(catan_neighcalc .OR. ion_dynflag .OR. cion_dynflag) THEN
     PRINT *, "Number of atoms of ion type: ", ioncnt
     PRINT *, "Number of atoms of cntion type: ", c_ioncnt
     
     ALLOCATE(ionarray(ioncnt,2),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ionarray"
     ALLOCATE(counterarray(c_ioncnt,2),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ionarray"
  
     ! Load ion array

     i = 0

     DO WHILE(dumsortarr(i+1,1) .NE. -1) 

        i = i + 1
        ionarray(i,1) = dumsortarr(i,1)
        ionarray(i,2) = dumsortarr(i,2)
        
     END DO

     IF(i .NE. ioncnt) THEN
        PRINT *, i, ioncnt
        STOP "Wrong total count in ionarray"
     END IF
     
     DO i = 1,ioncnt
     
        IF(ionarray(i,1) == -1 .OR. ionarray(i,2) == -1) THEN
           
           PRINT *, i,ionarray(i,1), ionarray(i,2)
           PRINT *, "Something wrong in assigning ionarray"
           STOP
           
        END IF
        
        IF(ionarray(i,2) .NE. iontype) THEN
           
           PRINT *, i,ionarray(i,1), ionarray(i,2)
           PRINT *, "Something wrong in ionarray type"
           STOP
           
        END IF
        
     END DO
     
     
     OPEN(unit = 93,file="iontypelist.txt",action="write",status="replace&
          &")
     
     WRITE(93,*) "Reference type/count: ", iontype, ioncnt
     
     DO i = 1,ioncnt
        WRITE(93,'(3(I0,1X))') i, ionarray(i,1), ionarray(i,2)
     END DO
     
     CLOSE(93)
     
     ! Load counterion array
     
     i = 0
     
     DO WHILE(dumcionarr(i+1,1) .NE. -1) 
        
        i = i + 1
        counterarray(i,1) = dumcionarr(i,1)
        counterarray(i,2) = dumcionarr(i,2)
        
     END DO
     
     IF(i .NE. c_ioncnt) THEN
        PRINT *, i, c_ioncnt
        STOP "Wrong total count in counterarray"
     END IF
     
     DO i = 1,c_ioncnt
        
        IF(counterarray(i,1) == -1 .OR. counterarray(i,2) == -1) THEN
           
           PRINT *, i,counterarray(i,1), counterarray(i,2)
           PRINT *, "Something wrong in assigning counterarray"
           STOP
           
        END IF
        
        IF(counterarray(i,2) .NE. c_iontype) THEN
           
           PRINT *, i,counterarray(i,1), counterarray(i,2)
           PRINT *, "Something wrong in counterionarray type"
           STOP
           
        END IF
        
     END DO
     
     
     OPEN(unit = 93,file="cntionlist.txt",action="write",status="repl&
          &ace")
     
     WRITE(93,*) "Reference type/count: ", c_iontype, c_ioncnt
     
     DO i = 1,c_ioncnt
        
        WRITE(93,'(3(I0,1X))') i, counterarray(i,1), counterarray(i,2)
        
     END DO
     
     CLOSE(93)
  
  ELSE

     ALLOCATE(ionarray(1,2),stat = AllocateStatus)
     DEALLOCATE(ionarray)
     ALLOCATE(counterarray(1,2),stat = AllocateStatus)
     DEALLOCATE(counterarray)

  END IF


  ! Polymerion array required only for diffusion systems

  IF (pion_dynflag) THEN

     PRINT *, "Number of atoms of polyion type: ",p_ioncnt

     ALLOCATE(polyionarray(p_ioncnt,2),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate polyionarray"

     i = 0

     DO WHILE(dumpionarr(i+1,1) .NE. -1) 

        i = i + 1
        polyionarray(i,1) = dumpionarr(i,1)
        polyionarray(i,2) = dumpionarr(i,2)
     
     END DO
     
     IF(i .NE. p_ioncnt) THEN
        PRINT *, i, p_ioncnt
        STOP "Wrong total count in counterarray"
     END IF
     
     DO i = 1,p_ioncnt
     
        IF(polyionarray(i,1) == -1 .OR. polyionarray(i,2) == -1) THEN
           
           PRINT *, i,polyionarray(i,1), polyionarray(i,2)
           PRINT *, "Something wrong in assigning polyionarray"
           STOP
           
        END IF
     
        IF(polyionarray(i,2) .NE. p_iontype) THEN
        
           PRINT *, i,polyionarray(i,1), polyionarray(i,2)
           PRINT *, "Something wrong in polyionarray"
           STOP
           
        END IF
     
     END DO
  
     
     OPEN(unit = 93,file="polyionlist.txt",action="write",status="repl&
          &ace")
  
     WRITE(93,*) "Reference type/count: ", p_iontype, p_ioncnt
     
     DO i = 1,p_ioncnt
        
        WRITE(93,'(3(I0,1X))') i,polyionarray(i,1), polyionarray(i,2)
        
     END DO
     
     CLOSE(93)

  ELSE

     ALLOCATE(polyionarray(1,2),stat = AllocateStatus)
     DEALLOCATE(polyionarray)

  END IF

END SUBROUTINE SORTALLARRAYS

!--------------------------------------------------------------------

SUBROUTINE MAP_REFTYPE(jin,atype,jout)
! Maps atomid into the corresponding place in array
  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER :: i
  INTEGER, INTENT(IN):: jin,atype
  INTEGER, INTENT(OUT) :: jout

  jout = -1

  IF(atype == iontype) THEN

     DO i = 1,ioncnt
        
        IF(jin == ionarray(i,1)) THEN

           jout = i

           EXIT

        END IF

     END DO

  ELSEIF(atype == c_iontype) THEN

     DO i = 1,c_ioncnt
        
        IF(jin == counterarray(i,1)) THEN

           jout = i

           EXIT

        END IF

     END DO

  ELSEIF(atype == p_iontype) THEN

     DO i = 1,p_ioncnt
        
        IF(jin == polyionarray(i,1)) THEN
           
           jout = i

           EXIT
           
        END IF
        
     END DO

  END IF
  
  IF(jout == -1) THEN
     
     PRINT *, jin, atype
     STOP "Could not find a match"

  END IF


END SUBROUTINE MAP_REFTYPE

!--------------------------------------------------------------------

SUBROUTINE DIFF_IONS()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,aid,atype,ierr,jout
  REAL    :: rxcm, rycm, rzcm
  REAL, DIMENSION(0:nframes-1) :: gxarr,gyarr,gzarr

  OPEN(unit = 91,file = "iondiffusivity.dat", status="replace",action&
       &="write",iostat = ierr)

  IF(ierr /= 0) STOP "Diffusion file not found"

  WRITE(91,*) ioncnt, iontype

  PRINT *, "Computing diffusivity of ions .. "

! Ion Diffusion Analysis

  DO i = 0,nframes-1

     gxarr(i) = 0.0
     gyarr(i) = 0.0
     gzarr(i) = 0.0

  END DO

!$OMP PARALLEL DO PRIVATE(tinc,ifin,i,tim,j,rxcm,rycm,rzcm,aid)&
!$OMP&  REDUCTION(+:gxarr,gyarr,gzarr)
  DO tinc = 0, nframes-1

     ifin = nframes - tinc

     DO i = 1,ifin

        tim = i + tinc

        DO j = 1,ioncnt

           rxcm = itrx_lmp(j,tim) - itrx_lmp(j,i)
           rycm = itry_lmp(j,tim) - itry_lmp(j,i)
           rzcm = itrz_lmp(j,tim) - itrz_lmp(j,i)

           gxarr(tinc) = gxarr(tinc) + rxcm**2
           gyarr(tinc) = gyarr(tinc) + rycm**2
           gzarr(tinc) = gzarr(tinc) + rzcm**2


        END DO

     END DO

     gxarr(tinc) = gxarr(tinc)/(ifin*ioncnt)
     gyarr(tinc) = gyarr(tinc)/(ifin*ioncnt)
     gzarr(tinc) = gzarr(tinc)/(ifin*ioncnt)


  END DO
!$OMP END PARALLEL DO

  DO i = 0, nframes-1

     WRITE(91,"(I10,1X,3(F14.5,1X))") i, gxarr(i),gyarr(i)&
          &,gzarr(i)

  END DO

  CLOSE(91)
END SUBROUTINE DIFF_IONS

!--------------------------------------------------------------------

SUBROUTINE DIFF_COUNTERIONS()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,tinc,ifin,tim,aid,atype,ierr,jout
  REAL    :: rxcm, rycm, rzcm
  REAL, DIMENSION(0:nframes-1) :: gxarr,gyarr,gzarr

  OPEN(unit = 91,file = "ciondiff.txt", status="replace"&
       &,action="write",iostat = ierr)

  IF(ierr /= 0) STOP "Diffusion file not found"

  WRITE(91,*) c_ioncnt, c_iontype

  PRINT *, "Computing diffusivity of counterions .. "

! Ion Diffusion Analysis

  DO i = 0,nframes-1

     gxarr(i) = 0.0
     gyarr(i) = 0.0
     gzarr(i) = 0.0

  END DO

! To do shifted time average for segmental diffusion

!$OMP PARALLEL DO PRIVATE(tinc,ifin,i,tim,j,rxcm,rycm,rzcm,aid)&
!$OMP&  REDUCTION(+:gxarr,gyarr,gzarr)
  DO tinc = 0, nframes-1

     ifin = nframes - tinc
     
     DO i = 1,ifin
        
        tim = i + tinc
      
        DO j = 1,c_ioncnt

           rxcm = ctrx_lmp(j,tim) - ctrx_lmp(j,i)
           rycm = ctry_lmp(j,tim) - ctry_lmp(j,i)
           rzcm = ctrz_lmp(j,tim) - ctrz_lmp(j,i)

           gxarr(tinc) = gxarr(tinc) + rxcm**2
           gyarr(tinc) = gyarr(tinc) + rycm**2
           gzarr(tinc) = gzarr(tinc) + rzcm**2
                      
        END DO

     END DO
     
     gxarr(tinc) = gxarr(tinc)/(REAL(ifin*c_ioncnt))
     gyarr(tinc) = gyarr(tinc)/(REAL(ifin*c_ioncnt))
     gzarr(tinc) = gzarr(tinc)/(REAL(ifin*c_ioncnt))

     
  END DO
!$OMP END PARALLEL DO

  DO i = 0, nframes-1

     WRITE(91,"(I10,1X,3(F14.5,1X))") i, gxarr(i),gyarr(i),gzarr(i)

  END DO

  CLOSE(91)

END SUBROUTINE DIFF_COUNTERIONS

!--------------------------------------------------------------------

SUBROUTINE CHECK_MOMENTUM(tval)

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tval
  INTEGER :: i,aid,atype
  REAL :: xmom, ymom, zmom
  REAL, PARAMETER :: tol_mom = 1e-5

  xmom = 0; ymom = 0; zmom = 0

  DO i = 1, ntotatoms
     
     aid = vel_xyz(i,1); atype = aidvals(aid,3)
     xmom = xmom + masses(atype,1)*vel_xyz(i,2)
     ymom = ymom + masses(atype,1)*vel_xyz(i,3)
     zmom = zmom + masses(atype,1)*vel_xyz(i,4)

  END DO

  IF( (abs(xmom) .GT. tol_mom) .OR. (abs(ymom) .GT. tol_mom) .OR.&
       & (abs(zmom) .GT. tol_mom) ) THEN

     PRINT *, "WARNING: Net momentum not zero: ", tval,xmom,ymom,zmom

  ELSE

     IF(tval == 0) THEN

        PRINT *, "Momentum conserved at the beginning: ", xmom, ymom,&
             & zmom

     END IF

  END IF

END SUBROUTINE CHECK_MOMENTUM

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_TOPO_ARRAYS()

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate LAMMPS structure

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(charge_lmp(ntotatoms,1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate charge_lmp"
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

! Allocate box details

  ALLOCATE(boxx_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxx_arr"
  ALLOCATE(boxy_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxy_arr"
  ALLOCATE(boxz_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxz_arr"

  PRINT *, "Successfully allocated memory for topology"

END SUBROUTINE ALLOCATE_TOPO_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS()

  USE ANALYZE_PARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate for statics

  IF(rdfcalc) THEN
     ALLOCATE(rdfarray(0:rmaxbin-1,npairs),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdfarray"
  ELSE
     ALLOCATE(rdfarray(1,1),stat = AllocateStatus)
     DEALLOCATE(rdfarray)
  END IF

! Allocate for dynamics 

  IF(ion_dynflag) THEN
     ALLOCATE(itrx_lmp(ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itrx_lmp"
     ALLOCATE(itry_lmp(ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itry_lmp"
     ALLOCATE(itrz_lmp(ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate itrz_lmp"
  ELSE
     ALLOCATE(itrx_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(itry_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(itrz_lmp(1,nframes),stat = AllocateStatus)
     DEALLOCATE(itrx_lmp)
     DEALLOCATE(itry_lmp)
     DEALLOCATE(itrz_lmp)
  END IF

  IF(cion_dynflag) THEN
     ALLOCATE(ctrx_lmp(c_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctrx_lmp"
     ALLOCATE(ctry_lmp(c_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctry_lmp"
     ALLOCATE(ctrz_lmp(c_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ctrz_lmp"
  ELSE
     ALLOCATE(ctrx_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(ctry_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(ctrz_lmp(1,nframes),stat = AllocateStatus)
     DEALLOCATE(ctrx_lmp)
     DEALLOCATE(ctry_lmp)
     DEALLOCATE(ctrz_lmp)
  END IF

  IF(pion_dynflag) THEN
     ALLOCATE(ptrx_lmp(p_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ptrx_lmp"
     ALLOCATE(ptry_lmp(p_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ptry_lmp"
     ALLOCATE(ptrz_lmp(p_ioncnt,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ptrz_lmp"
  ELSE
     ALLOCATE(ptrx_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(ptry_lmp(1,nframes),stat = AllocateStatus)
     ALLOCATE(ptrz_lmp(1,nframes),stat = AllocateStatus)
     DEALLOCATE(ptrx_lmp)
     DEALLOCATE(ptry_lmp)
     DEALLOCATE(ptrz_lmp)
  END IF
  
  IF(catan_neighcalc) THEN
     ALLOCATE(cat_an_neighavg(1:maxneighsize),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate cat_an_neighavg"
     ALLOCATE(an_cat_neighavg(1:maxneighsize),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate an_cat_neighavg"
  ELSE
     ALLOCATE(cat_an_neighavg(1),stat = AllocateStatus)
     ALLOCATE(an_cat_neighavg(1),stat = AllocateStatus)
     DEALLOCATE(cat_an_neighavg)
     DEALLOCATE(an_cat_neighavg)
  END IF

  PRINT *, "Successfully allocated memory for analyis"

END SUBROUTINE ALLOCATE_ANALYSIS_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE ANALYZE_PARAMS

  IMPLICIT NONE

  !Global arrays
  DEALLOCATE(aidvals)
  DEALLOCATE(rxyz_lmp)
  DEALLOCATE(vel_xyz)
  DEALLOCATE(boxx_arr)
  DEALLOCATE(boxy_arr)
  DEALLOCATE(boxz_arr)

  !Topo arrays
  IF(ntotbonds /= 0) DEALLOCATE(bond_lmp)
  IF(ntotangls /= 0) DEALLOCATE(angl_lmp)
  IF(ntotdihds /= 0) DEALLOCATE(dihd_lmp)
  IF(ntotimprs /= 0) DEALLOCATE(impr_lmp)
  
  !Statics calculations arrays
  IF(rdfcalc) DEALLOCATE(rdfarray)

  !Dynamic calculations arrays
  IF(ion_dynflag) THEN
     DEALLOCATE(itrx_lmp)
     DEALLOCATE(itry_lmp)
     DEALLOCATE(itrz_lmp)
  END IF

  IF(cion_dynflag) THEN
     DEALLOCATE(ctrx_lmp)
     DEALLOCATE(ctry_lmp)
     DEALLOCATE(ctrz_lmp)
  END IF

  IF(pion_dynflag) THEN
     DEALLOCATE(ptrx_lmp)
     DEALLOCATE(ptry_lmp)
     DEALLOCATE(ptrz_lmp)
  END IF

END SUBROUTINE DEALLOCATE_ARRAYS

!--------------------------------------------------------------------

