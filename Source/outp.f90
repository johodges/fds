#ifndef GITHASH_PP
#define GITHASH_PP "unknown"
#endif

!> \brief Routines for handling output

MODULE GET_DATA

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE OUTPUT_CLOCKS
USE MESH_POINTERS
USE DEVICE_VARIABLES
USE CONTROL_VARIABLES

USE OUTPUT_DATA
USE PROPERTY_DATA
!USE MESH_VARIABLES
USE COMPLEX_GEOMETRY, ONLY : WRITE_GEOM,WRITE_GEOM_ALL,CC_FGSC,CC_IDCF,CC_IDCC,CC_UNKZ,CC_UNKF,CC_FTYPE_RCGAS,&
                             CC_FTYPE_CFGAS,CC_FTYPE_CFINB,CC_SOLID,CC_CGSC,CC_IDRC,CC_CUTCFE,TRIANGULATE,&
                             CC_VGSC,CC_GASPHASE,MAKE_UNIQUE_VERT_ARRAY,AVERAGE_FACE_VALUES
USE CC_SCALARS, ONLY : GET_PRES_CFACE,GET_PRES_CFACE_TEST,GET_UVWGAS_CFACE,GET_MUDNS_CFACE
IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE
REAL(EB), POINTER, DIMENSION(:,:,:) :: WFX,WFY,WFZ
LOGICAL :: EX,DRY,OPN,FROM_BNDF=.FALSE.

TYPE(GEOMETRY_TYPE), POINTER :: G
TYPE (MESH_TYPE), POINTER :: M
TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (OBSTRUCTION_TYPE), POINTER :: OB
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (SPECIES_TYPE), POINTER :: SS
TYPE (REACTION_TYPE), POINTER :: RN
TYPE (SURFACE_TYPE),POINTER :: SF
TYPE (MATERIAL_TYPE),POINTER :: ML
TYPE (PROPERTY_TYPE), POINTER :: PY
TYPE (DEVICE_TYPE), POINTER :: DV, DV2
TYPE (SUBDEVICE_TYPE), POINTER :: SDV
TYPE (SLICE_TYPE), POINTER :: SL
TYPE (WALL_TYPE), POINTER :: WC
TYPE (THIN_WALL_TYPE), POINTER :: TW
TYPE (CFACE_TYPE), POINTER :: CFA
TYPE (BOUNDARY_FILE_TYPE), POINTER :: BF
TYPE (ISOSURFACE_FILE_TYPE), POINTER :: IS
TYPE (INITIALIZATION_TYPE), POINTER :: IN

PUBLIC GET_SMOKE3D_QQ,GET_BNDF_PACK,GET_GEOMVALS,GET_GEOMSIZES,GET_GEOMINFO,&
       GAS_PHASE_OUTPUT,SOLID_PHASE_OUTPUT

CONTAINS

SUBROUTINE GET_GEOMVALS(CC_INTERP2FACES,CC_CELL_CENTERED,SLICETYPE,&
                        I1,I2,J1,J2,K1,K2,NFACES,NFACES_CUTCELLS,VALS,&
                        IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,&
                        PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM,OPT_BNDF_INDEX)

USE PHYSICAL_FUNCTIONS, ONLY: GET_MASS_FRACTION

! copy data from QQ array into VALS(1:NFACES)

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: I1,I2,J1,J2,K1,K2,NFACES,NFACES_CUTCELLS,&
                       IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,NM
INTEGER, OPTIONAL,INTENT(IN) :: OPT_BNDF_INDEX
CHARACTER(*), INTENT(IN) :: SLICETYPE
LOGICAL, INTENT(IN) :: CC_INTERP2FACES,CC_CELL_CENTERED
REAL(FB), INTENT(OUT), DIMENSION(NFACES) :: VALS

INTEGER :: DIR, SLICE, IFACE
INTEGER :: I,J,K
CHARACTER(LEN=100) :: SLICETYPE_LOCAL
INTEGER :: CELLTYPE
INTEGER :: ICF, NVF, IFACECF, IVCF, IFACECUT

INTEGER :: X1AXIS, II, JJ, KK, ICC, JCC, NFC, ICCF, ICF2, IFACE2
REAL(EB):: VAL_CF

LOGICAL :: IS_RCFACE

SLICETYPE_LOCAL=TRIM(SLICETYPE) ! only generate CUTCELLS slice files if the immersed geometry option is turned on
IF (SLICETYPE=='INCLUDE_GEOM' .AND. .NOT.CC_IBM) SLICETYPE_LOCAL='IGNORE_GEOM'

CALL GETSLICEDIR(I1,I2,J1,J2,K1,K2,DIR,SLICE)
IF (SLICETYPE_LOCAL=='IGNORE_GEOM') THEN
   IFACE = 0
   IF (DIR==1) THEN
      DO K = K1+1, K2
         DO J = J1+1, J2
            IFACE = IFACE + 1
            VALS(IFACE) = QQ(SLICE,J,K,1)

            IFACE = IFACE + 1
            VALS(IFACE) = QQ(SLICE,J,K,1)
         ENDDO
      ENDDO
   ELSE IF (DIR==2) THEN
      DO K = K1+1, K2
         DO I = I1+1, I2
            IFACE = IFACE + 1
            VALS(IFACE) = QQ(I,SLICE,K,1)

            IFACE = IFACE + 1
            VALS(IFACE) = QQ(I,SLICE,K,1)
         ENDDO
      ENDDO
   ELSE
      DO J = J1+1, J2
         DO I = I1+1, I2
            IFACE = IFACE + 1
            VALS(IFACE) = QQ(I,J,SLICE,1)

            IFACE = IFACE + 1
            VALS(IFACE) = QQ(I,J,SLICE,1)
         ENDDO
      ENDDO
   ENDIF
ELSEIF (SLICETYPE_LOCAL=='INCLUDE_GEOM') THEN ! INTERP_C2F_FIELD
   X1AXIS = DIR
   IFACE = 0
   IFACECUT=NFACES-NFACES_CUTCELLS  ! start cutcell counter after 'regular' cells
   IF (DIR==1) THEN
      DO K = K1+1, K2
         DO J = J1+1, J2
            IF (ANY(CELL(CELL_INDEX(SLICE:SLICE+1,J,K))%SOLID)) CYCLE
            CELLTYPE = FCVAR(SLICE,J,K,CC_FGSC,IAXIS)
            IF (CELLTYPE == CC_CUTCFE) THEN
               ICF = FCVAR(SLICE,J,K,CC_IDCF,IAXIS) ! is a cut cell
               DO IFACECF=1,CUT_FACE(ICF)%NFACE
                  CALL GET_GASCUTFACE_SCALAR_SLICE(VAL_CF,X1AXIS,ICF,IFACECF,CC_INTERP2FACES,CC_CELL_CENTERED,&
                         IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
                  NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                  DO IVCF = 1, NVF-2
                     IFACECUT = IFACECUT + 1
                     VALS(IFACECUT) = REAL(VAL_CF,FB)
                  ENDDO
               ENDDO
               CALL GET_SOLIDCUTFACE_SCALAR_SLICE(X1AXIS,ICF,VAL_CF,&
                  IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
               DO IFACECF=CUT_FACE(ICF)%NFACE+1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE
                  NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                  DO IVCF = 1, NVF-2
                     IFACECUT = IFACECUT + 1
                     VALS(IFACECUT) = REAL(VAL_CF,FB)
                  ENDDO
               ENDDO
            ELSEIF(CELLTYPE == CC_SOLID) THEN
               CALL GET_SOLIDREGFACE_SCALAR_SLICE(X1AXIS,SLICE,J,K,VAL_CF,&
                  IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
               IFACE = IFACE + 1  ! is a solid or gas cell
               VALS(IFACE) = REAL(VAL_CF,FB)

               IFACE = IFACE + 1
               VALS(IFACE) = REAL(VAL_CF,FB)
            ELSE
               ! Check if FACE is TYPE RC face:
               IS_RCFACE = (CCVAR(SLICE,J,K,CC_CGSC)==CC_CUTCFE) .OR. (CCVAR(SLICE+1,J,K,CC_CGSC)==CC_CUTCFE)
               IF (IS_RCFACE) THEN
                  ! TO DO: Place holder to interpolate Slice Variable to RCFACE:
                  ! ..
                  IFACE = IFACE + 1  ! is a gas cell
                  VALS(IFACE) = QQ(SLICE,J,K,1)

                  IFACE = IFACE + 1
                  VALS(IFACE) = QQ(SLICE,J,K,1)

               ELSE
                  IFACE = IFACE + 1  ! is a gas cell
                  VALS(IFACE) = QQ(SLICE,J,K,1)

                  IFACE = IFACE + 1
                  VALS(IFACE) = QQ(SLICE,J,K,1)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ELSEIF (DIR==2) THEN
      DO K = K1+1, K2
         DO I = I1+1, I2
            IF (ANY(CELL(CELL_INDEX(I,SLICE:SLICE+1,K))%SOLID)) CYCLE
            CELLTYPE = FCVAR(I,SLICE,K,CC_FGSC,JAXIS)
            IF (CELLTYPE == CC_CUTCFE) THEN
               ICF = FCVAR(I,SLICE,K,CC_IDCF,JAXIS)
               DO IFACECF=1,CUT_FACE(ICF)%NFACE
                  CALL GET_GASCUTFACE_SCALAR_SLICE(VAL_CF,X1AXIS,ICF,IFACECF,CC_INTERP2FACES,CC_CELL_CENTERED,&
                         IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
                  NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                  DO IVCF = 1, NVF-2
                     IFACECUT = IFACECUT + 1
                     VALS(IFACECUT) = REAL(VAL_CF,FB)
                  ENDDO
               ENDDO
               CALL GET_SOLIDCUTFACE_SCALAR_SLICE(X1AXIS,ICF,VAL_CF,&
                  IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
               DO IFACECF=CUT_FACE(ICF)%NFACE+1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE
                  NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                  DO IVCF = 1, NVF-2
                     IFACECUT = IFACECUT + 1
                     VALS(IFACECUT) = REAL(VAL_CF,FB)
                  ENDDO
               ENDDO
            ELSEIF(CELLTYPE == CC_SOLID) THEN
               CALL GET_SOLIDREGFACE_SCALAR_SLICE(X1AXIS,I,SLICE,K,VAL_CF,&
                  IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
               IFACE = IFACE + 1  ! is a solid or gas cell
               VALS(IFACE) = REAL(VAL_CF,FB)

               IFACE = IFACE + 1
               VALS(IFACE) = REAL(VAL_CF,FB)
            ELSE
               ! Check if FACE is TYPE RC face:
               IS_RCFACE = (CCVAR(I,SLICE,K,CC_CGSC)==CC_CUTCFE) .OR. (CCVAR(I,SLICE+1,K,CC_CGSC)==CC_CUTCFE)
               IF (IS_RCFACE) THEN
                  ! TO DO: Place holder to interpolate Slice Variable to RCFACE:
                  ! ..
                  IFACE = IFACE + 1  ! is a gas cell
                  VALS(IFACE) = QQ(I,SLICE,K,1)

                  IFACE = IFACE + 1
                  VALS(IFACE) = QQ(I,SLICE,K,1)
               ELSE
                  IFACE = IFACE + 1
                  VALS(IFACE) = QQ(I,SLICE,K,1)

                  IFACE = IFACE + 1
                  VALS(IFACE) = QQ(I,SLICE,K,1)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ELSE
      DO J = J1+1, J2
         DO I = I1+1, I2
            IF (ANY(CELL(CELL_INDEX(I,J,SLICE:SLICE+1))%SOLID)) CYCLE
            CELLTYPE = FCVAR(I,J,SLICE,CC_FGSC,KAXIS)
            IF (CELLTYPE == CC_CUTCFE) THEN
               ICF = FCVAR(I,J,SLICE,CC_IDCF,KAXIS)
               DO IFACECF=1,CUT_FACE(ICF)%NFACE
                  CALL GET_GASCUTFACE_SCALAR_SLICE(VAL_CF,X1AXIS,ICF,IFACECF,CC_INTERP2FACES,CC_CELL_CENTERED,&
                         IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
                  NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                  DO IVCF = 1, NVF-2
                     IFACECUT = IFACECUT + 1
                     VALS(IFACECUT) = REAL(VAL_CF,FB)
                  ENDDO
               ENDDO
               CALL GET_SOLIDCUTFACE_SCALAR_SLICE(X1AXIS,ICF,VAL_CF,&
                  IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
               DO IFACECF=CUT_FACE(ICF)%NFACE+1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE
                  NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                  DO IVCF = 1, NVF-2
                     IFACECUT = IFACECUT + 1
                     VALS(IFACECUT) = REAL(VAL_CF,FB)
                  ENDDO
               ENDDO
            ELSEIF(CELLTYPE == CC_SOLID) THEN
               CALL GET_SOLIDREGFACE_SCALAR_SLICE(X1AXIS,I,J,SLICE,VAL_CF,&
                  IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)
               IFACE = IFACE + 1  ! is a solid or gas cell
               VALS(IFACE) = REAL(VAL_CF,FB)

               IFACE = IFACE + 1
               VALS(IFACE) = REAL(VAL_CF,FB)
            ELSE
               ! Check if FACE is TYPE RC face:
               IS_RCFACE = (CCVAR(I,J,SLICE,CC_CGSC)==CC_CUTCFE) .OR. (CCVAR(I,J,SLICE+1,CC_CGSC)==CC_CUTCFE)
               IF (IS_RCFACE) THEN
                  ! TO DO: Place holder to interpolate Slice Variable to RCFACE:
                  ! ..
                  IFACE = IFACE + 1  ! is a gas cell
                  VALS(IFACE) = QQ(I,J,SLICE,1)

                  IFACE = IFACE + 1
                  VALS(IFACE) = QQ(I,J,SLICE,1)
               ELSE
                  IFACE = IFACE + 1
                  VALS(IFACE) = QQ(I,J,SLICE,1)

                  IFACE = IFACE + 1
                  VALS(IFACE) = QQ(I,J,SLICE,1)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDIF
ELSEIF (SLICETYPE_LOCAL=='INBOUND_FACES') THEN
   IFACECUT=NFACES-NFACES_CUTCELLS  ! start cutcell counter after 'regular' cells
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
         IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
         IF (CCVAR(I,J,K,CC_IDCF) > 0) THEN
            ICF = CCVAR(I,J,K,CC_IDCF)
            DO IFACECF=1,CUT_FACE(ICF)%NFACE
               VAL_CF = SOLID_PHASE_OUTPUT(ABS(IND),T,NM,Y_INDEX,Z_INDEX,PART_INDEX,OPT_BNDF_INDEX=OPT_BNDF_INDEX, &
                                           OPT_CFACE_INDEX=CUT_FACE(ICF)%CFACE_INDEX(IFACECF))
               NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
               DO IVCF = 1, NVF-2 ! face is convex
                  IFACECUT = IFACECUT + 1
                  VALS(IFACECUT) = REAL(VAL_CF,FB)
               ENDDO
            ENDDO
         ENDIF
         ENDDO
      ENDDO
   ENDDO
ELSEIF (SLICETYPE_LOCAL=='CUT_CELLS') THEN
   IFACECUT=NFACES-NFACES_CUTCELLS
   VAL_CF=0._EB
   DO KK = 1, KBAR
      DO JJ = 1, JBAR
         DO II = 1, IBAR
            IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID) CYCLE
            IF (CCVAR(II,JJ,KK,CC_IDCC) <= 0) CYCLE
            ICC = CCVAR(II,JJ,KK,CC_IDCC)
            DO JCC=1,CUT_CELL(ICC)%NCELL
               NFC=CUT_CELL(ICC)%CCELEM(1,JCC)
               ! Loop on faces corresponding to cut-cell ICC2:
               DO ICCF=1,NFC
                  IFACE=CUT_CELL(ICC)%CCELEM(ICCF+1,JCC)
                  SELECT CASE(CUT_CELL(ICC)%FACE_LIST(1,IFACE))
                     CASE(CC_FTYPE_RCGAS) ! REGULAR GASPHASE
                        DO IVCF = 1,2
                           IFACECUT = IFACECUT + 1
                           VALS(IFACECUT) = REAL(VAL_CF,FB)
                        ENDDO

                     CASE(CC_FTYPE_CFGAS)
                        ICF2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
                        IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
                        NVF     = CUT_FACE(ICF2)%CFELEM(1,IFACE2)
                        DO IVCF = 1, NVF-2 ! for now assume face is convex
                           IFACECUT = IFACECUT + 1
                           VALS(IFACECUT) = REAL(VAL_CF,FB)
                        ENDDO

                     CASE(CC_FTYPE_CFINB)
                        ICF2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
                        IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
                        NVF     = CUT_FACE(ICF2)%CFELEM(1,IFACE2); DIR = 0
                        DO IVCF = 1, NVF-2 ! face is convex
                           IFACECUT = IFACECUT + 1
                           VALS(IFACECUT) = REAL(VAL_CF,FB)
                        ENDDO

                  END SELECT
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE GET_GEOMVALS


SUBROUTINE GET_GEOMSIZES(SLICETYPE,I1,I2,J1,J2,K1,K2,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS)

! determine NVERTS and NFACES for one of the following cases
!
! IGNORE_GEOM  - creates a slice file geometry file that ignores immersed geometric objects .  Triangles inside obstacle
!                regions (a solid) are tagged with a 1, triangles outside of obstacle regions (the gas) are tagged
!                with a 0 . Smokeview uses this information to show/hide these two regions
! INCLUDE_GEOM - creates a slice file geometry file that accounts for immersed geometric objects .  If there are no immersed
!                objects present then this slice type is equivalent to the 'IGNORE_GEOM' case.  Triangles completely inside a
!                solid are tagged with a 1, triangles completely in the gas are tagged with a 0 and triangles in a cutcell are
!                with a tagged 2.  As with the IGNORE_GEOM type, Smokeview uses this information to show/hide these regions

   CHARACTER(*), INTENT(IN) :: SLICETYPE
   INTEGER, INTENT(IN) :: I1,I2,J1,J2,K1,K2
   INTEGER, INTENT(OUT) :: NVERTS, NVERTS_CUTCELLS, NFACES, NFACES_CUTCELLS

   INTEGER :: DIR,SLICE
   INTEGER :: I, J, K
   INTEGER :: ICF, IFACE, NVF, ICC, JCC, ICF2, IFACE2, NFC, ICCF

   CHARACTER(LEN=100) :: SLICETYPE_LOCAL

   SLICETYPE_LOCAL=TRIM(SLICETYPE) ! only generate CUTCELLS slice files if the immersed geometry option is turned on
   IF (SLICETYPE=='INCLUDE_GEOM' .AND. .NOT.CC_IBM) SLICETYPE_LOCAL='IGNORE_GEOM'

   NVERTS=0
   NFACES=0
   NVERTS_CUTCELLS=0
   NFACES_CUTCELLS=0
   IF (SLICETYPE_LOCAL=='IGNORE_GEOM') THEN
      CALL GETSLICEDIR(I1,I2,J1,J2,K1,K2,DIR,SLICE)
      IF (DIR==1) THEN
        NVERTS = (J2 + 1 - J1)*(K2 + 1 - K1)
        NFACES = 2*(J2 - J1)*(K2 - K1)
      ELSE IF (DIR==2) THEN
        NVERTS = (I2 + 1 - I1)*(K2 + 1 - K1)
        NFACES = 2*(I2 - I1)*(K2 - K1)
      ELSE
        NVERTS = (I2 + 1 - I1)*(J2 + 1 - J1)
        NFACES = 2*(I2 - I1)*(J2 - J1)
      ENDIF
   ELSE IF (SLICETYPE_LOCAL=='INCLUDE_GEOM') THEN
      CALL GETSLICEDIR(I1,I2,J1,J2,K1,K2,DIR,SLICE)
      IF (DIR==1) THEN
         NVERTS = (J2 + 1 - J1)*(K2 + 1 - K1)
         NFACES = 0
         DO K = K1+1, K2
            DO J = J1+1, J2
               IF (ANY(CELL(CELL_INDEX(SLICE:SLICE+1,J,K))%SOLID)) CYCLE
               IF (FCVAR(SLICE,J,K,CC_FGSC,IAXIS) == CC_CUTCFE) THEN
                  ICF = FCVAR(SLICE,J,K,CC_IDCF,IAXIS) ! a cutcell so count number of faces
                  DO IFACE=1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE ! Adds also SOLID side faces.
                     NVF=CUT_FACE(ICF)%CFELEM(1,IFACE)
                     NFACES_CUTCELLS = NFACES_CUTCELLS + NVF - 2
                     NVERTS_CUTCELLS = NVERTS_CUTCELLS + NVF
                  ENDDO
               ELSE
                  NFACES = NFACES + 2 ! a gas or solid cell so add 2 to the number of faces
               ENDIF
            ENDDO
         ENDDO
      ELSE IF (DIR==2) THEN
         NVERTS = (I2 + 1 - I1)*(K2 + 1 - K1)
         DO K = K1+1, K2
            DO I = I1+1, I2
               IF(ANY(CELL(CELL_INDEX(I,SLICE:SLICE+1,K))%SOLID)) CYCLE
               IF (FCVAR(I,SLICE,K,CC_FGSC,JAXIS) == CC_CUTCFE) THEN
                  ICF = FCVAR(I,SLICE,K,CC_IDCF,JAXIS)
                  DO IFACE=1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE ! Adds also SOLID side faces.
                     NVF=CUT_FACE(ICF)%CFELEM(1,IFACE)
                     NFACES_CUTCELLS = NFACES_CUTCELLS + NVF - 2
                     NVERTS_CUTCELLS = NVERTS_CUTCELLS + NVF
                  ENDDO
               ELSE
                  NFACES = NFACES + 2
               ENDIF
            ENDDO
         ENDDO
      ELSE
         NVERTS = (I2 + 1 - I1)*(J2 + 1 - J1)
         DO I = I1+1, I2
            DO J = J1+1, J2
               IF(ANY(CELL(CELL_INDEX(I,J,SLICE:SLICE+1))%SOLID)) CYCLE
               IF (FCVAR(I,J,SLICE,CC_FGSC,KAXIS) == CC_CUTCFE) THEN
                  ICF = FCVAR(I,J,SLICE,CC_IDCF,KAXIS)
                  DO IFACE=1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE ! Adds also SOLID side faces.
                     NVF=CUT_FACE(ICF)%CFELEM(1,IFACE)
                     NFACES_CUTCELLS = NFACES_CUTCELLS + NVF - 2
                     NVERTS_CUTCELLS = NVERTS_CUTCELLS + NVF
                  ENDDO
               ELSE
                  NFACES = NFACES + 2
               ENDIF
            ENDDO
         ENDDO
      ENDIF
   ELSE IF (SLICETYPE_LOCAL=='INBOUND_FACES') THEN
      DO K = 1, KBAR
         DO J = 1, JBAR
            DO I = 1, IBAR
               IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
               IF (CCVAR(I,J,K,CC_IDCF) > 0) THEN ! There are INBOUNDARY cut-faces on this cell:
                  ICF = CCVAR(I,J,K,CC_IDCF)
                  DO IFACE=1,CUT_FACE(ICF)%NFACE ! Adds also SOLID side faces.
                     NVF=CUT_FACE(ICF)%CFELEM(1,IFACE)
                     NFACES_CUTCELLS = NFACES_CUTCELLS + NVF - 2
                     NVERTS_CUTCELLS = NVERTS_CUTCELLS + NVF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ELSE IF (SLICETYPE_LOCAL=='CUT_CELLS') THEN
      DO K = 1, KBAR
         DO J = 1, JBAR
            DO I = 1, IBAR
               IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
               IF (CCVAR(I,J,K,CC_IDCC) <= 0) CYCLE
               ICC = CCVAR(I,J,K,CC_IDCC)
               DO JCC=1,CUT_CELL(ICC)%NCELL
                  NFC=CUT_CELL(ICC)%CCELEM(1,JCC)
                  ! Loop on faces corresponding to cut-cell ICC2:
                  DO ICCF=1,NFC
                     IFACE=CUT_CELL(ICC)%CCELEM(ICCF+1,JCC)
                     SELECT CASE(CUT_CELL(ICC)%FACE_LIST(1,IFACE))
                     CASE(CC_FTYPE_RCGAS) ! REGULAR GASPHASE
                        NVF = 4
                        NFACES_CUTCELLS = NFACES_CUTCELLS + NVF - 2
                        NVERTS_CUTCELLS = NVERTS_CUTCELLS + NVF
                     CASE(CC_FTYPE_CFGAS)
                        ICF2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
                        IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
                        NVF=CUT_FACE(ICF2)%CFELEM(1,IFACE2)
                        NFACES_CUTCELLS = NFACES_CUTCELLS + NVF - 2
                        NVERTS_CUTCELLS = NVERTS_CUTCELLS + NVF
                     CASE(CC_FTYPE_CFINB)
                        ICF2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
                        IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
                        NVF=CUT_FACE(ICF2)%CFELEM(1,IFACE2)
                        NFACES_CUTCELLS = NFACES_CUTCELLS + NVF - 2
                        NVERTS_CUTCELLS = NVERTS_CUTCELLS + NVF
                     END SELECT
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   NFACES = NFACES + NFACES_CUTCELLS
   NVERTS = NVERTS + NVERTS_CUTCELLS
END SUBROUTINE GET_GEOMSIZES


SUBROUTINE GET_GEOMINFO(SLICETYPE,I1,I2,J1,J2,K1,K2,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS,&
                        VERTS,FACES,LOCATIONS,SURFIND,GEOMIND)

! generate VERTS(1:3*NVERTS) and FACES(1:3*NFACES) arrays

   CHARACTER(*), INTENT(IN) :: SLICETYPE
   INTEGER, INTENT(IN) :: I1,I2,J1,J2,K1,K2
   INTEGER, INTENT(IN) :: NVERTS, NVERTS_CUTCELLS, NFACES, NFACES_CUTCELLS
   INTEGER, INTENT(OUT), DIMENSION(3*NFACES), TARGET :: FACES
   INTEGER, INTENT(OUT), DIMENSION(NFACES) :: LOCATIONS
   INTEGER, OPTIONAL, INTENT(OUT), DIMENSION(NFACES) :: SURFIND,GEOMIND
   REAL(FB), INTENT(OUT), DIMENSION(3*NVERTS), TARGET :: VERTS

   INTEGER :: VERT_OFFSET
   INTEGER, POINTER, DIMENSION(:) :: FACEPTR
   REAL(FB), POINTER, DIMENSION(:) :: VERTPTR

   INTEGER :: DIR, SLICE
   INTEGER :: NI, NJ, NK
   INTEGER :: I, J, K
   INTEGER IFACE, IVERT, IVERTCUT, IFACECUT, IVERTCF, IFACECF
   INTEGER VERTBEG, VERTEND, FACEBEG, FACEEND
   LOGICAL IS_SOLID
   INTEGER :: ICF, NVF, IVCF, IADD, JADD, KADD, X1AXIS
   INTEGER :: II, JJ, KK, ICC, JCC, NFC, ICCF, LOWHIGH, ILH, ICF2, IFACE2
   INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCTYPE

   CHARACTER(LEN=100) :: SLICETYPE_LOCAL

   SLICETYPE_LOCAL=TRIM(SLICETYPE) ! only generate CUTCELLS slice files if the immersed geometry option is turned on
   IF (SLICETYPE=='INCLUDE_GEOM' .AND. .NOT.CC_IBM) SLICETYPE_LOCAL='IGNORE_GEOM'

   LOCATIONS = 0 ! initially assume triangles are in gas and tag with 0
   IF (SLICETYPE_LOCAL=='IGNORE_GEOM') THEN
      NI = I2 + 1 - I1
      NJ = J2 + 1 - J1
      NK = K2 + 1 - K1
      CALL GETSLICEDIR(I1,I2,J1,J2,K1,K2,DIR,SLICE)
      IVERT = 0
      IFACE = 0
      IF (DIR==1) THEN
         DO K=K1,K2
            DO J=J1,J2
               DO I = SLICE,SLICE
                  IVERT = IVERT + 1
                  VERTS(3*IVERT-2) = XPLT(SLICE)
                  VERTS(3*IVERT-1) = YPLT(J)
                  VERTS(3*IVERT)   = ZPLT(K)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,NK-1
            DO J=1,NJ-1
               IS_SOLID = CELL(CELL_INDEX(SLICE,J+J1,K+K1))%SOLID
               IFACE = IFACE + 1
               IF (IS_SOLID) LOCATIONS(IFACE) = 1 + 16 ! triangle is in a solid so tag with 1
               FACES(3*IFACE-2) = IJK(  J,  K,NJ)
               FACES(3*IFACE-1) = IJK(J+1,  K,NJ)
               FACES(3*IFACE)   = IJK(J+1,K+1,NJ)

               IFACE = IFACE + 1
               IF (IS_SOLID) LOCATIONS(IFACE) = 1 + 4 ! triangle is in a solid so tag with 1
               FACES(3*IFACE-2) = IJK(  J,  K,NJ)
               FACES(3*IFACE-1) = IJK(J+1,K+1,NJ)
               FACES(3*IFACE)   = IJK(  J,K+1,NJ)
            ENDDO
         ENDDO
      ELSE IF (DIR==2) THEN
         DO K=K1,K2
            DO J=SLICE,SLICE
               DO I = I1,I2
                  IVERT = IVERT + 1
                  VERTS(3*IVERT-2) = XPLT(I)
                  VERTS(3*IVERT-1) = YPLT(SLICE)
                  VERTS(3*IVERT)   = ZPLT(K)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,NK-1
            DO I=1,NI-1
               IS_SOLID = CELL(CELL_INDEX(I+I1,SLICE,K+K1))%SOLID
               IFACE = IFACE + 1
               IF (IS_SOLID) LOCATIONS(IFACE) = 1 + 16
               FACES(3*IFACE-2) = IJK(  I,  K,NI)
               FACES(3*IFACE-1) = IJK(I+1,  K,NI)
               FACES(3*IFACE)   = IJK(I+1,K+1,NI)

               IFACE = IFACE + 1
               IF (IS_SOLID) LOCATIONS(IFACE) = 1 + 4
               FACES(3*IFACE-2) = IJK(  I,  K,NI)
               FACES(3*IFACE-1) = IJK(I+1,K+1,NI)
               FACES(3*IFACE)   = IJK(  I,K+1,NI)
            ENDDO
         ENDDO
      ELSE
         DO K=SLICE,SLICE
            DO J=J1,J2
               DO I = I1,I2
                  IVERT = IVERT + 1
                  VERTS(3*IVERT-2) = XPLT(I)
                  VERTS(3*IVERT-1) = YPLT(J)
                  VERTS(3*IVERT)   = ZPLT(SLICE)
               ENDDO
            ENDDO
         ENDDO
         DO J=1,NJ-1
            DO I=1,NI-1
               IS_SOLID = CELL(CELL_INDEX(I+I1,J+J1,SLICE))%SOLID
               IFACE = IFACE + 1
               IF (IS_SOLID) LOCATIONS(IFACE) = 1 + 16
               FACES(3*IFACE-2) = IJK(  I,  J,NI)
               FACES(3*IFACE-1) = IJK(I+1,  J,NI)
               FACES(3*IFACE)   = IJK(I+1,J+1,NI)

               IFACE = IFACE + 1
               IF (IS_SOLID) LOCATIONS(IFACE) = 1 + 4
               FACES(3*IFACE-2) = IJK(  I,  J,NI)
               FACES(3*IFACE-1) = IJK(I+1,J+1,NI)
               FACES(3*IFACE)   = IJK(  I,J+1,NI)
            ENDDO
         ENDDO
      ENDIF
   ELSE IF (SLICETYPE_LOCAL=='INCLUDE_GEOM') THEN
      IVERTCUT=NVERTS-NVERTS_CUTCELLS ! start cutcell counters after 'regular' cells
      IFACECUT=NFACES-NFACES_CUTCELLS
      NI = I2 + 1 - I1
      NJ = J2 + 1 - J1
      NK = K2 + 1 - K1
      CALL GETSLICEDIR(I1,I2,J1,J2,K1,K2,DIR,SLICE)
      IVERT = 0
      IFACE = 0
      IF (DIR==1) THEN
         DO K=K1,K2
            DO J=J1,J2
               DO I = SLICE,SLICE
                  IVERT = IVERT + 1
                  VERTS(3*IVERT-2) = XPLT(SLICE)
                  VERTS(3*IVERT-1) = YPLT(J)
                  VERTS(3*IVERT)   = ZPLT(K)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,NK-1
            DO J=1,NJ-1
               IF (ANY(CELL(CELL_INDEX(SLICE:SLICE+1,J,K))%SOLID)) CYCLE
               IF (FCVAR(SLICE,J,K,CC_FGSC,IAXIS) == CC_CUTCFE) THEN
                  ICF = FCVAR(SLICE,J,K,CC_IDCF,IAXIS) ! store cutcell faces and vertices
                  DO IFACECF=1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE
                     NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                     VERTBEG = IVERTCUT + 1
                     VERTBEG = 3*VERTBEG - 2
                     VERTEND = IVERTCUT + NVF
                     VERTEND = 3*VERTEND
                     DO IVCF=1,NVF
                        IVERTCUT = IVERTCUT + 1
                        IVERTCF=CUT_FACE(ICF)%CFELEM(IVCF+1,IFACECF)
                        VERTS(3*IVERTCUT-2:3*IVERTCUT) = REAL(CUT_FACE(ICF)%XYZVERT(1:3,IVERTCF),FB)
                     ENDDO

                     FACEBEG = 3*(IFACECUT+1) - 2
                     FACEEND = FACEBEG + 3*(NVF-2) - 1
                     FACEPTR(1:3*(NVF-2))        =>FACES(FACEBEG:FACEEND)
                     VERTPTR(1:1+VERTEND-VERTBEG)=>VERTS(VERTBEG:VERTEND)
                     VERT_OFFSET = IVERTCUT - NVF
                     ALLOCATE(LOCTYPE(NVF-2))
                     CALL TRIANGULATE(DIR,VERTPTR,NVF,VERT_OFFSET,FACEPTR,LOCTYPE)
                     DO IVCF = 1, NVF-2 ! for now assume face is convex
                        ! vertex indices 1, 2, ..., NVF
                        ! faces (1,2,3), (1,3,4), ..., (1,NVF-1,NVF)
                        IFACECUT = IFACECUT + 1
                        LOCATIONS(IFACECUT) = 2 + LOCTYPE(IVCF)
                        IF(IFACECF > CUT_FACE(ICF)%NFACE) LOCATIONS(IFACECUT) = 1 + LOCTYPE(IVCF) ! Solid side cut-faces.
! after TRIANGULATE is verified remove the following 3 lines of code (and similar lines in 2 locations below)
!                        FACES(3*IFACECUT-2) = (IVERTCUT-NVF)+1
!                        FACES(3*IFACECUT-1) = (IVERTCUT-NVF)+1+IVCF
!                        FACES(3*IFACECUT)   = (IVERTCUT-NVF)+2+IVCF
                     ENDDO
                     DEALLOCATE(LOCTYPE)
                  ENDDO
               ELSE
                  IFACE = IFACE + 1 ! store solid and gas faces and vertices (2 faces per cell)
                  LOCATIONS(IFACE) = 0 + 16
                  IF ( FCVAR(SLICE,J,K,CC_FGSC,IAXIS) == CC_SOLID) LOCATIONS(IFACE)=1 + 16
                  FACES(3*IFACE-2) = IJK(  J,  K,NJ)
                  FACES(3*IFACE-1) = IJK(J+1,  K,NJ)
                  FACES(3*IFACE)   = IJK(J+1,K+1,NJ)

                  IFACE = IFACE + 1
                  LOCATIONS(IFACE) = 0 + 4
                  IF ( FCVAR(SLICE,J,K,CC_FGSC,IAXIS) == CC_SOLID) LOCATIONS(IFACE)=1 + 4
                  FACES(3*IFACE-2) = IJK(  J,  K,NJ)
                  FACES(3*IFACE-1) = IJK(J+1,K+1,NJ)
                  FACES(3*IFACE)   = IJK(  J,K+1,NJ)
               ENDIF
            ENDDO
         ENDDO
      ELSE IF (DIR==2) THEN
         DO K=K1,K2
            DO J=SLICE,SLICE
               DO I = I1,I2
                  IVERT = IVERT + 1
                  VERTS(3*IVERT-2) = XPLT(I)
                  VERTS(3*IVERT-1) = YPLT(SLICE)
                  VERTS(3*IVERT)   = ZPLT(K)
               ENDDO
            ENDDO
         ENDDO
         DO K=1,NK-1
            DO I=1,NI-1
               IF (ANY(CELL(CELL_INDEX(I,SLICE:SLICE+1,K))%SOLID)) CYCLE
               IF (FCVAR(I,SLICE,K,CC_FGSC,JAXIS) == CC_CUTCFE) THEN
                  ICF = FCVAR(I,SLICE,K,CC_IDCF,JAXIS)
                  DO IFACECF=1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE
                     NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                     VERTBEG = IVERTCUT + 1
                     VERTBEG = 3*VERTBEG - 2
                     VERTEND = IVERTCUT + NVF
                     VERTEND = 3*VERTEND
                     DO IVCF=1,NVF
                        IVERTCUT = IVERTCUT + 1
                        IVERTCF=CUT_FACE(ICF)%CFELEM(IVCF+1,IFACECF)
                        VERTS(3*IVERTCUT-2:3*IVERTCUT) = REAL(CUT_FACE(ICF)%XYZVERT(1:3,IVERTCF),FB)
                     ENDDO
                     FACEBEG = 3*(IFACECUT+1) - 2
                     FACEEND = FACEBEG + 3*(NVF-2) - 1
                     FACEPTR(1:3*(NVF-2))        =>FACES(FACEBEG:FACEEND)
                     VERTPTR(1:1+VERTEND-VERTBEG)=>VERTS(VERTBEG:VERTEND)
                     VERT_OFFSET = IVERTCUT - NVF
                     ALLOCATE(LOCTYPE(NVF-2))
                     CALL TRIANGULATE(DIR,VERTPTR,NVF,VERT_OFFSET,FACEPTR,LOCTYPE)
                     DO IVCF = 1, NVF-2 ! for now assume face is convex
                        IFACECUT = IFACECUT + 1
                        LOCATIONS(IFACECUT) = 2 + LOCTYPE(IVCF)
                        IF(IFACECF > CUT_FACE(ICF)%NFACE) LOCATIONS(IFACECUT) = 1 + LOCTYPE(IVCF) ! Solid side cut-faces.
!                        FACES(3*IFACECUT-2) = IVERTCUT-NVF+1
!                        FACES(3*IFACECUT-1) = IVERTCUT-NVF+1+IVCF
!                        FACES(3*IFACECUT)   = IVERTCUT-NVF+1+IVCF+1
                     ENDDO
                     DEALLOCATE(LOCTYPE)
                  ENDDO
               ELSE
                  IFACE = IFACE + 1
                  LOCATIONS(IFACE) = 0 + 16
                  IF ( FCVAR(I,SLICE,K,CC_FGSC,JAXIS) == CC_SOLID) LOCATIONS(IFACE)=1 + 16
                  FACES(3*IFACE-2) = IJK(  I,  K,NI)
                  FACES(3*IFACE-1) = IJK(I+1,  K,NI)
                  FACES(3*IFACE)   = IJK(I+1,K+1,NI)

                  IFACE = IFACE + 1
                  LOCATIONS(IFACE) = 0 + 4
                  IF ( FCVAR(I,SLICE,K,CC_FGSC,JAXIS) == CC_SOLID) LOCATIONS(IFACE)=1 + 4
                  FACES(3*IFACE-2) = IJK(  I,  K,NI)
                  FACES(3*IFACE-1) = IJK(I+1,K+1,NI)
                  FACES(3*IFACE)   = IJK(  I,K+1,NI)
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO K=SLICE,SLICE
            DO J=J1,J2
               DO I = I1,I2
                  IVERT = IVERT + 1
                  VERTS(3*IVERT-2) = XPLT(I)
                  VERTS(3*IVERT-1) = YPLT(J)
                  VERTS(3*IVERT)   = ZPLT(SLICE)
               ENDDO
            ENDDO
         ENDDO
         DO J=1,NJ-1
            DO I=1,NI-1
               IF (ANY(CELL(CELL_INDEX(I,J,SLICE:SLICE+1))%SOLID)) CYCLE
               IF (FCVAR(I,J,SLICE,CC_FGSC,KAXIS) == CC_CUTCFE) THEN
                  ICF = FCVAR(I,J,SLICE,CC_IDCF,KAXIS)
                  DO IFACECF=1,CUT_FACE(ICF)%NFACE+CUT_FACE(ICF)%NSFACE
                     NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                     VERTBEG = IVERTCUT + 1
                     VERTBEG = 3*VERTBEG - 2
                     VERTEND = IVERTCUT + NVF
                     VERTEND = 3*VERTEND
                     DO IVCF=1,NVF
                        IVERTCUT = IVERTCUT + 1
                        IVERTCF=CUT_FACE(ICF)%CFELEM(IVCF+1,IFACECF)
                        VERTS(3*IVERTCUT-2:3*IVERTCUT) = REAL(CUT_FACE(ICF)%XYZVERT(1:3,IVERTCF),FB)
                     ENDDO
                     FACEBEG = 3*(IFACECUT+1) - 2
                     FACEEND = FACEBEG + 3*(NVF-2) - 1
                     FACEPTR(1:3*(NVF-2))        =>FACES(FACEBEG:FACEEND)
                     VERTPTR(1:1+VERTEND-VERTBEG)=>VERTS(VERTBEG:VERTEND)
                     VERT_OFFSET = IVERTCUT - NVF
                     ALLOCATE(LOCTYPE(NVF-2))
                     CALL TRIANGULATE(DIR,VERTPTR,NVF,VERT_OFFSET,FACEPTR,LOCTYPE)
                     DO IVCF = 1, NVF-2 ! for now assume face is convex
                        IFACECUT = IFACECUT + 1
                        LOCATIONS(IFACECUT) = 2 + LOCTYPE(IVCF)
                        IF(IFACECF > CUT_FACE(ICF)%NFACE) LOCATIONS(IFACECUT) = 1 + LOCTYPE(IVCF) ! Solid side cut-faces.
!                        FACES(3*IFACECUT-2) = IVERTCUT-NVF+1
!                        FACES(3*IFACECUT-1) = IVERTCUT-NVF+1+IVCF
!                        FACES(3*IFACECUT)   = IVERTCUT-NVF+1+IVCF+1
                     ENDDO
                     DEALLOCATE(LOCTYPE)
                  ENDDO
               ELSE
                  IFACE = IFACE + 1
                  LOCATIONS(IFACE) = 0 + 16
                  IF ( FCVAR(I,J,SLICE,CC_FGSC,KAXIS) == CC_SOLID) LOCATIONS(IFACE)=1 + 16
                  FACES(3*IFACE-2) = IJK(  I,  J,NI)
                  FACES(3*IFACE-1) = IJK(I+1,  J,NI)
                  FACES(3*IFACE)   = IJK(I+1,J+1,NI)

                  IFACE = IFACE + 1
                  LOCATIONS(IFACE) = 0 + 4
                  IF ( FCVAR(I,J,SLICE,CC_FGSC,KAXIS) == CC_SOLID) LOCATIONS(IFACE)=1 + 4
                  FACES(3*IFACE-2) = IJK(  I,  J,NI)
                  FACES(3*IFACE-1) = IJK(I+1,J+1,NI)
                  FACES(3*IFACE)   = IJK(  I,J+1,NI)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
   ELSE IF (SLICETYPE_LOCAL=='INBOUND_FACES') THEN
      DIR   = 0
      IVERTCUT=NVERTS-NVERTS_CUTCELLS ! start cutcell counters after 'regular' cells
      IFACECUT=NFACES-NFACES_CUTCELLS
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
            IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
            IF (CCVAR(I,J,K,CC_IDCF) > 0) THEN
               ICF = CCVAR(I,J,K,CC_IDCF)
               DO IFACECF=1,CUT_FACE(ICF)%NFACE
                  NVF=CUT_FACE(ICF)%CFELEM(1,IFACECF)
                  VERTBEG = IVERTCUT + 1
                  VERTBEG = 3*VERTBEG - 2
                  VERTEND = IVERTCUT + NVF
                  VERTEND = 3*VERTEND
                  DO IVCF=1,NVF
                     IVERTCUT = IVERTCUT + 1
                     IVERTCF=CUT_FACE(ICF)%CFELEM(IVCF+1,IFACECF)
                     VERTS(3*IVERTCUT-2:3*IVERTCUT) = REAL(CUT_FACE(ICF)%XYZVERT(1:3,IVERTCF),FB)
                  ENDDO
                  IF(PRESENT(SURFIND)) SURFIND(IFACECUT+1:IFACECUT+NVF-2) = CUT_FACE(ICF)%SURF_INDEX(IFACECF)
                  IF(PRESENT(GEOMIND)) GEOMIND(IFACECUT+1:IFACECUT+NVF-2) = CUT_FACE(ICF)%  BODTRI(1,IFACECF)
                  FACEBEG = 3*(IFACECUT+1) - 2
                  FACEEND = FACEBEG + 3*(NVF-2) - 1
                  FACEPTR(1:3*(NVF-2))        =>FACES(FACEBEG:FACEEND)
                  VERTPTR(1:1+VERTEND-VERTBEG)=>VERTS(VERTBEG:VERTEND)
                  VERT_OFFSET = IVERTCUT - NVF
                  ALLOCATE(LOCTYPE(NVF-2))
                  CALL TRIANGULATE(DIR,VERTPTR,NVF,VERT_OFFSET,FACEPTR,LOCTYPE)
                  DO IVCF = 1, NVF-2 ! for now assume face is convex
                     IFACECUT = IFACECUT + 1
                     LOCATIONS(IFACECUT) = 1 + LOCTYPE(IVCF) ! Consider them as SOLID.
                  ENDDO
                  DEALLOCATE(LOCTYPE)
               ENDDO
            ENDIF
            ENDDO
         ENDDO
      ENDDO
   ELSE IF (SLICETYPE_LOCAL=='CUT_CELLS') THEN
      IVERTCUT=NVERTS-NVERTS_CUTCELLS ! start cutcell counters after 'regular' cells
      IFACECUT=NFACES-NFACES_CUTCELLS
      DO KK = 1, KBAR
         DO JJ = 1, JBAR
            DO II = 1, IBAR
               IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID) CYCLE
               IF (CCVAR(II,JJ,KK,CC_IDCC) <= 0) CYCLE
               ICC = CCVAR(II,JJ,KK,CC_IDCC)
               DO JCC=1,CUT_CELL(ICC)%NCELL
                  NFC=CUT_CELL(ICC)%CCELEM(1,JCC)
                  ! Loop on faces corresponding to cut-cell ICC2:
                  DO ICCF=1,NFC
                     IFACE=CUT_CELL(ICC)%CCELEM(ICCF+1,JCC)
                     SELECT CASE(CUT_CELL(ICC)%FACE_LIST(1,IFACE))
                     CASE(CC_FTYPE_RCGAS) ! REGULAR GASPHASE
                        LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
                        X1AXIS  = CUT_CELL(ICC)%FACE_LIST(3,IFACE)
                        ILH     = LOWHIGH - 1
                        I=II; J=JJ; K=KK;
                        SELECT CASE(X1AXIS)
                        CASE(IAXIS)
                           I=II-1+ILH
                           DO KADD=-1,0
                              DO JADD=-1,0
                                 IVERTCUT = IVERTCUT + 1
                                 VERTS(3*IVERTCUT-2) = REAL(X(I     ),FB)
                                 VERTS(3*IVERTCUT-1) = REAL(Y(J+JADD),FB)
                                 VERTS(3*IVERTCUT)   = REAL(Z(K+KADD),FB)
                              ENDDO
                           ENDDO
                        CASE(JAXIS)
                           J=JJ-1+ILH
                           DO IADD=-1,0
                              DO KADD=-1,0
                                 IVERTCUT = IVERTCUT + 1
                                 VERTS(3*IVERTCUT-2) = REAL(X(I+IADD),FB)
                                 VERTS(3*IVERTCUT-1) = REAL(Y(J     ),FB)
                                 VERTS(3*IVERTCUT)   = REAL(Z(K+KADD),FB)
                              ENDDO
                           ENDDO
                        CASE(KAXIS)
                           K=KK-1+ILH
                           DO JADD=-1,0
                              DO IADD=-1,0
                                 IVERTCUT = IVERTCUT + 1
                                 VERTS(3*IVERTCUT-2) = REAL(X(I+IADD),FB)
                                 VERTS(3*IVERTCUT-1) = REAL(Y(J+JADD),FB)
                                 VERTS(3*IVERTCUT)   = REAL(Z(K     ),FB)
                              ENDDO
                           ENDDO
                        END SELECT
                        IFACECUT = IFACECUT + 1
                        LOCATIONS(IFACECUT) = 0 + 16
                        FACES(3*IFACECUT-2:3*IFACECUT) = (/ IVERTCUT-3, IVERTCUT-2, IVERTCUT   /) ! Local Nodes 1, 2, 4

                        IFACECUT = IFACECUT + 1
                        LOCATIONS(IFACECUT) = 0 + 16
                        FACES(3*IFACECUT-2:3*IFACECUT) = (/ IVERTCUT  , IVERTCUT-1, IVERTCUT-3 /) ! Local Nodes 4, 3, 1
                     CASE(CC_FTYPE_CFGAS)
                        ICF2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
                        IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
                        X1AXIS  = CUT_FACE(ICF2)%IJK(KAXIS+1); DIR = X1AXIS
                        NVF     = CUT_FACE(ICF2)%CFELEM(1,IFACE2)
                        VERTBEG = IVERTCUT + 1
                        VERTBEG = 3*VERTBEG - 2
                        VERTEND = IVERTCUT + NVF
                        VERTEND = 3*VERTEND
                        DO IVCF=1,NVF
                           IVERTCUT = IVERTCUT + 1
                           IVERTCF=CUT_FACE(ICF2)%CFELEM(IVCF+1,IFACE2)
                           VERTS(3*IVERTCUT-2:3*IVERTCUT) = REAL(CUT_FACE(ICF2)%XYZVERT(1:3,IVERTCF),FB)
                        ENDDO
                        FACEBEG = 3*(IFACECUT+1) - 2
                        FACEEND = FACEBEG + 3*(NVF-2) - 1
                        FACEPTR(1:3*(NVF-2))        =>FACES(FACEBEG:FACEEND)
                        VERTPTR(1:1+VERTEND-VERTBEG)=>VERTS(VERTBEG:VERTEND)
                        VERT_OFFSET = IVERTCUT - NVF
                        ALLOCATE(LOCTYPE(NVF-2))
                        CALL TRIANGULATE(DIR,VERTPTR,NVF,VERT_OFFSET,FACEPTR,LOCTYPE)
                        DO IVCF = 1, NVF-2 ! for now assume face is convex
                           IFACECUT = IFACECUT + 1
                           LOCATIONS(IFACECUT) = 2 + LOCTYPE(IVCF)
                           IF(IFACE2 > CUT_FACE(ICF2)%NFACE) LOCATIONS(IFACECUT) = 1 + LOCTYPE(IVCF) ! Solid side.
                        ENDDO
                        DEALLOCATE(LOCTYPE)
                     CASE(CC_FTYPE_CFINB)
                        ICF2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
                        IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
                        NVF     = CUT_FACE(ICF2)%CFELEM(1,IFACE2); DIR = 0
                        VERTBEG = IVERTCUT + 1
                        VERTBEG = 3*VERTBEG - 2
                        VERTEND = IVERTCUT + NVF
                        VERTEND = 3*VERTEND
                        DO IVCF=1,NVF
                           IVERTCUT = IVERTCUT + 1
                           IVERTCF=CUT_FACE(ICF2)%CFELEM(IVCF+1,IFACE2)
                           VERTS(3*IVERTCUT-2:3*IVERTCUT) = REAL(CUT_FACE(ICF2)%XYZVERT(1:3,IVERTCF),FB)
                        ENDDO
                        FACEBEG = 3*(IFACECUT+1) - 2
                        FACEEND = FACEBEG + 3*(NVF-2) - 1
                        FACEPTR(1:3*(NVF-2))        =>FACES(FACEBEG:FACEEND)
                        VERTPTR(1:1+VERTEND-VERTBEG)=>VERTS(VERTBEG:VERTEND)
                        VERT_OFFSET = IVERTCUT - NVF
                        ALLOCATE(LOCTYPE(NVF-2))
                        CALL TRIANGULATE(DIR,VERTPTR,NVF,VERT_OFFSET,FACEPTR,LOCTYPE)
                        DO IVCF = 1, NVF-2 ! for now assume face is convex
                           IFACECUT = IFACECUT + 1
                           LOCATIONS(IFACECUT) = 1 + LOCTYPE(IVCF) ! Consider them as SOLID.
                        ENDDO
                        DEALLOCATE(LOCTYPE)
                     END SELECT
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF
END SUBROUTINE GET_GEOMINFO


SUBROUTINE GETSLICEDIR(I1,I2,J1,J2,K1,K2,DIR,SLICE)
INTEGER, INTENT(IN) :: I1, I2, J1, J2, K1, K2
INTEGER, INTENT(OUT) :: DIR, SLICE

IF (ABS(K1-K2)<MIN(ABS(I1-I2),ABS(J1-J2))) THEN
   DIR=3
   SLICE = K1
ELSE IF (ABS(J1-J2)<MIN(ABS(I1-I2),ABS(K1-K2))) THEN
   DIR=2
   SLICE = J1
ELSE
   DIR=1
   SLICE = I1
ENDIF
RETURN

END SUBROUTINE GETSLICEDIR


INTEGER FUNCTION IJK(I,J,NI)
INTEGER, INTENT(IN) :: I, J, NI
IJK = I + (J-1)*NI
END FUNCTION IJK


SUBROUTINE GET_GASCUTFACE_SCALAR_SLICE(VAL_CF,X1AXIS,ICF,IFACE,CC_INTERP2FACES,CC_CELL_CENTERED,&
                         IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)

USE PHYSICAL_FUNCTIONS, ONLY: GET_MASS_FRACTION

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: X1AXIS,ICF,IFACE,&
                       IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,NM
LOGICAL, INTENT(IN) :: CC_INTERP2FACES,CC_CELL_CENTERED
REAL(EB),INTENT(OUT):: VAL_CF

! Local Variables:
REAL(EB) :: X1F, IDX, CCM1, CCP1, VAL_LOC(LOW_IND:HIGH_IND)
INTEGER  :: ISIDE, ICC, JCC, LOCAL_IND, II, JJ, KK
REAL(EB) :: Y_SPECIES(LOW_IND:HIGH_IND)
! REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES)

! Point to mesh has been called for MESHES(NM):

Y_SPECIES(LOW_IND:HIGH_IND) = 1._EB

! Here interpolate values from cut-cell centers:
X1F= CUT_FACE(ICF)%XYZCEN(X1AXIS,IFACE)
IDX= 1._EB/ ( CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE) - &
              CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE) )
CCM1= IDX*(CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE)-X1F)
CCP1= IDX*(X1F-CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE))
LOCAL_IND=HIGH_IND

IF (.NOT.CC_INTERP2FACES .AND. CC_CELL_CENTERED) THEN
   CCM1=1._EB
   CCP1=0._EB
   LOCAL_IND=LOW_IND
ENDIF

VAL_LOC(LOW_IND:HIGH_IND)= 0._EB
DO ISIDE=LOW_IND,LOCAL_IND
   SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE,IFACE))
   CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use value from CUT_CELL data struct:
      ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE,IFACE)
      JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE,IFACE)
      II = CUT_CELL(ICC)%IJK(IAXIS)
      JJ = CUT_CELL(ICC)%IJK(JAXIS)
      KK = CUT_CELL(ICC)%IJK(KAXIS)
      VAL_LOC(ISIDE) = GAS_PHASE_OUTPUT(T,DT,NM,II,JJ,KK,&
                       IND,IND2,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,&
                       PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,ICC,JCC)
   END SELECT
ENDDO
VAL_CF = CCM1*VAL_LOC(LOW_IND) + CCP1*VAL_LOC(HIGH_IND)

RETURN
END SUBROUTINE GET_GASCUTFACE_SCALAR_SLICE


SUBROUTINE GET_SOLIDCUTFACE_SCALAR_SLICE(X1AXIS,ICF,VAL_CF, &
              IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)

USE PHYSICAL_FUNCTIONS, ONLY: GET_MASS_FRACTION

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: X1AXIS,ICF,IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,NM
REAL(EB),INTENT(OUT):: VAL_CF

! Local Variables:
INTEGER :: II_LO,II_HI,JJ_LO,JJ_HI,KK_LO,KK_HI,IJK(IAXIS:KAXIS),IJK2(IAXIS:KAXIS,16),ICELL,II,JJ,KK
LOGICAL :: FOUND
REAL(EB):: Y_SPECIES

! Point to mesh has been called for MESHES(NM): This routine searches for a REGULAR SOLID cell in the
! vicinity of the SOLID cut-face and assigns to the latter the scalar value of the former.

VAL_CF    = 0._EB
Y_SPECIES = 1._EB

IJK(IAXIS:KAXIS)=CUT_FACE(ICF)%IJK(IAXIS:KAXIS)

SELECT CASE(X1AXIS)
CASE(IAXIS)
   II_LO=IJK(IAXIS);   II_HI=IJK(IAXIS)+1
   JJ_LO=IJK(JAXIS)-1; JJ_HI=IJK(JAXIS)+1
   KK_LO=IJK(KAXIS)-1; KK_HI=IJK(KAXIS)+1

   IJK2(IAXIS:KAXIS, 1) = (/ II_LO, JJ_LO, IJK(KAXIS) /)
   IJK2(IAXIS:KAXIS, 2) = (/ II_LO, JJ_HI, IJK(KAXIS) /)
   IJK2(IAXIS:KAXIS, 3) = (/ II_LO, IJK(JAXIS), KK_LO /)
   IJK2(IAXIS:KAXIS, 4) = (/ II_LO, IJK(JAXIS), KK_HI /)
   IJK2(IAXIS:KAXIS, 5) = (/ II_HI, JJ_LO, IJK(KAXIS) /)
   IJK2(IAXIS:KAXIS, 6) = (/ II_HI, JJ_HI, IJK(KAXIS) /)
   IJK2(IAXIS:KAXIS, 7) = (/ II_HI, IJK(JAXIS), KK_LO /)
   IJK2(IAXIS:KAXIS, 8) = (/ II_HI, IJK(JAXIS), KK_HI /)
   IJK2(IAXIS:KAXIS, 9) = (/ II_LO, JJ_LO, KK_LO /)
   IJK2(IAXIS:KAXIS,10) = (/ II_LO, JJ_HI, KK_LO /)
   IJK2(IAXIS:KAXIS,11) = (/ II_LO, JJ_LO, KK_HI /)
   IJK2(IAXIS:KAXIS,12) = (/ II_LO, JJ_HI, KK_HI /)
   IJK2(IAXIS:KAXIS,13) = (/ II_HI, JJ_LO, KK_LO /)
   IJK2(IAXIS:KAXIS,14) = (/ II_HI, JJ_HI, KK_LO /)
   IJK2(IAXIS:KAXIS,15) = (/ II_HI, JJ_LO, KK_HI /)
   IJK2(IAXIS:KAXIS,16) = (/ II_HI, JJ_HI, KK_HI /)

CASE(JAXIS)
   II_LO=IJK(IAXIS)-1; II_HI=IJK(IAXIS)+1
   JJ_LO=IJK(JAXIS);   JJ_HI=IJK(JAXIS)+1
   KK_LO=IJK(KAXIS)-1; KK_HI=IJK(KAXIS)+1

   IJK2(IAXIS:KAXIS, 1) = (/ IJK(IAXIS), JJ_LO, KK_LO /)
   IJK2(IAXIS:KAXIS, 2) = (/ IJK(IAXIS), JJ_LO, KK_HI /)
   IJK2(IAXIS:KAXIS, 3) = (/ II_LO, JJ_LO, IJK(KAXIS) /)
   IJK2(IAXIS:KAXIS, 4) = (/ II_HI, JJ_LO, IJK(KAXIS) /)
   IJK2(IAXIS:KAXIS, 5) = (/ IJK(IAXIS), JJ_HI, KK_LO /)
   IJK2(IAXIS:KAXIS, 6) = (/ IJK(IAXIS), JJ_HI, KK_HI /)
   IJK2(IAXIS:KAXIS, 7) = (/ II_LO, JJ_HI, IJK(KAXIS) /)
   IJK2(IAXIS:KAXIS, 8) = (/ II_HI, JJ_HI, IJK(KAXIS) /)
   IJK2(IAXIS:KAXIS, 9) = (/ II_LO, JJ_LO, KK_LO /)
   IJK2(IAXIS:KAXIS,10) = (/ II_LO, JJ_LO, KK_HI /)
   IJK2(IAXIS:KAXIS,11) = (/ II_HI, JJ_LO, KK_LO /)
   IJK2(IAXIS:KAXIS,12) = (/ II_HI, JJ_LO, KK_HI /)
   IJK2(IAXIS:KAXIS,13) = (/ II_LO, JJ_HI, KK_LO /)
   IJK2(IAXIS:KAXIS,14) = (/ II_LO, JJ_HI, KK_HI /)
   IJK2(IAXIS:KAXIS,15) = (/ II_HI, JJ_HI, KK_LO /)
   IJK2(IAXIS:KAXIS,16) = (/ II_HI, JJ_HI, KK_HI /)

CASE(KAXIS)
   II_LO=IJK(IAXIS)-1; II_HI=IJK(IAXIS)+1
   JJ_LO=IJK(JAXIS)-1; JJ_HI=IJK(JAXIS)+1
   KK_LO=IJK(KAXIS);   KK_HI=IJK(KAXIS)+1

   IJK2(IAXIS:KAXIS, 1) = (/ II_LO, IJK(JAXIS), KK_LO /)
   IJK2(IAXIS:KAXIS, 2) = (/ II_HI, IJK(JAXIS), KK_LO /)
   IJK2(IAXIS:KAXIS, 3) = (/ IJK(IAXIS), JJ_LO, KK_LO /)
   IJK2(IAXIS:KAXIS, 4) = (/ IJK(IAXIS), JJ_HI, KK_LO /)
   IJK2(IAXIS:KAXIS, 5) = (/ II_LO, IJK(JAXIS), KK_HI /)
   IJK2(IAXIS:KAXIS, 6) = (/ II_HI, IJK(JAXIS), KK_HI /)
   IJK2(IAXIS:KAXIS, 7) = (/ IJK(IAXIS), JJ_LO, KK_HI /)
   IJK2(IAXIS:KAXIS, 8) = (/ IJK(IAXIS), JJ_HI, KK_HI /)
   IJK2(IAXIS:KAXIS, 9) = (/ II_LO, JJ_LO, KK_LO /)
   IJK2(IAXIS:KAXIS,10) = (/ II_HI, JJ_LO, KK_LO /)
   IJK2(IAXIS:KAXIS,11) = (/ II_LO, JJ_HI, KK_LO /)
   IJK2(IAXIS:KAXIS,12) = (/ II_HI, JJ_HI, KK_LO /)
   IJK2(IAXIS:KAXIS,13) = (/ II_LO, JJ_LO, KK_HI /)
   IJK2(IAXIS:KAXIS,14) = (/ II_HI, JJ_LO, KK_HI /)
   IJK2(IAXIS:KAXIS,15) = (/ II_LO, JJ_HI, KK_HI /)
   IJK2(IAXIS:KAXIS,16) = (/ II_HI, JJ_HI, KK_HI /)

END SELECT

FOUND=.FALSE.
DO ICELL=1,16
   ! Look only for internal cells:
   II=IJK2(IAXIS,ICELL)
   IF(II < 1 .OR. II > IBAR) CYCLE
   JJ=IJK2(JAXIS,ICELL)
   IF(JJ < 1 .OR. JJ > JBAR) CYCLE
   KK=IJK2(KAXIS,ICELL)
   IF(KK < 1 .OR. KK > KBAR) CYCLE
   IF (CCVAR(II,JJ,KK,CC_CGSC) /= CC_SOLID) CYCLE
   FOUND=.TRUE.
   EXIT
ENDDO

IF(.NOT.FOUND) THEN ! This is a thin object. Use first gas cut-cell value:
   DO ICELL=1,16
      ! Look only for internal cells:
      II=IJK2(IAXIS,ICELL)
      IF(II < 1 .OR. II > IBAR) CYCLE
      JJ=IJK2(JAXIS,ICELL)
      IF(JJ < 1 .OR. JJ > JBAR) CYCLE
      KK=IJK2(KAXIS,ICELL)
      IF(KK < 1 .OR. KK > KBAR) CYCLE
      IF (CCVAR(II,JJ,KK,CC_CGSC) /= CC_CUTCFE) CYCLE
      FOUND=.TRUE.
      EXIT
   ENDDO
ENDIF

! Use closest solid Cell values for SOLID cut-face:
IF (FOUND) THEN
   VAL_CF = GAS_PHASE_OUTPUT(T,DT,NM,II,JJ,KK,&
                             IND,IND2,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX)
ENDIF

RETURN
END SUBROUTINE GET_SOLIDCUTFACE_SCALAR_SLICE


SUBROUTINE GET_SOLIDREGFACE_SCALAR_SLICE(X1AXIS,I,J,K,VAL_CF,&
              IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,T,DT,NM)

USE PHYSICAL_FUNCTIONS, ONLY: GET_MASS_FRACTION

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: X1AXIS,I,J,K,IND,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX,NM
REAL(EB),INTENT(OUT):: VAL_CF

! Local Variables:
INTEGER :: II_LO,II_HI,JJ_LO,JJ_HI,KK_LO,KK_HI,SOLID_LO,SOLID_HI
REAL(EB):: CC1(LOW_IND:HIGH_IND),CCSUM
REAL(EB) :: Y_SPECIES,VAL_CF_LO,VAL_CF_HI

VAL_CF    = 0._EB
Y_SPECIES = 1._EB
VAL_CF_LO = 0._EB
VAL_CF_HI = 0._EB

SELECT CASE(X1AXIS)
CASE(IAXIS)
   II_LO=I; II_HI=I+1
   JJ_LO=J; JJ_HI=J
   KK_LO=K; KK_HI=K
CASE(JAXIS)
   II_LO=I; II_HI=I
   JJ_LO=J; JJ_HI=J+1
   KK_LO=K; KK_HI=K
CASE(KAXIS)
   II_LO=I; II_HI=I
   JJ_LO=J; JJ_HI=J
   KK_LO=K; KK_HI=K+1
END SELECT

SOLID_LO = CCVAR(II_LO,JJ_LO,KK_LO,CC_CGSC)
SOLID_HI = CCVAR(II_HI,JJ_HI,KK_HI,CC_CGSC)

! This discards interpolation from Adjacent cut-cells:
CC1(LOW_IND:HIGH_IND) = 0._EB
IF(SOLID_LO == CC_SOLID) CC1( LOW_IND)= 1._EB
IF(SOLID_HI == CC_SOLID) CC1(HIGH_IND)= 1._EB

! Interpolation coefficients:
CCSUM = SUM(CC1(LOW_IND:HIGH_IND))
IF( CCSUM > 0._EB ) CC1(LOW_IND:HIGH_IND)=CC1(LOW_IND:HIGH_IND)/CCSUM

IF (CC1( LOW_IND)>TWENTY_EPSILON_EB) THEN
   VAL_CF_LO = GAS_PHASE_OUTPUT(T,DT,NM,II_LO,JJ_LO,KK_LO,&
                                IND,IND2,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX)
ENDIF

IF (CC1(HIGH_IND)>TWENTY_EPSILON_EB) THEN
   VAL_CF_HI = GAS_PHASE_OUTPUT(T,DT,NM,II_HI,JJ_HI,KK_HI,&
                                IND,IND2,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX,MATL_INDEX)
ENDIF

VAL_CF = CC1(LOW_IND)*VAL_CF_LO + CC1(HIGH_IND)*VAL_CF_HI

RETURN
END SUBROUTINE GET_SOLIDREGFACE_SCALAR_SLICE


!> \brief Get values for boundary output
!>
!> \param IP Patch number
!> \param NM Mesh number

SUBROUTINE GET_BNDF_PACK(T,NM,NF,IP,PP,PPN)
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NF,NM,IP
REAL(FB), POINTER, DIMENSION(:,:), INTENT(IN) :: PP,PPN
TYPE(PATCH_TYPE), POINTER :: PA
INTEGER :: IND,ISUM,I,J,K,IC,IW,L,L1,L2,N,N1,N2,NC

PA => PATCH(IP)

PP  = REAL(OUTPUT_QUANTITY(-IND)%AMBIENT_VALUE,FB)
PPN = 0._FB
IBK = 0

! Adjust PATCH indices depending on orientation

SELECT CASE(ABS(PA%IOR))
   CASE(1) ; L1=PA%JG1 ; L2=PA%JG2 ; N1=PA%KG1 ; N2=PA%KG2
   CASE(2) ; L1=PA%IG1 ; L2=PA%IG2 ; N1=PA%KG1 ; N2=PA%KG2
   CASE(3) ; L1=PA%IG1 ; L2=PA%IG2 ; N1=PA%JG1 ; N2=PA%JG2
END SELECT

! Evaluate the given boundary quantity at each cell of the current PATCH

DO K=PA%KG1,PA%KG2
   DO J=PA%JG1,PA%JG2
      DO I=PA%IG1,PA%IG2
         IC = CELL_INDEX(I,J,K)
         IW = CELL(IC)%WALL_INDEX(-PA%IOR) ; IF (IW==0) CYCLE
         SELECT CASE(ABS(PA%IOR))
            CASE(1) ; L=J ; N=K
            CASE(2) ; L=I ; N=K
            CASE(3) ; L=I ; N=J
         END SELECT
         IF (WALL(IW)%BOUNDARY_TYPE/=NULL_BOUNDARY .AND. &
            WALL(IW)%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. .NOT.CELL(IC)%SOLID) THEN
            IBK(L,N) = 1
            PP(L,N)  = REAL(SOLID_PHASE_OUTPUT(IND,T,NM,BF%Y_INDEX,BF%Z_INDEX,BF%PART_INDEX,OPT_WALL_INDEX=IW,&
                                               OPT_BNDF_INDEX=NF),FB)
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Integrate the boundary quantity in time

IF (BNDF_COUNTER(NM)>0 .AND. BF%TIME_INTEGRAL_INDEX>0) THEN
   DO N=N1,N2
      DO L=L1,L2
         NC = NC + 1
         BNDF_TIME_INTEGRAL(NC,BF%TIME_INTEGRAL_INDEX) = BNDF_TIME_INTEGRAL(NC,BF%TIME_INTEGRAL_INDEX) + &
                            PP(L,N)*REAL(BNDF_CLOCK(BNDF_COUNTER(NM))-BNDF_CLOCK(BNDF_COUNTER(NM)-1),FB)
         PP(L,N) = BNDF_TIME_INTEGRAL(NC,BF%TIME_INTEGRAL_INDEX)
      ENDDO
   ENDDO
ENDIF

! Interpolate the boundary quantity PP at cell corners, PPN

IF (.NOT.BF%CELL_CENTERED) THEN

   ! Dont include undetermined values in interpolation for FIRE ARRIVAL TIME
   IF (OUTPUT_QUANTITY(BF%INDEX)%NAME=='FIRE ARRIVAL TIME') THEN
      WHERE(PP>9.E5_FB) IBK=0
    ENDIF

   DO N=N1-1,N2
      DO L=L1-1,L2
         IF (IBK(L,N)==1)     PPN(L,N) = PPN(L,N) + PP(L,N)
         IF (IBK(L+1,N)==1)   PPN(L,N) = PPN(L,N) + PP(L+1,N)
         IF (IBK(L,N+1)==1)   PPN(L,N) = PPN(L,N) + PP(L,N+1)
         IF (IBK(L+1,N+1)==1) PPN(L,N) = PPN(L,N) + PP(L+1,N+1)
         ISUM = IBK(L,N)+IBK(L,N+1)+IBK(L+1,N)+IBK(L+1,N+1)
         IF (ISUM>0) THEN
            PPN(L,N) = PPN(L,N)/REAL(ISUM,FB)
         ELSE
            PPN(L,N) = REAL(SOLID_PHASE_OUTPUT(IND,T,NM,BF%Y_INDEX,BF%Z_INDEX,BF%PART_INDEX,OPT_WALL_INDEX=0,&
                                               OPT_BNDF_INDEX=NF),FB)
         ENDIF
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE GET_BNDF_PACK


!> \brief Get values for smoke3d output
!>
!> \param S3 Current smoke3d object
!> \param T Current simulation time (s)
!> \param DT Current time step size (s)
!> \param NM Mesh number
!> \param FF cell centered values
!> \param QQ node interpolated values

SUBROUTINE GET_SMOKE3D_QQ(S3,T,DT,NM,FF,QQ)
USE ISOSMOKE, ONLY: SMOKE3D_TO_FILE
TYPE(SMOKE3D_TYPE), POINTER, INTENT(IN) :: S3
REAL(FB), POINTER, DIMENSION(:,:,:,:), INTENT(INOUT) :: QQ
REAL(EB), POINTER, DIMENSION(:,:,:), INTENT(INOUT) :: FF
REAL(EB), INTENT(IN) :: T,DT
INTEGER,  INTENT(IN) :: NM
INTEGER  :: I,J,K
REAL(EB) :: FR_C

DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         FF(I,J,K)=GAS_PHASE_OUTPUT(T,DT,NM,I,J,K,S3%QUANTITY_INDEX,0,S3%Y_INDEX,0,S3%Z_INDEX,0,0,0,0,0,0)
      ENDDO
   ENDDO
ENDDO

! Adjust the temperature as it is used in the expression for the radiation source term

IF (S3%DISPLAY_TYPE=='TEMPERATURE' .AND. RTE_SOURCE_CORRECTION) THEN
   FR_C = RTE_SOURCE_CORRECTION_FACTOR**0.25_EB
   WHERE (CHI_R*Q>QR_CLIP) FF = (FF+TMPM)*FR_C - TMPM
ENDIF

! Interpolate data to cell nodes
DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         QQ(I,J,K,1) = REAL((FF(I,J,K)  +FF(I+1,J,K)  +FF(I,J,K+1)  +FF(I+1,J,K+1)+ &
                             FF(I,J+1,K)+FF(I+1,J+1,K)+FF(I,J+1,K+1)+FF(I+1,J+1,K+1))*0.125_FB,FB)
      ENDDO
   ENDDO
ENDDO

IF (CC_IBM) THEN
   DO K=0,KBAR
      DO J=0,JBAR
         DO I=0,IBAR
            IF(MESHES(NM)%VERTVAR(I,J,K,CC_VGSC) /= CC_SOLID) CYCLE
            QQ(I,J,K,1) = 0._FB
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE GET_SMOKE3D_QQ


! \brief Compute the integrals needed for layer height, average upper and lower layer temperatures

SUBROUTINE GET_LAYER_HEIGHT_INTEGRALS(II,JJ,K_LO,K_HI,Z_INT,Z_LO,I_1,I_2,I_3,I_4,TMP_LOW)

INTEGER, INTENT(IN) :: II,JJ,K_LO,K_HI
REAL(EB), INTENT(OUT) :: I_1,I_2,I_3,I_4,TMP_LOW
REAL(EB), INTENT(IN)  :: Z_LO,Z_INT
INTEGER :: K

I_1 = 0._EB
I_2 = 0._EB
I_3 = 0._EB
I_4 = 0._EB
DO K=K_LO,K_HI
   IF (CELL(CELL_INDEX(II,JJ,K))%SOLID) CYCLE
   I_1 = I_1 + DZ(K)*TMP(II,JJ,K)
   I_2 = I_2 + DZ(K)/TMP(II,JJ,K)
   I_4 = I_4 + DZ(K)
   IF (Z(K-1)-Z_LO>=Z_INT) THEN
      I_3 = I_3 + TMP(II,JJ,K)*DZ(K)
   ELSEIF (Z(K)-Z_LO>Z_INT) THEN
      I_3 = I_3 + TMP(II,JJ,K)  *(Z(K)-Z_LO-Z_INT)
   ELSE
   ENDIF
ENDDO
TMP_LOW = TMP(II,JJ,K_LO)

END SUBROUTINE GET_LAYER_HEIGHT_INTEGRALS


REAL(EB) FUNCTION WAVELET_ERROR_MEASURE(II,JJ,KK,IND,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,DT,NM)
REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: II,JJ,KK,IND,NM,VELO_INDEX,Y_INDEX,Z_INDEX,PART_INDEX
REAL(EB) :: SS(4)

! wavelet error measure
WAVELET_ERROR_MEASURE = 0._EB

SS(1) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,MAX(0,II-2),JJ,KK,              IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
SS(2) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,MAX(0,II-1),JJ,KK,              IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
SS(3) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,JJ,KK,                       IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
SS(4) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,MIN(MESHES(NM)%IBP1,II+1),JJ,KK,IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
WAVELET_ERROR_MEASURE = WAVELET_ERROR(SS)

IF (.NOT.TWO_D) THEN
   SS(1) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,MAX(0,JJ-2),KK,              IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
   SS(2) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,MAX(0,JJ-1),KK,              IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
   SS(3) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,JJ,KK,                       IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
   SS(4) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,MIN(MESHES(NM)%JBP1,JJ+1),KK,IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
   WAVELET_ERROR_MEASURE = MAX(WAVELET_ERROR_MEASURE,WAVELET_ERROR(SS))
ENDIF

SS(1) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,JJ,MAX(0,KK-2),              IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
SS(2) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,JJ,MAX(0,KK-1),              IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
SS(3) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,JJ,KK,                       IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
SS(4) = GAS_PHASE_OUTPUT(T_BEGIN,DT,NM,II,JJ,MIN(MESHES(NM)%KBP1,KK+1),IND,0,Y_INDEX,Z_INDEX,0,PART_INDEX,VELO_INDEX,0,0,0,0)
WAVELET_ERROR_MEASURE = MAX(WAVELET_ERROR_MEASURE,WAVELET_ERROR(SS))

END FUNCTION WAVELET_ERROR_MEASURE


REAL(EB) FUNCTION WAVELET_ERROR(S)

INTEGER, PARAMETER :: M=2 ! only need two level transform, but could be generalized
REAL(EB), INTENT(IN) :: S(2*M)
REAL(EB) :: SS(2*M),A(M,M)=0._EB,C(M,M)=0._EB,C1,C2,SMIN,SMAX,DS
INTEGER :: I,J,K,N

! Comments: This function generates a normalized error measure WAVELET_ERROR based on coefficients
! from a simple Haar wavelet transform.  The function requires the input of 4 scalar values.  The
! error is estimated at the point of the value S(3) based on a piece-wise constant reconstruction
! of the underlying function.  For example...
!
!     |<---------- interval --------->|
!
!            S(2)
!             o-------       S(4)
!    S(1)                     o-------
!     o-------
!                    S(3)
!                     o-------
!                     ^
!                     |
!             error computed here

! normalize signal
SMAX=MAXVAL(S)
SMIN=MINVAL(S)
DS=SMAX-SMIN
IF (DS<1.E-6) THEN
   WAVELET_ERROR = 0._EB
   RETURN
ELSE
   SS=(S-SMIN)/DS
ENDIF

! discrete Haar wavelet transform
N=M
DO I=1,M
   DO J=1,N
      K=2*J-1
      IF (I==1) THEN
         A(I,J) = 0.5_EB*(SS(K)+SS(K+1))
         C(I,J) = 0.5_EB*(SS(K)-SS(K+1))
      ELSE
         A(I,J) = 0.5_EB*(A(I-1,K)+A(I-1,K+1))
         C(I,J) = 0.5_EB*(A(I-1,K)-A(I-1,K+1))
      ENDIF
   ENDDO
   N=N/2;
ENDDO

C1 = SUM(C(1,:))
C2 = SUM(C(2,:))

WAVELET_ERROR = ABS(C1-C2)

END FUNCTION WAVELET_ERROR


!> \brief Back out k_sgs (subgrid kinetic energy per unit mass) from Deardorff eddy viscosity

REAL(EB) FUNCTION SUBGRID_KINETIC_ENERGY(MU_TURB,RHO,C_NU,DELTA)

REAL(EB), INTENT(IN) :: MU_TURB,RHO,C_NU,DELTA
REAL(EB) :: DENOM

DENOM = RHO*C_NU*DELTA
IF (DENOM>TWENTY_EPSILON_EB) THEN
   SUBGRID_KINETIC_ENERGY = (MAX(MU_TURB,0._EB)/DENOM)**2
ELSE
   SUBGRID_KINETIC_ENERGY = 0._EB
ENDIF

END FUNCTION SUBGRID_KINETIC_ENERGY


!> \brief Compute gas phase output quantities
!>
!> \param T Current simulation time (s)
!> \param DT Current time step size (s)
!> \param NM Current mesh
!> \param II Cell index in \f$ x \f$ direction
!> \param JJ Cell index in \f$ y \f$ direction
!> \param KK Cell index in \f$ z \f$ direction
!> \param IND Index of the output quantity
!> \param IND2 Index of the sometimes needed second output quantity
!> \param Y_INDEX Index of the primitive gas species
!> \param Z_INDEX Index of the gas species mixture
!> \param ELEM_INDX Index of the chemical element
!> \param PART_INDEX Index of the Lagrangian particle class
!> \param VELO_INDEX Index of the velocity component, x=1, y=2, z=3
!> \param PIPE_INDEX Index of the pipe branch
!> \param PROP_INDEX Index of the PROPerty group parameters
!> \param REAC_INDEX Index of the REACtion
!> \param MATL_INDEX Index of the Material
!> \param ICC_IN,JCC_IN Optional indexes of cut-cell.

REAL(EB) RECURSIVE FUNCTION GAS_PHASE_OUTPUT(T,DT,NM,II,JJ,KK,IND,IND2,Y_INDEX,Z_INDEX,ELEM_INDX,PART_INDEX,VELO_INDEX,PIPE_INDEX,&
                                           PROP_INDEX,REAC_INDEX,MATL_INDEX,ICC_IN,JCC_IN) RESULT(GAS_PHASE_OUTPUT_RES)

USE MEMORY_FUNCTIONS, ONLY: REALLOCATE
USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D,INTERPOLATE1D_UNIFORM,UPDATE_HISTOGRAM
USE PHYSICAL_FUNCTIONS, ONLY: GET_MASS_FRACTION,FED,FIC,GET_SPECIFIC_HEAT,RELATIVE_HUMIDITY, &
                              GET_CONDUCTIVITY,GET_MOLECULAR_WEIGHT,GET_MASS_FRACTION_ALL,GET_ENTHALPY,GET_SENSIBLE_ENTHALPY, &
                              GET_VISCOSITY,GET_POTENTIAL_TEMPERATURE,GET_SPECIFIC_GAS_CONSTANT,&
                              SURFACE_DENSITY, CALC_EQUIV_RATIO,FORCED_CONVECTION_MODEL
USE COMP_FUNCTIONS, ONLY : CURRENT_TIME,SYSTEM_MEM_USAGE
USE TURBULENCE, ONLY: K_SGS_POPE
USE RADCONS, ONLY: WL_LOW, WL_HIGH, RADTMP
USE RAD, ONLY: BLACKBODY_FRACTION
USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_P_3,VD2D_MMS_H_3
USE CC_SCALARS, ONLY: CC_CUTCELL_VELOCITY

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: II,JJ,KK,IND,IND2,NM,Y_INDEX,Z_INDEX,ELEM_INDX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,REAC_INDEX, &
                       MATL_INDEX
INTEGER, INTENT(IN), OPTIONAL :: ICC_IN,JCC_IN
REAL(EB) :: H_TC,TMP_TC,RE_D,NUSSELT,VEL,K_G,MU_G,COSTHETA,FAC,&
            Q_SUM,TMP_G,UU,VV,WW,VEL2,Y_MF_INT,PATHLENGTH,EXT_COEF,MASS_EXT_COEF,ZZ_FUEL,ZZ_OX,&
            VELSR,WATER_VOL_FRAC,RHS,DT_C,DT_E,T_RATIO,Y_E_LAG, H_G,H_G_SUM,CPBAR,CP,ZZ_GET(1:N_TRACKED_SPECIES),RCON,&
            EXPON,Y_SPECIES,MEC,Y_SPECIES2,Y_H2O,R_Y_H2O,R_DN,SGN,Y_ALL(N_SPECIES),H_S,D_Z_N(0:I_MAX_TEMP),&
            DISSIPATION_RATE,S11,S22,S33,S12,S13,S23,DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,ONTHDIV,SS,ETA,DELTA,R_DX2,&
            UVW,UODX,VODY,WODZ,XHAT,ZHAT,BBF,GAMMA_LOC,VC,VOL,PHI,GAS_PHASE_OUTPUT_CC,&
            GAS_PHASE_OUTPUT_CFA,CFACE_AREA,VELOCITY_COMPONENT(1:3),ATOTV(1:3),TMP_F,R_D,MW,PROBE_TMP,PROBE_RHO,PROBE_DELTA_P
INTEGER :: N,I,J,K,NN,IL,III,JJJ,KKK,IP,JP,KP,FED_ACTIVITY,IP1,JP1,KP1,IM1,JM1,KM1,IIM1,JJM1,KKM1,NR,NS,RAM,&
           ICC,JCC,NCELL,AXIS,ICF,NFACE,JCF,JCC_LO,JCC_HI,PDPA_FORMULA,IC
REAL(FB) :: RN
REAL(EB), PARAMETER :: EPS=1.E-10_EB
REAL :: CPUTIME
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D

! Get species mass fraction if necessary

Y_H2O     = 0._EB
R_Y_H2O   = 0._EB
Y_SPECIES = 1._EB

IF (Z_INDEX > 0) THEN
   Y_SPECIES = ZZ(II,JJ,KK,Z_INDEX)
   RCON = SPECIES_MIXTURE(Z_INDEX)%RCON
ELSEIF (Y_INDEX > 0) THEN
   ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
   RCON = SPECIES(Y_INDEX)%RCON
   CALL GET_MASS_FRACTION(ZZ_GET,Y_INDEX,Y_SPECIES)
ENDIF
IF (DRY .AND. H2O_INDEX > 0) THEN
   ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
   CALL GET_MASS_FRACTION(ZZ_GET,H2O_INDEX,Y_H2O)
   R_Y_H2O = SPECIES(H2O_INDEX)%RCON * Y_H2O
   IF (Y_INDEX==H2O_INDEX) Y_SPECIES=0._EB
ENDIF

! Get desired output value

IND_SELECT: SELECT CASE(IND)
   CASE DEFAULT  ! SMOKE/WATER
      GAS_PHASE_OUTPUT_RES = 0._EB
   CASE( 1)  ! DENSITY
      GAS_PHASE_OUTPUT_RES = RHO(II,JJ,KK)*Y_SPECIES
   CASE( 2)  ! F_X
      GAS_PHASE_OUTPUT_RES = FVX(II,JJ,KK)
   CASE( 3)  ! F_Y
      GAS_PHASE_OUTPUT_RES = FVY(II,JJ,KK)
   CASE( 4)  ! F_Z
      GAS_PHASE_OUTPUT_RES = FVZ(II,JJ,KK)
   CASE( 5)  ! TEMPERATURE
      GAS_PHASE_OUTPUT_RES = TMP(II,JJ,KK) - TMPM
   CASE( 6)  ! U-VELOCITY
      GAS_PHASE_OUTPUT_RES = U(II,JJ,KK)
   CASE( 7)  ! V-VELOCITY
      GAS_PHASE_OUTPUT_RES = V(II,JJ,KK)
   CASE( 8)  ! W-VELOCITY
      GAS_PHASE_OUTPUT_RES = W(II,JJ,KK)
   CASE( 9)  ! PRESSURE
      GAS_PHASE_OUTPUT_RES = PBAR(KK,PRESSURE_ZONE(II,JJ,KK)) + &
                             RHO(II,JJ,KK)*(0.5_EB*(H(II,JJ,KK)+HS(II,JJ,KK))-KRES(II,JJ,KK)) - P_0(KK)
   CASE(10)  ! VELOCITY
      SELECT CASE(ABS(VELO_INDEX))
         CASE DEFAULT
            SGN = 1._EB
         CASE(1)
            SGN = SIGN(1._EB,U(II,JJ,KK))*SIGN(1,VELO_INDEX)
         CASE(2)
            SGN = SIGN(1._EB,V(II,JJ,KK))*SIGN(1,VELO_INDEX)
         CASE(3)
            SGN = SIGN(1._EB,W(II,JJ,KK))*SIGN(1,VELO_INDEX)
      END SELECT
      GAS_PHASE_OUTPUT_RES = SGN*SQRT(0.25_EB*((U(MAX(0,II-1),JJ,KK)+U(MIN(IBAR,II),JJ,KK))**2+&
                                               (V(II,MAX(0,JJ-1),KK)+V(II,MIN(JBAR,JJ),KK))**2+&
                                               (W(II,JJ,MAX(0,KK-1))+W(II,JJ,MIN(KBAR,KK)))**2))
   CASE(11)  ! HRRPUV
      GAS_PHASE_OUTPUT_RES = Q(II,JJ,KK)*0.001_EB
   CASE(12)  ! H
      GAS_PHASE_OUTPUT_RES = 0.5_EB*(HS(II,JJ,KK)+H(II,JJ,KK))
   CASE(13)  ! MIXTURE FRACTION
      ! requires FUEL + AIR --> PROD (SIMPLE_CHEMISTRY, N_SIMPLE_CHEMISTRY_REACTIONS=1)
      ! f = Z_FUEL + Z_PROD/(1+S), where S is the mass stoichiometric coefficient for AIR
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NR=1,N_REACTIONS
         IF (REACTION(NR)%SIMPLE_CHEMISTRY .AND. REACTION(NR)%N_SIMPLE_CHEMISTRY_REACTIONS > 0) THEN
            ! Unburned fuel
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + ZZ(II,JJ,KK,REACTION(NR)%FUEL_SMIX_INDEX)
            IF (REACTION(NR)%N_SIMPLE_CHEMISTRY_REACTIONS == 1) THEN
               ! Single step products
                GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + ZZ(II,JJ,KK,REACTION(NR)%PROD_SMIX_INDEX)/(1._EB+REACTION(NR)%S)
            ELSE
               ! Two step first intermediate products
                GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + ZZ(II,JJ,KK,REACTION(NR)%PROD_SMIX_INDEX)/(1._EB+REACTION(NR)%S)
                ! Two step second products
                GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + ZZ(II,JJ,KK,REACTION(REACTION(NR)%PAIR_INDEX)%PROD_SMIX_INDEX)/ &
                   ((1._EB+REACTION(NR)%S)*(1._EB+REACTION(REACTION(NR)%PAIR_INDEX)%S))
            ENDIF
         ENDIF
         IF (.NOT. REACTION(NR)%SIMPLE_CHEMISTRY) &
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + ZZ(II,JJ,KK,REACTION(NR)%FUEL_SMIX_INDEX) + &
                                   ZZ(II,JJ,KK,REACTION(NR)%PROD_SMIX_INDEX)/(1._EB+REACTION(NR)%S)
      ENDDO
   CASE(14)  ! DIVERGENCE
      GAS_PHASE_OUTPUT_RES = D(II,JJ,KK)
   CASE(15)  ! MIXING TIME
      GAS_PHASE_OUTPUT_RES = MIX_TIME(II,JJ,KK)
   CASE(16)  ! ABSORPTION COEFFICIENT
      III = MAX(1,MIN(II,IBAR))
      JJJ = MAX(1,MIN(JJ,JBAR))
      KKK = MAX(1,MIN(KK,KBAR))
      GAS_PHASE_OUTPUT_RES = KAPPA_GAS(III,JJJ,KKK)
   CASE(17)  ! VISCOSITY
      GAS_PHASE_OUTPUT_RES = MU(II,JJ,KK)
   CASE(18)  ! INTEGRATED INTENSITY
      GAS_PHASE_OUTPUT_RES = UII(II,JJ,KK)*0.001_EB
   CASE(19)  ! RADIATION LOSS
      GAS_PHASE_OUTPUT_RES = QR(II,JJ,KK)*0.001_EB
   CASE(20)  ! PARTICLE RADIATION LOSS
      IF (N_LP_ARRAY_INDICES>0) THEN
         GAS_PHASE_OUTPUT_RES = QR_W(II,JJ,KK)*0.001_EB
      ELSE
         GAS_PHASE_OUTPUT_RES = 0._EB
      ENDIF
   CASE(21)  ! RELATIVE HUMIDITY
      IF (H2O_INDEX<=0) THEN
         GAS_PHASE_OUTPUT_RES = 0._EB
      ELSE
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
         CALL GET_MASS_FRACTION(ZZ_GET,H2O_INDEX,Y_H2O)
         IF (H2O_SMIX_INDEX > 0) THEN
            IF (SPECIES_MIXTURE(H2O_SMIX_INDEX)%CONDENSATION_SMIX_INDEX > 0) &
               Y_H2O = Y_H2O - ZZ_GET(SPECIES_MIXTURE(H2O_SMIX_INDEX)%CONDENSATION_SMIX_INDEX)
         ENDIF
         GAS_PHASE_OUTPUT_RES = RELATIVE_HUMIDITY(Y_H2O,TMP(II,JJ,KK),PBAR(KK,PRESSURE_ZONE(II,JJ,KK)))
      ENDIF
   CASE(22)  ! HS
      GAS_PHASE_OUTPUT_RES = HS(II,JJ,KK)
   CASE(23)  ! KINETIC ENERGY (per unit mass) -- do not average because this operation is dissipative
      UU   = U(MIN(IBAR,II),JJ,KK)
      VV   = V(II,MIN(JBAR,JJ),KK)
      WW   = W(II,JJ,MIN(KBAR,KK))
      GAS_PHASE_OUTPUT_RES  = 0.5_EB*( UU**2 + VV**2 + WW**2 )

   CASE(24)  ! STRAIN RATE X
      III = MAX(1,MIN(II,IBAR))
      GAS_PHASE_OUTPUT_RES = (W(III,JJ+1,KK)-W(III,JJ,KK))*RDYN(JJ) + (V(III,JJ,KK+1)-V(III,JJ,KK))*RDZN(KK)
   CASE(25)  ! STRAIN RATE Y
      JJJ = MAX(1,MIN(JJ,JBAR))
      GAS_PHASE_OUTPUT_RES = (U(II,JJJ,KK+1)-U(II,JJJ,KK))*RDZN(KK) + (W(II+1,JJJ,KK)-W(II,JJJ,KK))*RDXN(II)
   CASE(26)  ! STRAIN RATE Z
      KKK = MAX(1,MIN(KK,KBAR))
      GAS_PHASE_OUTPUT_RES = (V(II+1,JJ,KKK)-V(II,JJ,KKK))*RDXN(II) + (U(II,JJ+1,KKK)-U(II,JJ,KKK))*RDYN(JJ)
   CASE(27)  ! VORTICITY X
      III = MAX(1,MIN(II,IBAR))
      GAS_PHASE_OUTPUT_RES = (W(III,JJ+1,KK)-W(III,JJ,KK))*RDYN(JJ) - (V(III,JJ,KK+1)-V(III,JJ,KK))*RDZN(KK)
   CASE(28)  ! VORTICITY Y
      JJJ = MAX(1,MIN(JJ,JBAR))
      GAS_PHASE_OUTPUT_RES = (U(II,JJJ,KK+1)-U(II,JJJ,KK))*RDZN(KK) - (W(II+1,JJJ,KK)-W(II,JJJ,KK))*RDXN(II)
   CASE(29)  ! VORTICITY Z
      KKK = MAX(1,MIN(KK,KBAR))
      GAS_PHASE_OUTPUT_RES = (V(II+1,JJ,KKK)-V(II,JJ,KKK))*RDXN(II) - (U(II,JJ+1,KKK)-U(II,JJ,KKK))*RDYN(JJ)

   CASE(30)  ! C_SMAG
      GAS_PHASE_OUTPUT_RES = 0._EB
      SELECT CASE (TURB_MODEL)
         CASE (CONSMAG,DYNSMAG)
            III = MAX(1,MIN(II,IBAR))
            JJJ = MAX(1,MIN(JJ,JBAR))
            KKK = MAX(1,MIN(KK,KBAR))
            DELTA = LES_FILTER_WIDTH(III,JJJ,KKK)
            GAS_PHASE_OUTPUT_RES = SQRT(CSD2(III,JJJ,KKK))/DELTA
      END SELECT
   CASE(31)  ! SPECIFIC HEAT
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CP,TMP(II,JJ,KK))
      GAS_PHASE_OUTPUT_RES = CP*0.001_EB

   CASE(32)  ! ORIENTED VELOCITY
      GAS_PHASE_OUTPUT_RES = U(II,JJ,KK)*ORIENTATION_VECTOR(1,DV%ORIENTATION_INDEX) + &
                             V(II,JJ,KK)*ORIENTATION_VECTOR(2,DV%ORIENTATION_INDEX) + &
                             W(II,JJ,KK)*ORIENTATION_VECTOR(3,DV%ORIENTATION_INDEX)

   CASE(33,50)  ! CONDUCTIVITY, MOLECULAR CONDUCTIVITY
      IF (SIM_MODE==DNS_MODE .OR. IND==50) THEN
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
         CALL GET_CONDUCTIVITY(ZZ_GET,GAS_PHASE_OUTPUT_RES,TMP(II,JJ,KK))
      ELSE
         GAS_PHASE_OUTPUT_RES = MU(II,JJ,KK)*CPOPR
      ENDIF

   CASE(34)  ! BACKGROUND PRESSURE
      GAS_PHASE_OUTPUT_RES = PBAR(KK,PRESSURE_ZONE(II,JJ,KK))

   CASE(35)  ! MOLECULAR WEIGHT
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_MOLECULAR_WEIGHT(ZZ_GET,GAS_PHASE_OUTPUT_RES)

   CASE(36)  ! POTENTIAL TEMPERATURE
      GAS_PHASE_OUTPUT_RES = GET_POTENTIAL_TEMPERATURE(TMP(II,JJ,KK),ZC(KK))

   CASE(37)  ! DIFFUSIVITY
      SELECT CASE (SIM_MODE)
         CASE DEFAULT
            GAS_PHASE_OUTPUT_RES = MU(II,JJ,KK)*RSC_T/RHO(II,JJ,KK)
         CASE (DNS_MODE)
            D_Z_N = D_Z(:,Z_INDEX)
            CALL INTERPOLATE1D_UNIFORM(LBOUND(D_Z_N,1),D_Z_N,TMP(II,JJ,KK),GAS_PHASE_OUTPUT_RES)
      END SELECT

   CASE(38)  ! RTE SOURCE CORRECTION FACTOR
      GAS_PHASE_OUTPUT_RES = RTE_SOURCE_CORRECTION_FACTOR
   CASE(39)  ! RAM (non-standard. You must uncomment GETPID in func.f90/SYSTEM_MEM_USAGE to use this quantity.)
      CALL SYSTEM_MEM_USAGE(RAM)
      GAS_PHASE_OUTPUT_RES = REAL(RAM,EB)/1000._EB
   CASE(40)  ! TIME
      GAS_PHASE_OUTPUT_RES = T_BEGIN + (T-T_BEGIN)*TIME_SHRINK_FACTOR
   CASE(41)  ! TIME STEP
      GAS_PHASE_OUTPUT_RES = DT
   CASE(42)  ! WALL CLOCK TIME
      GAS_PHASE_OUTPUT_RES = CURRENT_TIME() - WALL_CLOCK_START
   CASE(43)  ! WALL CLOCK TIME ITERATIONS
      IF (INITIALIZATION_PHASE) THEN
         GAS_PHASE_OUTPUT_RES = 0._EB
      ELSE
         GAS_PHASE_OUTPUT_RES = CURRENT_TIME() - WALL_CLOCK_START_ITERATIONS
      ENDIF
   CASE(44)  ! CPU TIME
      CALL CPU_TIME(CPUTIME)
      GAS_PHASE_OUTPUT_RES = CPUTIME - CPU_TIME_START
   CASE(45)  ! ITERATION
      GAS_PHASE_OUTPUT_RES = ICYC

   CASE(46:47)  ! SPECIFIC ENTHALPY and ENTHALPY
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_ENTHALPY(ZZ_GET,H_G,TMP(II,JJ,KK))
      IF (IND==46) GAS_PHASE_OUTPUT_RES = H_G*0.001_EB
      IF (IND==47) GAS_PHASE_OUTPUT_RES = RHO(II,JJ,KK)*H_G*0.001_EB

   CASE(48:49)  ! SPECIFIC SENSIBLE ENTHALPY and SENSIBLE ENTHALPY
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(II,JJ,KK))
      IF (IND==48) GAS_PHASE_OUTPUT_RES = H_S*0.001_EB
      IF (IND==49) GAS_PHASE_OUTPUT_RES = RHO(II,JJ,KK)*H_S*0.001_EB

   CASE(51)  ! RESOLVED KINETIC ENERGY (per unit mass)
      GAS_PHASE_OUTPUT_RES = KRES(II,JJ,KK)

   CASE(52)  ! WAVELET ERROR (wavelet error measure)
      GAS_PHASE_OUTPUT_RES = WAVELET_ERROR_MEASURE(II,JJ,KK,IND2,Y_INDEX,Z_INDEX,PART_INDEX,VELO_INDEX,DT,NM)

   CASE(53)  ! CELL U
      III = MAX(1,MIN(II,IBAR))
      GAS_PHASE_OUTPUT_RES = 0.5_EB*(U(III,JJ,KK)+U(MAX(1,III-1),JJ,KK))
   CASE(54)  ! CELL V
      JJJ = MAX(1,MIN(JJ,JBAR))
      GAS_PHASE_OUTPUT_RES = 0.5_EB*(V(II,JJJ,KK)+V(II,MAX(1,JJJ-1),KK))
   CASE(55)  ! CELL W
      KKK = MAX(1,MIN(KK,KBAR))
      GAS_PHASE_OUTPUT_RES = 0.5_EB*(W(II,JJ,KKK)+W(II,JJ,MAX(1,KKK-1)))

   CASE(56)  ! SUBGRID KINETIC ENERGY (per unit mass)
      DELTA = LES_FILTER_WIDTH(II,JJ,KK)
      IF (TEST_NEW_KSGS_MODEL) THEN
         GAS_PHASE_OUTPUT_RES = K_SGS_POPE(MU(II,JJ,KK)/RHO(II,JJ,KK),STRAIN_RATE(II,JJ,KK)**2,DELTA)
      ELSE
         GAS_PHASE_OUTPUT_RES = SUBGRID_KINETIC_ENERGY(MU(II,JJ,KK)-MU_DNS(II,JJ,KK),RHO(II,JJ,KK),C_DEARDORFF,DELTA)
      ENDIF

   CASE(57)  ! MAXIMUM VELOCITY ERROR
      GAS_PHASE_OUTPUT_RES = MAXVAL(VELOCITY_ERROR_MAX)

   CASE(58)  ! PRESSURE ITERATIONS
      GAS_PHASE_OUTPUT_RES = PRESSURE_ITERATIONS

   CASE(59)  ! OPEN NOZZLES
      GAS_PHASE_OUTPUT_RES = DEVC_PIPE_OPERATING(PIPE_INDEX)

   CASE(60)  ! ACTUATED SPRINKLERS
      GAS_PHASE_OUTPUT_RES = N_ACTUATED_SPRINKLERS

   CASE(61)  ! DRAG FORCE X
      GAS_PHASE_OUTPUT_RES = -0.5_EB*(RHO(II,JJ,KK)+RHO(II+1,JJ,KK))*FVX_D(II,JJ,KK)
   CASE(62)  ! DRAG FORCE Y
      GAS_PHASE_OUTPUT_RES = -0.5_EB*(RHO(II,JJ,KK)+RHO(II,JJ+1,KK))*FVY_D(II,JJ,KK)
   CASE(63)  ! DRAG FORCE Z
      GAS_PHASE_OUTPUT_RES = -0.5_EB*(RHO(II,JJ,KK)+RHO(II,JJ,KK+1))*FVZ_D(II,JJ,KK)

   CASE(64)  ! EFFECTIVE FLAME TEMPERATURE
      III = MAX(1,MIN(II,IBAR))
      JJJ = MAX(1,MIN(JJ,JBAR))
      KKK = MAX(1,MIN(KK,KBAR))
      IF (CHI_R(III,JJJ,KKK)*Q(II,JJ,KK)>QR_CLIP) THEN
         GAS_PHASE_OUTPUT_RES = TMP(II,JJ,KK)*RTE_SOURCE_CORRECTION_FACTOR**0.25_EB - TMPM
      ELSE
         GAS_PHASE_OUTPUT_RES = TMP(II,JJ,KK) - TMPM
      ENDIF

   CASE(68:69)  ! SPECIFIC INTERNAL ENERGY and INTERNAL ENERGY (per unit volume)
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_ENTHALPY(ZZ_GET,H_G,TMP(II,JJ,KK))
      IF (IND==68) GAS_PHASE_OUTPUT_RES = ( H_G - PBAR(KK,PRESSURE_ZONE(II,JJ,KK))/RHO(II,JJ,KK) )*0.001_EB
      IF (IND==69) GAS_PHASE_OUTPUT_RES = ( RHO(II,JJ,KK)*H_G - PBAR(KK,PRESSURE_ZONE(II,JJ,KK)) )*0.001_EB

   CASE(70)  ! CFL
      IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID) THEN
         GAS_PHASE_OUTPUT_RES = 0._EB
      ELSE
         IIM1 = MAX(II-1,0)
         JJM1 = MAX(JJ-1,0)
         KKM1 = MAX(KK-1,0)
         UODX = MAXVAL(ABS(US(IIM1:II,JJ,KK)))*RDX(II)
         VODY = MAXVAL(ABS(VS(II,JJM1:JJ,KK)))*RDY(JJ)
         WODZ = MAXVAL(ABS(WS(II,JJ,KKM1:KK)))*RDZ(KK)
         SELECT CASE (CFL_VELOCITY_NORM)
            CASE(0) ; UVW = MAX(UODX,VODY,WODZ) + ABS(DS(II,JJ,KK))
            CASE(1) ; UVW = UODX + VODY + WODZ  + ABS(DS(II,JJ,KK))
            CASE(2) ; UVW = SQRT(UODX**2+VODY**2+WODZ**2) + ABS(DS(II,JJ,KK))
            CASE(3) ; UVW = MAX(UODX,VODY,WODZ)
         END SELECT
         GAS_PHASE_OUTPUT_RES = DT*UVW
      ENDIF

   CASE(71)  ! VN
      IF (TWO_D) THEN
         R_DX2 = RDX(II)**2 + RDZ(KK)**2
      ELSE
         R_DX2 = RDX(II)**2 + RDY(JJ)**2 + RDZ(KK)**2
      ENDIF
      GAS_PHASE_OUTPUT_RES = DT*2._EB*R_DX2*MAX(D_Z_MAX(II,JJ,KK),MAX(RPR_T,RSC_T)*MU(II,JJ,KK)/RHO(II,JJ,KK))

   CASE(72)  ! CFL MAX
      GAS_PHASE_OUTPUT_RES = CFL
   CASE(73)  ! VN MAX
      GAS_PHASE_OUTPUT_RES = VN
   CASE(74)  ! POISSON ERROR
      GAS_PHASE_OUTPUT_RES = POIS_ERR
   CASE(75)  ! DIVERGENCE ERROR
      GAS_PHASE_OUTPUT_RES = RESMAX
   CASE(76)  ! RADIAL VELOCITY
      GAS_PHASE_OUTPUT_RES = ( XC(II)*0.5_EB*(U(II,JJ,KK)+U(II-1,JJ,KK)) + YC(JJ)*0.5_EB*(V(II,JJ,KK)+V(II,JJ-1,KK)) )/ &
                             SQRT(XC(II)**2+YC(JJ)**2)

   CASE(77)  ! LEVEL SET VALUE
      GAS_PHASE_OUTPUT_RES = PHI_LS(II,JJ)
   CASE(78)  ! RADIATION EMISSION
      GAS_PHASE_OUTPUT_RES = RADIATION_EMISSION(II,JJ,KK)*0.001_EB
   CASE(79)  ! RADIATION ABSORPTION
      GAS_PHASE_OUTPUT_RES = RADIATION_ABSORPTION(II,JJ,KK)*0.001_EB
   CASE(80)  ! CELL INDEX I
      GAS_PHASE_OUTPUT_RES = REAL(II,EB)
   CASE(81)  ! CELL INDEX J
      GAS_PHASE_OUTPUT_RES = REAL(JJ,EB)
   CASE(82)  ! CELL INDEX K
      GAS_PHASE_OUTPUT_RES = REAL(KK,EB)

   CASE(83)  ! Q CRITERION : Q = 1/2 (tr(Dij)^2 - tr(Dij^2))
      GAS_PHASE_OUTPUT_RES = 0._EB
      III=II; JJJ=JJ; KKK=KK
      IF (II == 0   ) III = II+1
      IF (II == IBP1) III = II-1
      IF (JJ == 0   ) JJJ = JJ+1
      IF (JJ == JBP1) JJJ = JJ-1
      IF (KK == 0   ) KKK = KK+1
      IF (KK == KBP1) KKK = KK-1
      IM1 = III-1
      JM1 = JJJ-1
      KM1 = KKK-1
      IIM1 = MAX(1,III-1)
      JJM1 = MAX(1,JJJ-1)
      KKM1 = MAX(1,KKK-1)
      IP1 = III+1
      JP1 = JJJ+1
      KP1 = KKK+1
      DUDX = RDX(III)*(U(III,JJJ,KKK)-U(IM1,JJJ,KKK))
      DUDY = 0.25_EB*RDY(JJJ)*(U(III,JP1,KKK)-U(III,JJM1,KKK)+U(IM1,JP1,KKK)-U(IM1,JJM1,KKK))
      DUDZ = 0.25_EB*RDZ(KKK)*(U(III,JJJ,KP1)-U(III,JJJ,KKM1)+U(IM1,JJJ,KP1)-U(IM1,JJJ,KKM1))
      DVDX = 0.25_EB*RDX(III)*(V(IP1,JJJ,KKK)-V(IIM1,JJJ,KKK)+V(IP1,JM1,KKK)-V(IIM1,JM1,KKK))
      DVDY = RDY(JJJ)*(V(III,JJJ,KKK)-V(III,JM1,KKK))
      DVDZ = 0.25_EB*RDZ(KKK)*(V(III,JJJ,KP1)-V(III,JJJ,KKM1)+V(III,JM1,KP1)-V(III,JM1,KKM1))
      DWDX = 0.25_EB*RDX(III)*(W(IP1,JJJ,KKK)-W(IIM1,JJJ,KKK)+W(IP1,JJJ,KM1)-W(IIM1,JJJ,KM1))
      DWDY = 0.25_EB*RDY(JJJ)*(W(III,JP1,KKK)-W(III,JJM1,KKK)+W(III,JP1,KM1)-W(III,JJM1,KM1))
      DWDZ = RDZ(KKK)*(W(III,JJJ,KKK)-W(III,JJJ,KM1))

      ! Q = 1/2 (tr(Dij)^2 - tr(Dij^2))
      GAS_PHASE_OUTPUT_RES = 0.5_EB*( (DUDX+DVDY+DWDZ)**2._EB            - &  ! tr(Dij)^2
                                      (DUDX*DUDX + DUDY*DVDX + DUDZ*DWDX + &  ! tr(Dij^2) = Dik*Dki
                                       DVDX*DUDY + DVDY*DVDY + DVDZ*DWDY + &
                                       DWDX*DUDZ + DWDY*DVDZ + DWDZ*DWDZ))
   CASE(84)  ! STRAIN RATE
      IM1 = MAX(0,II-1)
      JM1 = MAX(0,JJ-1)
      KM1 = MAX(0,KK-1)
      IIM1 = MAX(1,II-1)
      JJM1 = MAX(1,JJ-1)
      KKM1 = MAX(1,KK-1)
      IP1 = MIN(IBAR,II+1)
      JP1 = MIN(JBAR,JJ+1)
      KP1 = MIN(KBAR,KK+1)
      DUDX = RDX(II)*(U(II,JJ,KK)-U(IM1,JJ,KK))
      DVDY = RDY(JJ)*(V(II,JJ,KK)-V(II,JM1,KK))
      DWDZ = RDZ(KK)*(W(II,JJ,KK)-W(II,JJ,KM1))
      ONTHDIV = ONTH*(DUDX+DVDY+DWDZ)
      S11 = DUDX - ONTHDIV
      S22 = DVDY - ONTHDIV
      S33 = DWDZ - ONTHDIV
      DUDY = 0.25_EB*RDY(JJ)*(U(II,JP1,KK)-U(II,JJM1,KK)+U(IM1,JP1,KK)-U(IM1,JJM1,KK))
      DUDZ = 0.25_EB*RDZ(KK)*(U(II,JJ,KP1)-U(II,JJ,KKM1)+U(IM1,JJ,KP1)-U(IM1,JJ,KKM1))
      DVDX = 0.25_EB*RDX(II)*(V(IP1,JJ,KK)-V(IIM1,JJ,KK)+V(IP1,JM1,KK)-V(IIM1,JM1,KK))
      DVDZ = 0.25_EB*RDZ(KK)*(V(II,JJ,KP1)-V(II,JJ,KKM1)+V(II,JM1,KP1)-V(II,JM1,KKM1))
      DWDX = 0.25_EB*RDX(II)*(W(IP1,JJ,KK)-W(IIM1,JJ,KK)+W(IP1,JJ,KM1)-W(IIM1,JJ,KM1))
      DWDY = 0.25_EB*RDY(JJ)*(W(II,JP1,KK)-W(II,JJM1,KK)+W(II,JP1,KM1)-W(II,JJM1,KM1))
      S12 = 0.5_EB*(DUDY+DVDX)
      S13 = 0.5_EB*(DUDZ+DWDX)
      S23 = 0.5_EB*(DVDZ+DWDY)
      GAS_PHASE_OUTPUT_RES = SQRT(2._EB*(S11**2 + S22**2 + S33**2 + 2._EB*(S12**2 + S13**2 + S23**2)))
   CASE(85)  ! KOLMOGOROV LENGTH SCALE
      SS = GAS_PHASE_OUTPUT(T,DT,NM,II,JJ,KK,84,IND2,Y_INDEX,Z_INDEX,ELEM_INDX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,&
                            REAC_INDEX,MATL_INDEX)
      DISSIPATION_RATE = MU(II,JJ,KK)/RHO(II,JJ,KK)*SS**2
      GAS_PHASE_OUTPUT_RES = ((MU_DNS(II,JJ,KK)/RHO(II,JJ,KK))**3/(DISSIPATION_RATE+EPS))**0.25_EB
   CASE(86)  ! CELL REYNOLDS NUMBER
      III = MAX(1,MIN(II,IBAR))
      JJJ = MAX(1,MIN(JJ,JBAR))
      KKK = MAX(1,MIN(KK,KBAR))
      DELTA = LES_FILTER_WIDTH(III,JJJ,KKK)
      ETA = GAS_PHASE_OUTPUT(T,DT,NM,II,JJ,KK,85,IND2,Y_INDEX,Z_INDEX,ELEM_INDX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,&
                             REAC_INDEX,MATL_INDEX)
      GAS_PHASE_OUTPUT_RES = DELTA/(ETA+EPS)
   CASE(87)  ! MOLECULAR VISCOSITY
      GAS_PHASE_OUTPUT_RES = MU_DNS(II,JJ,KK)
   CASE(88)  ! DISSIPATION RATE
      SS = GAS_PHASE_OUTPUT(T,DT,NM,II,JJ,KK,84,IND2,Y_INDEX,Z_INDEX,ELEM_INDX,PART_INDEX,VELO_INDEX,PIPE_INDEX,PROP_INDEX,&
                            REAC_INDEX,MATL_INDEX)
      GAS_PHASE_OUTPUT_RES = MU(II,JJ,KK)/RHO(II,JJ,KK)*SS**2
   CASE(89)  ! KINEMATIC VISCOSITY
      GAS_PHASE_OUTPUT_RES = MU(II,JJ,KK)/RHO(II,JJ,KK)
   CASE(90)  ! MASS FRACTION
      GAS_PHASE_OUTPUT_RES = Y_SPECIES/(1._EB-Y_H2O)

   CASE(91:93) ! MASS FLUX
      IP=II ; JP=JJ ; KP=KK
      SELECT CASE(IND)
         CASE(91) ; IP=II+1 ; VEL=U(II,JJ,KK)  ! MASS FLUX X
         CASE(92) ; JP=JJ+1 ; VEL=V(II,JJ,KK)  ! MASS FLUX Y
         CASE(93) ; KP=KK+1 ; VEL=W(II,JJ,KK)  ! MASS FLUX Z
      END SELECT
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IP,JP,KP,1:N_TRACKED_SPECIES)
      Y_SPECIES2 = 1.0_EB
      IF (Z_INDEX > 0) THEN
         Y_SPECIES2 = ZZ_GET(Z_INDEX)
      ELSEIF (Y_INDEX > 0) THEN
         CALL GET_MASS_FRACTION(ZZ_GET,Y_INDEX,Y_SPECIES2)
      ENDIF
      GAS_PHASE_OUTPUT_RES = 0.5_EB*(RHO(II,JJ,KK)*Y_SPECIES+RHO(IP,JP,KP)*Y_SPECIES2)*VEL

   CASE(94)  ! VOLUME FRACTION
      GAS_PHASE_OUTPUT_RES =  RCON*Y_SPECIES/RSUM(II,JJ,KK)/(1._EB-R_Y_H2O/RSUM(II,JJ,KK))
   CASE(95)  ! VISIBILITY
      IF (Z_INDEX>0) THEN
         MEC = SPECIES_MIXTURE(Z_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (Y_INDEX>0) THEN
         MEC = SPECIES(Y_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ENDIF
      EXT_COEF = Y_SPECIES*RHO(II,JJ,KK)*MEC
      GAS_PHASE_OUTPUT_RES = VISIBILITY_FACTOR/MAX(EC_LL,EXT_COEF)
   CASE(96)  ! AEROSOL VOLUME FRACTION
      IF (Z_INDEX >0) THEN
         GAS_PHASE_OUTPUT_RES = Y_SPECIES*RHO(II,JJ,KK)/SPECIES(SPECIES_MIXTURE(Z_INDEX)%SINGLE_SPEC_INDEX)%DENSITY_SOLID
      ELSEIF (Y_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = Y_SPECIES*RHO(II,JJ,KK)/SPECIES(Y_INDEX)%DENSITY_SOLID
      ENDIF
   CASE(97)  ! EXTINCTION COEFFICIENT
      IF (Z_INDEX>0) THEN
         MEC = SPECIES_MIXTURE(Z_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (Y_INDEX>0) THEN
         MEC = SPECIES(Y_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ENDIF
      EXT_COEF = Y_SPECIES*RHO(II,JJ,KK)*MEC
      GAS_PHASE_OUTPUT_RES = Y_SPECIES*RHO(II,JJ,KK)*MEC
   CASE(98)  ! OPTICAL DENSITY
      IF (Z_INDEX>0) THEN
         MEC = SPECIES_MIXTURE(Z_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (Y_INDEX>0) THEN
         MEC = SPECIES(Y_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ENDIF
      GAS_PHASE_OUTPUT_RES = Y_SPECIES*RHO(II,JJ,KK)*MEC/2.3_EB

   CASE(99)  ! PRESSURE POISSON RESIDUAL
      GAS_PHASE_OUTPUT_RES = PP_RESIDUAL(II,JJ,KK)
   CASE(100) ! PRESSURE ZONE
      GAS_PHASE_OUTPUT_RES = PRESSURE_ZONE(II,JJ,KK)

   CASE(101)  ! FIC
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      GAS_PHASE_OUTPUT_RES = FIC(ZZ_GET,RSUM(II,JJ,KK))

   CASE(102)  ! BULK DENSITY
      IC = CELL_INDEX(II,JJ,KK)
      IF (.NOT.CELL(IC)%SOLID .OR. CELL(IC)%OBST_INDEX<1) THEN
         GAS_PHASE_OUTPUT_RES = 0._EB
      ELSE
         GAS_PHASE_OUTPUT_RES = OBSTRUCTION(CELL(IC)%OBST_INDEX)%MASS*RDX(II)*RRN(II)*RDY(JJ)*RDZ(KK)
      ENDIF

   CASE(105:107) ! Hot Gas Layer Reduction
      CALL GET_LAYER_HEIGHT_INTEGRALS(SDV%I1,SDV%J1,SDV%K1,SDV%K2,DV%Z_INT,DV%Z1,SDV%VALUE_1,SDV%VALUE_2,SDV%VALUE_3,&
                                      SDV%VALUE_4,DV%TMP_LOW)
      GAS_PHASE_OUTPUT_RES = SDV%VALUE_1

   CASE(109)  ! FED
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      IF (PROP_INDEX>0) THEN
         FED_ACTIVITY = PROPERTY(PROP_INDEX)%FED_ACTIVITY
      ELSE
         FED_ACTIVITY = 2
      ENDIF
      GAS_PHASE_OUTPUT_RES = FED(ZZ_GET,RSUM(II,JJ,KK),FED_ACTIVITY)

   CASE(110)  ! THERMOCOUPLE
      IF (T > T_BEGIN) THEN
         TMP_G = TMP(II,JJ,KK)
         IF (PY%HEAT_TRANSFER_COEFFICIENT<0._EB) THEN
            UU      = U(II,JJ,KK)
            VV      = V(II,JJ,KK)
            WW      = W(II,JJ,KK)
            VEL2    = UU**2+VV**2+WW**2
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
            CALL GET_VISCOSITY(ZZ_GET,MU_G,TMP(II,JJ,KK))
            CALL GET_CONDUCTIVITY(ZZ_GET,K_G,TMP(II,JJ,KK))
            RE_D    = RHO(II,JJ,KK)*SQRT(VEL2)*PY%DIAMETER/MU_G
            CALL FORCED_CONVECTION_MODEL(NUSSELT,RE_D,PR_ONTH,SURF_SPHERICAL)
            H_TC    = NUSSELT*K_G/PY%DIAMETER
         ELSE
            H_TC    = PY%HEAT_TRANSFER_COEFFICIENT
         ENDIF
         FAC = 6._EB/(PY%DENSITY*PY%SPECIFIC_HEAT(NINT(DV%TMP_L))*PY%DIAMETER)*DT
         DV%TMP_L = (DV%TMP_L + FAC*(H_TC*(TMP_G-0.5_EB*DV%TMP_L) + &
                     PY%EMISSIVITY*(0.25_EB*UII(II,JJ,KK)+SIGMA*DV%TMP_L**4))) / &
                    (1._EB + FAC*(0.5_EB*H_TC+2._EB*PY%EMISSIVITY*SIGMA*DV%TMP_L**3))
      ENDIF
      GAS_PHASE_OUTPUT_RES = DV%TMP_L - TMPM

   CASE(111:113)  ! ENTHALPY FLUX
      IP=II ; JP=JJ ; KP=KK
      SELECT CASE(IND)
         CASE(111) ; IP=II+1 ; VEL=U(II,JJ,KK) ; R_DN=RDXN(II)  ! ENTHALPY FLUX X
         CASE(112) ; JP=JJ+1 ; VEL=V(II,JJ,KK) ; R_DN=RDYN(JJ)  ! ENTHALPY FLUX Y
         CASE(113) ; KP=KK+1 ; VEL=W(II,JJ,KK) ; R_DN=RDZN(KK)  ! ENTHALPY FLUX Z
      END SELECT
      TMP_TC = 0.5_EB*(TMP(II,JJ,KK)+TMP(IP,JP,KP))
      ZZ_GET(1:N_TRACKED_SPECIES) = 0.5_EB*(ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)+ZZ(IP,JP,KP,1:N_TRACKED_SPECIES))
      CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_G_SUM,TMP_TC)
      CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_G,TMPA)
      GAS_PHASE_OUTPUT_RES = VEL*0.5_EB*(RHO(II,JJ,KK)+RHO(IP,JP,KP))*(H_G_SUM-H_G)
      IF (SIM_MODE==DNS_MODE) THEN
         CALL GET_CONDUCTIVITY(ZZ_GET,K_G,TMP(II,JJ,KK))
      ELSE
         K_G = MU(II,JJ,KK)*CPOPR
      ENDIF
      GAS_PHASE_OUTPUT_RES = (GAS_PHASE_OUTPUT_RES - K_G*(TMP(IP,JP,KP)-TMP(II,JJ,KK))*R_DN)*0.001

   CASE(114) ! BI-DIRECTIONAL PROBE
      ! Fits taken from
      ! McCaffrey and Heskestad, A Robust Bidirectional Low-Velocity Probe for Flame and Fire Application
      ! Combustion and Flame, 26, 125 - 127, (1976).
      IF (PY%TC) THEN
         PROBE_TMP = GAS_PHASE_OUTPUT(T,DT,NM,II,JJ,KK,110,IND2,Y_INDEX,Z_INDEX,ELEM_INDX,PART_INDEX,VELO_INDEX,PIPE_INDEX,&
                                           PROP_INDEX,REAC_INDEX,MATL_INDEX,ICC_IN,JCC_IN) + TMPM
      ELSE
         PROBE_TMP = TMP(II,JJ,KK)
      ENDIF
      PROBE_RHO = MW_AIR*P_STP/(R0*PROBE_TMP)
      UU = 0.5_EB*(U(MAX(0,II-1),JJ,KK)+U(MIN(IBAR,II),JJ,KK))
      VV = 0.5_EB*(V(II,MAX(0,JJ-1),KK)+V(II,MIN(JBAR,JJ),KK))
      WW = 0.5_EB*(W(II,JJ,MAX(0,KK-1))+W(II,JJ,MIN(KBAR,KK)))
      VEL2 = UU**2+VV**2+WW**2
      VEL = SQRT(VEL2)
      ! Adjust for effect of flow direction on measured pressure
      COSTHETA = (UU*ORIENTATION_VECTOR(1,DV%ORIENTATION_INDEX)+VV*ORIENTATION_VECTOR(2,DV%ORIENTATION_INDEX)+ &
                  WW*ORIENTATION_VECTOR(3,DV%ORIENTATION_INDEX))/VEL
      FAC = MAX(0._EB,-2.308_EB*ABS(COSTHETA)**3 + 2.533_EB*ABS(COSTHETA)**2 + 0.7847_EB*ABS(COSTHETA) - 0.0097_EB)
      VEL = FAC*VEL
      ! Adjust for effect of Re number on measured pressure
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_VISCOSITY(ZZ_GET,MU_G,TMP(II,JJ,KK))
      RE_D = MIN(3800._EB,MAX(40._EB,RHO(II,JJ,KK)*VEL*PY%PROBE_DIAMETER/MU_G))
      FAC = 1.533_EB-0.001366_EB*RE_D+0.000001688_EB*RE_D**2-0.0000000009706_EB*RE_D**3+&
            0.0000000000002555_EB*RE_D**4-2.484E-17_EB*RE_D**5
      IF (PY%CALIBRATION_CONSTANT > 0._EB) THEN
         GAS_PHASE_OUTPUT_RES = SIGN(1._EB,COSTHETA)*VEL*PY%CALIBRATION_CONSTANT*FAC*SQRT(RHO(II,JJ,KK)/PROBE_RHO)
      ELSE
         PROBE_DELTA_P = (VEL*FAC)**2*RHO(II,JJ,KK)*0.5_EB
         ! LJ AIR viscosity fit 100 K to 5000 K
         PROBE_TMP = MIN(5000._EB,MAX(100._EB,PROBE_TMP))
         MU_G = 1.5205E-22_EB*PROBE_TMP**5 - 2.1417E-18_EB*PROBE_TMP**4 + 1.1402E-14_EB*PROBE_TMP**3 - &
                  2.9846E-11_EB*PROBE_TMP**2 + 5.9898E-8_EB*PROBE_TMP + 0.000002352_EB
         JJJ = 1
         BP_LOOP: DO
            RE_D = MIN(3800._EB,MAX(40._EB,PROBE_RHO*VEL*PY%PROBE_DIAMETER/MU_G))
            FAC = 1.533_EB-0.001366_EB*RE_D+0.000001688_EB*RE_D**2-0.0000000009706_EB*RE_D**3+&
                  0.0000000000002555_EB*RE_D**4-2.484E-17_EB*RE_D**5
            GAS_PHASE_OUTPUT_RES = SIGN(1._EB,COSTHETA)*1._EB/FAC*SQRT(2._EB*PROBE_DELTA_P/PROBE_RHO)
            IF (JJJ > 9 .OR. ABS(VEL-GAS_PHASE_OUTPUT_RES)/(GAS_PHASE_OUTPUT_RES+TWENTY_EPSILON_EB) < 0.001_EB) EXIT BP_LOOP
            VEL = 0.2_EB*VEL+0.8_EB*GAS_PHASE_OUTPUT_RES
            JJJ = JJJ + 1
         ENDDO BP_LOOP
      ENDIF

   CASE(130) ! EXTINCTION
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      ZZ_FUEL = 0._EB
      ZZ_OX = 0._EB
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NR=1,N_REACTIONS
         DO NS=1,N_TRACKED_SPECIES
            IF (REACTION(NR)%NU(NS) < 0._EB) THEN
               IF (NS == REACTION(NR)%FUEL_SMIX_INDEX) ZZ_FUEL = ZZ_FUEL + ZZ_GET(NS)
               IF (NS /= REACTION(NR)%FUEL_SMIX_INDEX .AND. NR == 1) ZZ_OX = ZZ_GET(NS)
            ENDIF
         ENDDO
      ENDDO
      IF (ZZ_FUEL<TWENTY_EPSILON_EB .OR. ZZ_OX<TWENTY_EPSILON_EB .OR. Q(II,JJ,KK)<TWENTY_EPSILON_EB) GAS_PHASE_OUTPUT_RES = -1._EB
      IF (ZZ_FUEL>ZZ_MIN_GLOBAL .AND. ZZ_OX > ZZ_MIN_GLOBAL .AND. Q(II,JJ,KK) < TWENTY_EPSILON_EB) GAS_PHASE_OUTPUT_RES = 1._EB

   CASE(131) ! CHEMISTRY SUBITERATIONS
      GAS_PHASE_OUTPUT_RES = CHEM_SUBIT(II,JJ,KK)

   CASE(132) ! REAC SOURCE TERM
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = REAC_SOURCE_TERM(II,JJ,KK,Z_INDEX)
      ELSEIF (Y_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DOT_PRODUCT(Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),REAC_SOURCE_TERM(II,JJ,KK,1:N_TRACKED_SPECIES))
      ENDIF

   CASE(133) ! SUM LUMPED MASS FRACTIONS
      GAS_PHASE_OUTPUT_RES = SUM(ZZ(II,JJ,KK,1:N_TRACKED_SPECIES))

   CASE(134) ! SUM PRIMITIVE MASS FRACTIONS
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_MASS_FRACTION_ALL(ZZ_GET,Y_ALL)
      GAS_PHASE_OUTPUT_RES = SUM(Y_ALL)

   CASE(135) ! MACH NUMBER
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CP,TMP(II,JJ,KK))
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RCON)
      GAMMA_LOC = CP/(CP-RCON)
      GAS_PHASE_OUTPUT_RES = SQRT(2._EB*KRES(II,JJ,KK))/SQRT(RCON*TMP(II,JJ,KK)*GAMMA_LOC)

   CASE(136) ! UNMIXED FRACTION
      GAS_PHASE_OUTPUT_RES = INITIAL_UNMIXED_FRACTION*EXP(-DT/MIX_TIME(II,JJ,KK))

   CASE(138) ! HRRPUV REAC
      GAS_PHASE_OUTPUT_RES = Q_REAC(II,JJ,KK,REAC_INDEX)*0.001_EB

   CASE(140) ! FVX_B
      GAS_PHASE_OUTPUT_RES = FVX_B(II,JJ,KK)
   CASE(141) ! FVY_B
      GAS_PHASE_OUTPUT_RES = FVY_B(II,JJ,KK)
   CASE(142) ! FVZ_B
      GAS_PHASE_OUTPUT_RES = FVZ_B(II,JJ,KK)

   CASE(143) ! COMBUSTION EFFICIENCY
      IF (Q(II,JJ,KK)>TWENTY_EPSILON_EB) THEN
         GAS_PHASE_OUTPUT_RES = MIN(DT/MIX_TIME(II,JJ,KK),1._EB)
      ELSE
         GAS_PHASE_OUTPUT_RES = 0._EB
      ENDIF
   CASE(144)  ! ELEMENT MASS FRACTION
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NS=1,N_TRACKED_SPECIES
         GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + &
             ZZ_GET(NS)*SPECIES_MIXTURE(NS)%ATOMS(ELEM_INDX)*ELEMENT(ELEM_INDX)%MASS/SPECIES_MIXTURE(NS)%MW
      ENDDO
   CASE(145)  ! EQUIVALENCE RATIO
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      GAS_PHASE_OUTPUT_RES = 0._EB
      CALL CALC_EQUIV_RATIO(ZZ_GET(1:N_TRACKED_SPECIES), GAS_PHASE_OUTPUT_RES)
   CASE(150) ! SUM LUMPED VOLUME FRACTIONS
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_MOLECULAR_WEIGHT(ZZ_GET,MW)
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO N=1,N_TRACKED_SPECIES
         GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + ZZ(II,JJ,KK,N)/SPECIES_MIXTURE(N)%MW*MW
      ENDDO

   CASE(153) ! NOZZLE FLOW RATE
      GAS_PHASE_OUTPUT_RES = PY%FLOW_RATE

   CASE(154:155) ! TRANSMISSION, PATH OBSCURATION
      EXT_COEF   = 0._EB
      IF (PY%Y_INDEX>0) THEN
         MASS_EXT_COEF = SPECIES(PY%Y_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (PY%Z_INDEX>0) THEN
         MASS_EXT_COEF = SPECIES_MIXTURE(PY%Z_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (SOOT_INDEX>0) THEN
         MASS_EXT_COEF = SPECIES(SOOT_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSE
         MASS_EXT_COEF = 0._EB
      ENDIF
      DO NN=1,SDV%N_PATH
         I = SDV%I_PATH(NN)
         J = SDV%J_PATH(NN)
         K = SDV%K_PATH(NN)
         IF (PY%Y_INDEX>0) THEN
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_MASS_FRACTION(ZZ_GET,PY%Y_INDEX,Y_MF_INT)
         ELSEIF (PY%Z_INDEX>0) THEN
            Y_MF_INT = ZZ(I,J,K,PY%Z_INDEX)
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_MASS_FRACTION(ZZ_GET,SOOT_INDEX,Y_MF_INT)
         ENDIF
         EXT_COEF = EXT_COEF + Y_MF_INT*RHO(I,J,K)*SDV%D_PATH(NN)
      ENDDO
      GAS_PHASE_OUTPUT_RES = MASS_EXT_COEF*EXT_COEF  ! This output is only a component of the actual output QUANTITY

   CASE(156) ! SPRINKLER LINK TEMPERATURE
      I = DV%I(1)
      J = DV%J(1)
      K = DV%K(1)
      TMP_G = TMP(I,J,K)
      VEL2  = 0.25_EB*( (U(I,J,K)+U(I-1,J,K))**2 +(V(I,J,K)+V(I,J-1,K))**2 + (W(I,J,K)+W(I,J,K-1))**2 )
      VEL   = SQRT(VEL2)
      VELSR = SQRT(VEL)
      WATER_VOL_FRAC = 0._EB
      IF (H2O_INDEX > 0) THEN
         DO NN = 1,N_LAGRANGIAN_CLASSES
            IF (LAGRANGIAN_PARTICLE_CLASS(NN)%Y_INDEX==H2O_INDEX) WATER_VOL_FRAC = WATER_VOL_FRAC + &
               AVG_DROP_DEN(I,J,K,LAGRANGIAN_PARTICLE_CLASS(NN)%ARRAY_INDEX)/LAGRANGIAN_PARTICLE_CLASS(NN)%DENSITY
         ENDDO
      ENDIF
      RHS      = (VELSR*(TMP_G-DV%TMP_L) - PY%C_FACTOR*(DV%TMP_L-PY%INITIAL_TEMPERATURE) - C_DIMARZO*VEL*WATER_VOL_FRAC)/PY%RTI
      DV%TMP_L = MAX(MIN(TMP_G,PY%INITIAL_TEMPERATURE) , DV%TMP_L + DT*RHS)
      GAS_PHASE_OUTPUT_RES = DV%TMP_L - TMPM

   CASE(157) ! LINK TEMPERATURE
      I = DV%I(1)
      J = DV%J(1)
      K = DV%K(1)
      TMP_G = TMP(I,J,K)
      VEL2  = 0.25_EB*( (U(I,J,K)+U(I-1,J,K))**2 + (V(I,J,K)+V(I,J-1,K))**2 + (W(I,J,K)+W(I,J,K-1))**2 )
      VEL   = SQRT(VEL2)
      VELSR = SQRT(VEL)
      DV%TMP_L  = DV%TMP_L + DT*VELSR*(TMP_G-DV%TMP_L)/PY%RTI
      GAS_PHASE_OUTPUT_RES       = DV%TMP_L - TMPM

   CASE(158) ! CHAMBER OBSCURATION
      IF (Y_INDEX > 0) THEN
         MASS_EXT_COEF = SPECIES(Y_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (Z_INDEX>0) THEN
         MASS_EXT_COEF = SPECIES_MIXTURE(Z_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (SOOT_INDEX>0) THEN
         MASS_EXT_COEF = SPECIES(SOOT_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSE
         MASS_EXT_COEF = 0._EB
      ENDIF
      I = DV%I(1)
      J = DV%J(1)
      K = DV%K(1)
      VEL2 = 0.25_EB*( (U(I,J,K)+U(I-1,J,K))**2 + (V(I,J,K)+V(I,J-1,K))**2 + (W(I,J,K)+W(I,J,K-1))**2 )
      VEL  = MAX(SQRT(VEL2),1.0E-10_EB)
      IF (DV%N_T_E>=UBOUND(DV%T_E,1)) THEN
         DV%T_E => REALLOCATE(DV%T_E,0,DV%N_T_E+1000)
         DV%Y_E => REALLOCATE(DV%Y_E,0,DV%N_T_E+1000)
      ENDIF
      DV%N_T_E = DV%N_T_E + 1
      DV%Y_E(DV%N_T_E) = Y_SPECIES
      DV%T_E(DV%N_T_E) = T
      DT_C = PY%ALPHA_C*VEL**PY%BETA_C
      DT_E = PY%ALPHA_E*VEL**PY%BETA_E
      Y_E_LAG = 0._EB
      LAG_LOOP: DO IL=DV%N_T_E-1,0,-1
         IF (DV%T_E(IL) > T-DT_E) CYCLE LAG_LOOP
         T_RATIO = (T-DT_E-DV%T_E(IL))/(DV%T_E(IL+1)-DV%T_E(IL))
         Y_E_LAG = MAX(0._EB,DV%Y_E(IL) + T_RATIO*(DV%Y_E(IL+1)-DV%Y_E(IL)))
         EXIT LAG_LOOP
      ENDDO LAG_LOOP
      DV%Y_C = MAX(0._EB,DV%Y_C + DT*(Y_E_LAG - DV%Y_C)/DT_C)
      GAS_PHASE_OUTPUT_RES = (1._EB-EXP(-MASS_EXT_COEF*RHO(I,J,K)*DV%Y_C))*100._EB  ! Obscuration

   CASE(159) ! CONTROL VALUE
      GAS_PHASE_OUTPUT_RES = CONTROL(DV%CTRL_INDEX)%INSTANT_VALUE

   CASE(160) ! CONTROL
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (CONTROL(DV%CTRL_INDEX)%CURRENT_STATE) GAS_PHASE_OUTPUT_RES = 1._EB

   CASE(161) ! ASPIRATION
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (T >= DV%T) THEN
         ! If enough time has passed shift soot density array
         DV%T = T + DV%DT
         DV%TIME_ARRAY(0:99) = DV%TIME_ARRAY(1:100)
         DV%YY_SOOT(:,0:99) = DV%YY_SOOT(:,1:100)
         DV%YY_SOOT(:,100) = 0._EB
      ENDIF
      DV%TIME_ARRAY(100) = T
      DO N = 1, DV%N_INPUTS
         ! Update soot density array
         DV2 => DEVICE(DV%DEVC_INDEX(N))
         IF (ABS(DV%T - T - DV%DT)<=SPACING(DV%T)) THEN
            DV%YY_SOOT(N,100) = DV2%INSTANT_VALUE
         ELSE
            DV%YY_SOOT(N,100) = (DV%YY_SOOT(N,100) * (T - DV%TIME_ARRAY(99) - DT) +  DT * DV2%INSTANT_VALUE) / &
                                (T - DV%TIME_ARRAY(99))
         END IF
         ! Sum soot densities weighted by flow rate
         CALL INTERPOLATE1D(DV%TIME_ARRAY,DV%YY_SOOT(N,:),T-DV2%DELAY,Y_SPECIES)
         GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + DV2%FLOWRATE * Y_SPECIES
      ENDDO
      ! Complete weighting and compute % obs
      GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES / DV%TOTAL_FLOWRATE
      IF (DV2%Y_INDEX > 0) THEN
         MASS_EXT_COEF = SPECIES(DV2%Y_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (DV2%Z_INDEX>0) THEN
         MASS_EXT_COEF = SPECIES_MIXTURE(DV2%Z_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSEIF (SOOT_INDEX>0) THEN
         MASS_EXT_COEF = SPECIES(SOOT_INDEX)%MASS_EXTINCTION_COEFFICIENT
      ELSE
         MASS_EXT_COEF = 0._EB
      ENDIF
      GAS_PHASE_OUTPUT_RES = (1._EB-EXP(-MASS_EXT_COEF*GAS_PHASE_OUTPUT_RES))*100._EB  ! Obscuration

   CASE(163) ! PATHLENGTH
      PATHLENGTH = 0._EB
      DO NN=1,SDV%N_PATH
         PATHLENGTH = PATHLENGTH + SDV%D_PATH(NN)
      ENDDO
      GAS_PHASE_OUTPUT_RES = PATHLENGTH

   CASE(164) ! FIRE DEPTH
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NN=1,SDV%N_PATH
         I = SDV%I_PATH(NN)
         J = SDV%J_PATH(NN)
         K = SDV%K_PATH(NN)
         IF (Q(I,J,K)>(1.E3_EB*DV%SETPOINT)) THEN
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + SDV%D_PATH(NN)
         ENDIF
      ENDDO

   CASE(170) ! MPUV
      LPC => LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)
      GAS_PHASE_OUTPUT_RES = AVG_DROP_DEN(II,JJ,KK,LPC%ARRAY_INDEX)

   CASE(171) ! ADD
      LPC => LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)
      GAS_PHASE_OUTPUT_RES = AVG_DROP_RAD(II,JJ,KK,LPC%ARRAY_INDEX)*2.E6_EB

   CASE(172) ! ADT
      LPC => LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)
      GAS_PHASE_OUTPUT_RES = AVG_DROP_TMP(II,JJ,KK,LPC%ARRAY_INDEX) - TMPM

   CASE(173) ! ADA
      LPC => LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)
      GAS_PHASE_OUTPUT_RES = AVG_DROP_AREA(II,JJ,KK,LPC%ARRAY_INDEX)

   CASE(174) ! QABS
      GAS_PHASE_OUTPUT_RES = 0._EB
      LPC => LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)
      IF (ABS(AVG_DROP_AREA(II,JJ,KK,LPC%ARRAY_INDEX))>TWENTY_EPSILON_EB) THEN
         DO N = 1,NUMBER_SPECTRAL_BANDS
            IF (NUMBER_SPECTRAL_BANDS==1) THEN
               BBF = 1._EB
            ELSE
               BBF = BLACKBODY_FRACTION(WL_LOW(N),WL_HIGH(N),RADTMP)
            ENDIF
            CALL INTERPOLATE1D(LPC%R50,LPC%WQABS(:,N),AVG_DROP_RAD(II,JJ,KK,LPC%ARRAY_INDEX),Q_SUM)
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + BBF*Q_SUM
         ENDDO
         GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES/REAL(NUMBER_SPECTRAL_BANDS,EB)
      ENDIF

   CASE(175) ! QSCA
      GAS_PHASE_OUTPUT_RES = 0._EB
      LPC => LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)
      IF (ABS(AVG_DROP_AREA(II,JJ,KK,LPC%ARRAY_INDEX))>TWENTY_EPSILON_EB) THEN
         DO N = 1,NUMBER_SPECTRAL_BANDS
            IF (NUMBER_SPECTRAL_BANDS==1) THEN
               BBF = 1._EB
            ELSE
               BBF = BLACKBODY_FRACTION(WL_LOW(N),WL_HIGH(N),RADTMP)
            ENDIF
            CALL INTERPOLATE1D(LPC%R50,LPC%WQSCA(:,N),AVG_DROP_RAD(II,JJ,KK,LPC%ARRAY_INDEX),Q_SUM)
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + BBF*Q_SUM
         ENDDO
         GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES/REAL(NUMBER_SPECTRAL_BANDS,EB)
      ENDIF

   CASE(176) ! PARTICLE FLUX X
      GAS_PHASE_OUTPUT_RES = WFX(II,JJ,KK)

   CASE(177) ! PARTICLE FLUX Y
      GAS_PHASE_OUTPUT_RES = WFY(II,JJ,KK)

   CASE(178) ! PARTICLE FLUX Z
      GAS_PHASE_OUTPUT_RES = WFZ(II,JJ,KK)

   CASE(179) ! MPUV_Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NN = 1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(NN)
         IF (LPC%LIQUID_DROPLET .AND. LPC%Y_INDEX==Y_INDEX) &
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + AVG_DROP_DEN(II,JJ,KK,LPC%ARRAY_INDEX)
      ENDDO

   CASE(180) ! ADD_Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NN = 1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(NN)
         IF (LPC%LIQUID_DROPLET .AND. LPC%Y_INDEX==Y_INDEX) &
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + AVG_DROP_RAD(II,JJ,KK,LAGRANGIAN_PARTICLE_CLASS(NN)%ARRAY_INDEX)
      ENDDO

   CASE(181) ! ADT_Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NN = 1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(NN)
         IF (LPC%LIQUID_DROPLET .AND. LPC%Y_INDEX==Y_INDEX) &
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES +  AVG_DROP_TMP(II,JJ,KK,LAGRANGIAN_PARTICLE_CLASS(NN)%ARRAY_INDEX)-TMPM
      ENDDO

   CASE(182) ! ADA_Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NN = 1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(NN)
         IF (LPC%LIQUID_DROPLET .AND. LPC%Y_INDEX==Y_INDEX) &
            GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + AVG_DROP_AREA(II,JJ,KK,LAGRANGIAN_PARTICLE_CLASS(NN)%ARRAY_INDEX)
      ENDDO

   CASE(183) ! QABS_Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NN = 1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(NN)
         IF (LPC%LIQUID_DROPLET .AND. LPC%Y_INDEX==Y_INDEX) THEN
            IF (ABS(AVG_DROP_AREA(II,JJ,KK,LPC%ARRAY_INDEX))>TWENTY_EPSILON_EB) THEN
               DO N = 1,NUMBER_SPECTRAL_BANDS
                  IF (NUMBER_SPECTRAL_BANDS==1) THEN
                     BBF = 1._EB
                  ELSE
                     BBF = BLACKBODY_FRACTION(WL_LOW(N),WL_HIGH(N),RADTMP)
                  ENDIF
                  CALL INTERPOLATE1D(LPC%R50,LPC%WQABS(:,N),AVG_DROP_RAD(II,JJ,KK,LPC%ARRAY_INDEX),Q_SUM)
                  GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + BBF*Q_SUM
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES/REAL(NUMBER_SPECTRAL_BANDS,EB)

   CASE(184) ! QSCA_Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      DO NN = 1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(NN)
         IF (LPC%LIQUID_DROPLET .AND. LPC%Y_INDEX==Y_INDEX) THEN
            IF (ABS(AVG_DROP_AREA(II,JJ,KK,LPC%ARRAY_INDEX))>TWENTY_EPSILON_EB) THEN
               DO N = 1,NUMBER_SPECTRAL_BANDS
                  IF (NUMBER_SPECTRAL_BANDS==1) THEN
                     BBF = 1._EB
                  ELSE
                     BBF = BLACKBODY_FRACTION(WL_LOW(N),WL_HIGH(N),RADTMP)
                  ENDIF
                  CALL INTERPOLATE1D(LPC%R50,LPC%WQSCA(:,N),AVG_DROP_RAD(II,JJ,KK,LPC%ARRAY_INDEX),Q_SUM)
                  GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES + BBF*Q_SUM
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_RES/REAL(NUMBER_SPECTRAL_BANDS,EB)

   CASE(185) ! NUMBER OF PARTICLES
      GAS_PHASE_OUTPUT_RES = NLP

   CASE(186) ! DROPLET VOLUME FRACTION
      LPC => LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)
      GAS_PHASE_OUTPUT_RES = MIN(1._EB,AVG_DROP_DEN(II,JJ,KK,LPC%ARRAY_INDEX)/LPC%DENSITY)

   CASE(190) ! CELL PHASE
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID) GAS_PHASE_OUTPUT_RES=1._EB

   CASE(191) ! SCALAR UNKNOWN NUMBER
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (CC_IBM) GAS_PHASE_OUTPUT_RES = REAL(CCVAR(II,JJ,KK,CC_UNKZ),EB)

   CASE(192) ! F_X UNKNOWN NUMBER
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (CC_IBM) THEN
         GAS_PHASE_OUTPUT_RES = REAL(FCVAR(II,JJ,KK,CC_UNKF,IAXIS),EB)
         IF(FCVAR(II,JJ,KK,CC_IDRC,IAXIS)>0) GAS_PHASE_OUTPUT_RES = REAL(RC_FACE(FCVAR(II,JJ,KK,CC_IDRC,IAXIS))%UNKF,EB)
      ENDIF

   CASE(193) ! F_Y UNKNOWN NUMBER
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (CC_IBM) THEN
         GAS_PHASE_OUTPUT_RES = REAL(FCVAR(II,JJ,KK,CC_UNKF,JAXIS),EB)
         IF(FCVAR(II,JJ,KK,CC_IDRC,JAXIS)>0) GAS_PHASE_OUTPUT_RES = REAL(RC_FACE(FCVAR(II,JJ,KK,CC_IDRC,JAXIS))%UNKF,EB)
      ENDIF

   CASE(194) ! F_Z UNKNOWN NUMBER
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (CC_IBM) THEN
         GAS_PHASE_OUTPUT_RES = REAL(FCVAR(II,JJ,KK,CC_UNKF,KAXIS),EB)
         IF(FCVAR(II,JJ,KK,CC_IDRC,KAXIS)>0) GAS_PHASE_OUTPUT_RES = REAL(RC_FACE(FCVAR(II,JJ,KK,CC_IDRC,KAXIS))%UNKF,EB)
      ENDIF

   CASE(230) ! RANDOM NUMBER
      CALL RANDOM_NUMBER(RN)
      GAS_PHASE_OUTPUT_RES = REAL(RN,EB)

   CASE(231) ! PDPA
      GAS_PHASE_OUTPUT_RES = 0._EB

      PDPA_IF: IF ( (PY%PDPA_START<=T .AND. T<=PY%PDPA_END) .OR. .NOT.PY%PDPA_INTEGRATE ) THEN

         IF (.NOT.PY%PDPA_INTEGRATE) THEN
            DV%PDPA_NUMER = 0._EB
            DV%PDPA_DENOM = 0._EB
         ENDIF

         PDPA_FORMULA_SELECT: SELECT CASE(PY%QUANTITY)
            ! see user guide table: output quantities available for PDPA
            CASE DEFAULT;                  PDPA_FORMULA = 1
            CASE('ENTHALPY');              PDPA_FORMULA = 2
            CASE('PARTICLE FLUX X');       PDPA_FORMULA = 2
            CASE('PARTICLE FLUX Y');       PDPA_FORMULA = 2
            CASE('PARTICLE FLUX Z');       PDPA_FORMULA = 2
            CASE('U-VELOCITY');            PDPA_FORMULA = 1
            CASE('V-VELOCITY');            PDPA_FORMULA = 1
            CASE('W-VELOCITY');            PDPA_FORMULA = 1
            CASE('VELOCITY');              PDPA_FORMULA = 1
            CASE('TEMPERATURE');           PDPA_FORMULA = 1
            CASE('MASS CONCENTRATION');    PDPA_FORMULA = 2
            CASE('NUMBER CONCENTRATION');  PDPA_FORMULA = 2
         END SELECT PDPA_FORMULA_SELECT

         SELECT CASE(PDPA_FORMULA)
            CASE(1)
               IF (PY%PDPA_M-PY%PDPA_N==0) THEN
                  EXPON = 1._EB
               ELSE
                  EXPON = 1._EB/(PY%PDPA_M-PY%PDPA_N)
               ENDIF
            CASE(2)
               EXPON = 1._EB
               IF (PY%PDPA_NORMALIZE) THEN
                  DV%PDPA_DENOM = DV%PDPA_DENOM + FOTHPI*PY%PDPA_RADIUS**3
               ELSE
                  DV%PDPA_DENOM = 1._EB
               ENDIF
         END SELECT

         PDPA_PARTICLE_LOOP: DO I=1,NLP
            LP=>LAGRANGIAN_PARTICLE(I)
            LPC=>LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
            IF (PY%PART_INDEX/=LP%CLASS_INDEX .AND. PY%PART_INDEX/=-1) CYCLE PDPA_PARTICLE_LOOP
            BC => BOUNDARY_COORD(LP%BC_INDEX)
            IF ((BC%X-DV%X)**2+(BC%Y-DV%Y)**2+(BC%Z-DV%Z)**2 > PY%PDPA_RADIUS**2) CYCLE PDPA_PARTICLE_LOOP
            IF (.NOT.LPC%MASSLESS_TRACER) THEN
               B1 => BOUNDARY_PROP1(LP%B1_INDEX)
               R_D = LP%RADIUS
               TMP_F = B1%TMP_F
            ELSE
               R_D = 1._EB
               TMP_F = TMPA
            ENDIF
            ! see Table 20.1 in FDS User Guide
            PDPA_QUANTITY_SELECT: SELECT CASE(PY%QUANTITY)
               CASE DEFAULT;                  PHI = 1._EB
               CASE('ENTHALPY');              PHI = 0._EB
                  IF (LPC%SURF_INDEX==DROPLET_SURF_INDEX) THEN
                     CALL INTERPOLATE1D_UNIFORM(LBOUND(SPECIES(LPC%Y_INDEX)%C_P_L_BAR,1),&
                                                SPECIES(LPC%Y_INDEX)%C_P_L_BAR,TMP_F,CPBAR)
                     PHI = 0.001_EB*LPC%FTPR*R_D**3*CPBAR*TMP_F ! kJ
                  ELSEIF (LPC%SURF_INDEX>0) THEN
                     SF => SURFACE(LPC%SURF_INDEX)
                     IF (SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
                        ! SURFACE_DENSITY with MODE=3 returns energy density kJ/(m3-initial)
                        ! here VOL multiplies by the initial volume
                        SELECT CASE(SF%GEOMETRY)
                           CASE(SURF_CARTESIAN);   VOL = SF%LENGTH * SF%WIDTH * 2._EB*SF%THICKNESS
                           CASE(SURF_CYLINDRICAL); VOL = SF%LENGTH * PI*(SF%INNER_RADIUS+SF%THICKNESS)**2
                           CASE(SURF_SPHERICAL);   VOL = FOTHPI*(SF%INNER_RADIUS+SF%THICKNESS)**3
                        END SELECT
                        ONE_D => BOUNDARY_ONE_D(LP%OD_INDEX)
                        PHI = 0.001_EB*SURFACE_DENSITY(3,SF,ONE_D) * VOL ! kJ
                     ENDIF
                  ENDIF
               CASE('PARTICLE FLUX X');       PHI = LPC%FTPR*R_D**3*LP%U
               CASE('PARTICLE FLUX Y');       PHI = LPC%FTPR*R_D**3*LP%V
               CASE('PARTICLE FLUX Z');       PHI = LPC%FTPR*R_D**3*LP%W
               CASE('U-VELOCITY');            PHI = LP%U
               CASE('V-VELOCITY');            PHI = LP%V
               CASE('W-VELOCITY');            PHI = LP%W
               CASE('VELOCITY');              PHI = SQRT(LP%U**2 + LP%V**2 + LP%W**2)
               CASE('TEMPERATURE');           PHI = TMP_F - TMPM
               CASE('MASS CONCENTRATION');    PHI = LPC%FTPR*R_D**3
               CASE('NUMBER CONCENTRATION');  PHI = 1._EB
            END SELECT PDPA_QUANTITY_SELECT

            SELECT CASE(PDPA_FORMULA)
               CASE(1)
                  DV%PDPA_NUMER = DV%PDPA_NUMER + LP%PWT*(2._EB*R_D)**PY%PDPA_M * PHI
                  DV%PDPA_DENOM = DV%PDPA_DENOM + LP%PWT*(2._EB*R_D)**PY%PDPA_N
               CASE(2)
                  DV%PDPA_NUMER = DV%PDPA_NUMER + LP%PWT*PHI
            END SELECT

            IF (PY%HISTOGRAM)  CALL UPDATE_HISTOGRAM(PY%HISTOGRAM_NBINS,PY%HISTOGRAM_LIMITS,DV%HISTOGRAM_COUNTS,&
                                              (2._EB*R_D)**PY%PDPA_M * PHI,LP%PWT*R_D**PY%PDPA_N)

         ENDDO PDPA_PARTICLE_LOOP

         IF (DV%PDPA_DENOM>TWENTY_EPSILON_EB) GAS_PHASE_OUTPUT_RES = (DV%PDPA_NUMER/DV%PDPA_DENOM)**EXPON

      ENDIF PDPA_IF
   CASE(251)  ! WIND CHILL INDEX
      ! Wind speed at head height m/s, temperature Celsius
      ! WCT = 13.12 + 0.6215*TMP_G - 13.956*VEL_10m**(0.16) + 0.4867*TMP_G*VEL_10m**(0.16)
      ! Canada: Speed at head height = 2/3 * speed at 10 m height, v_10m = 1.5*v_head
      TMP_G = TMP(II,JJ,KK) - TMPM ! Temperature as Celsius
      VEL = 1.5_EB*SQRT(2._EB*KRES(II,JJ,KK)) ! Flow (wind) speed as m/s at 10 m height
      GAS_PHASE_OUTPUT_RES = MIN(13.12_EB+0.6215_EB*TMP_G-13.956_EB*VEL**(0.16_EB)+0.4867_EB*TMP_G*VEL**(0.16_EB),TMP_G)

   CASE(253)  ! ZONE PRESSURE SOLVER TYPE
      GAS_PHASE_OUTPUT_RES = REAL(PRES_FLAG,EB)
      IF (PRES_FLAG==ULMAT_FLAG) THEN
         IF (ZONE_MESH(ZONE_MESH(PRESSURE_ZONE(II,JJ,KK))%CONNECTED_ZONE_PARENT)%USE_FFT) THEN
            GAS_PHASE_OUTPUT_RES = REAL(FFT_FLAG,EB)
         ELSE
            ! uses PARDISO solver per mesh zone
            GAS_PHASE_OUTPUT_RES = REAL(ULMAT_FLAG,EB)
         ENDIF
      ENDIF

   CASE(254)  ! PRESSURE ZONE PARENT
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (PRES_FLAG==ULMAT_FLAG) THEN
         GAS_PHASE_OUTPUT_RES = REAL(ZONE_MESH(PRESSURE_ZONE(II,JJ,KK))%CONNECTED_ZONE_PARENT,EB)
      ENDIF

   CASE(500)  ! PRESSURE MMS
      XHAT = XC(II) - UF_MMS*T
      ZHAT = ZC(KK) - WF_MMS*T
      GAS_PHASE_OUTPUT_RES = VD2D_MMS_P_3(XHAT,ZHAT,T)
   CASE(501)  ! H MMS
      XHAT = XC(II) - UF_MMS*T
      ZHAT = ZC(KK) - WF_MMS*T
      GAS_PHASE_OUTPUT_RES = VD2D_MMS_H_3(XHAT,ZHAT,T)
   CASE(502)  ! CHI_R
      III = MAX(1,MIN(II,IBAR))
      JJJ = MAX(1,MIN(JJ,JBAR))
      KKK = MAX(1,MIN(KK,KBAR))
      GAS_PHASE_OUTPUT_RES = CHI_R(III,JJJ,KKK)
   CASE(504)  ! CFL 1
      IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID) THEN
         GAS_PHASE_OUTPUT_RES = 0._EB
      ELSE
         UODX = MAXVAL(ABS(US(II-1:II,JJ,KK)))*RDX(II)
         VODY = MAXVAL(ABS(VS(II,JJ-1:JJ,KK)))*RDY(JJ)
         WODZ = MAXVAL(ABS(WS(II,JJ,KK-1:KK)))*RDZ(KK)
         UVW = UODX + VODY + WODZ  + ABS(DS(II,JJ,KK)) ! CFL_VELOCITY_NORM=1
         GAS_PHASE_OUTPUT_RES = DT*UVW
      ENDIF
   CASE(505)  ! CFL 3
      IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID) THEN
         GAS_PHASE_OUTPUT_RES = 0._EB
      ELSE
         UODX = MAXVAL(ABS(US(II-1:II,JJ,KK)))*RDX(II)
         VODY = MAXVAL(ABS(VS(II,JJ-1:JJ,KK)))*RDY(JJ)
         WODZ = MAXVAL(ABS(WS(II,JJ,KK-1:KK)))*RDZ(KK)
         UVW = MAX(UODX,VODY,WODZ)                     ! CFL_VELOCITY_NORM=3
         GAS_PHASE_OUTPUT_RES = DT*UVW
      ENDIF

   CASE(508)  ! IDEAL GAS PRESSURE
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RCON)
      GAS_PHASE_OUTPUT_RES = RHO(II,JJ,KK)*RCON*TMP(II,JJ,KK)

   CASE(513)  ! DHDX
      GAS_PHASE_OUTPUT_RES = RDXN(II)*(HS(II+1,JJ,KK)-HS(II,JJ,KK))
   CASE(514)  ! DHDY
      GAS_PHASE_OUTPUT_RES = RDYN(JJ)*(HS(II,JJ+1,KK)-HS(II,JJ,KK))
   CASE(515)  ! DHDZ
      GAS_PHASE_OUTPUT_RES = RDZN(KK)*(HS(II,JJ,KK+1)-HS(II,JJ,KK))

   CASE(523)  ! ABSOLUTE PRESSURE
      GAS_PHASE_OUTPUT_RES  = PBAR(KK,PRESSURE_ZONE(II,JJ,KK)) + RHO(II,JJ,KK)*(H(II,JJ,KK)-KRES(II,JJ,KK))
   CASE(528)  ! ADVECTIVE MASS FLUX X
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = ADV_FX(II,JJ,KK,Z_INDEX)
      ELSEIF (Y_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES), ADV_FX(II,JJ,KK,1:N_TRACKED_SPECIES) )
      ENDIF
   CASE(529)  ! ADVECTIVE MASS FLUX Y
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = ADV_FY(II,JJ,KK,Z_INDEX)
      ELSEIF (Y_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES), ADV_FY(II,JJ,KK,1:N_TRACKED_SPECIES) )
      ENDIF
   CASE(530)  ! ADVECTIVE MASS FLUX Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = ADV_FZ(II,JJ,KK,Z_INDEX)
      ELSEIF (Y_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES), ADV_FZ(II,JJ,KK,1:N_TRACKED_SPECIES) )
      ENDIF
   CASE(531)  ! DIFFUSIVE MASS FLUX X
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DIF_FX(II,JJ,KK,Z_INDEX)
      ELSEIF (Y_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES), DIF_FX(II,JJ,KK,1:N_TRACKED_SPECIES) )
      ENDIF
   CASE(532)  ! DIFFUSIVE MASS FLUX Y
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DIF_FY(II,JJ,KK,Z_INDEX)
      ELSEIF (Y_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES), DIF_FY(II,JJ,KK,1:N_TRACKED_SPECIES) )
      ENDIF
   CASE(533)  ! DIFFUSIVE MASS FLUX Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DIF_FZ(II,JJ,KK,Z_INDEX)
      ELSEIF (Y_INDEX>0) THEN
         GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES), DIF_FZ(II,JJ,KK,1:N_TRACKED_SPECIES) )
      ENDIF
   CASE(534)  ! TOTAL MASS FLUX X
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX<=0 .AND. Y_INDEX<=0) THEN
         GAS_PHASE_OUTPUT_RES = SUM(ADV_FX(II,JJ,KK,:) + DIF_FX(II,JJ,KK,:))
      ELSE
         IF (Z_INDEX>0) THEN
            GAS_PHASE_OUTPUT_RES = ADV_FX(II,JJ,KK,Z_INDEX) + DIF_FX(II,JJ,KK,Z_INDEX)
         ELSEIF (Y_INDEX>0) THEN
            GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                                   (ADV_FX(II,JJ,KK,1:N_TRACKED_SPECIES) + DIF_FX(II,JJ,KK,1:N_TRACKED_SPECIES)) )
         ENDIF
      ENDIF
   CASE(535)  ! TOTAL MASS FLUX Y
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX<=0 .AND. Y_INDEX<=0) THEN
         GAS_PHASE_OUTPUT_RES = SUM(ADV_FY(II,JJ,KK,:) + DIF_FY(II,JJ,KK,:))
      ELSE
         IF (Z_INDEX>0) THEN
            GAS_PHASE_OUTPUT_RES = ADV_FY(II,JJ,KK,Z_INDEX) + DIF_FY(II,JJ,KK,Z_INDEX)
         ELSEIF (Y_INDEX>0) THEN
            GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                                   (ADV_FY(II,JJ,KK,1:N_TRACKED_SPECIES) + DIF_FY(II,JJ,KK,1:N_TRACKED_SPECIES)) )
         ENDIF
      ENDIF
   CASE(536)  ! TOTAL MASS FLUX Z
      GAS_PHASE_OUTPUT_RES = 0._EB
      IF (Z_INDEX<=0 .AND. Y_INDEX<=0) THEN
         GAS_PHASE_OUTPUT_RES = SUM(ADV_FZ(II,JJ,KK,:) + DIF_FZ(II,JJ,KK,:))
      ELSE
         IF (Z_INDEX>0) THEN
            GAS_PHASE_OUTPUT_RES = ADV_FZ(II,JJ,KK,Z_INDEX) + DIF_FZ(II,JJ,KK,Z_INDEX)
         ELSEIF (Y_INDEX>0) THEN
            GAS_PHASE_OUTPUT_RES = DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                                   (ADV_FZ(II,JJ,KK,1:N_TRACKED_SPECIES) + DIF_FZ(II,JJ,KK,1:N_TRACKED_SPECIES)) )
         ENDIF
      ENDIF

   CASE(550) ! CUTCELL VELOCITY DIVERGENCE
      GAS_PHASE_OUTPUT_RES = CARTVELDIV(II,JJ,KK)
   CASE(551) ! CARTESIAN VELOCITY DIVERGENCE
      GAS_PHASE_OUTPUT_RES = CARTVELDIV(II,JJ,KK)

   CASE(552) ! U_LS
      GAS_PHASE_OUTPUT_RES = U_LS(II,JJ)
   CASE(553) ! V_LS
      GAS_PHASE_OUTPUT_RES = V_LS(II,JJ)

 END SELECT IND_SELECT

! Fill GAS_PHASE_OUTPUT for CUT_CELLs.
! Some variables have already been filled in fire.f90
! Below we fill the values allocated in CC_CUTCELL_TYPE in type.f90

CC_IBM_IF: IF (CC_IBM) THEN

   IF (CCVAR(II,JJ,KK,CC_CGSC) == CC_SOLID .OR. CELL(CELL_INDEX(II,JJ,KK))%SOLID) EXIT CC_IBM_IF

   CCVAR_IF: IF (CCVAR(II,JJ,KK,CC_IDCC) > 0) THEN ! we have a cutcell
      ! cell centered quantities
      GAS_PHASE_OUTPUT_CC = 0._EB
      VC = 0._EB
      IF (PRESENT(ICC_IN)) THEN
         ICC = ICC_IN
         JCC_LO = JCC_IN
         JCC_HI = JCC_IN
      ELSE
         ICC=CCVAR(II,JJ,KK,CC_IDCC)
         NCELL=CUT_CELL(ICC)%NCELL
         JCC_LO = 1
         JCC_HI = NCELL
      ENDIF
      CC_LOOP: DO JCC=JCC_LO,JCC_HI
         ! Get species mass fraction if necessary
         Y_H2O     = 0._EB
         R_Y_H2O   = 0._EB
         Y_SPECIES = 1._EB
         IF (Z_INDEX > 0) THEN
            Y_SPECIES = CUT_CELL(ICC)%ZZ(Z_INDEX,JCC)
            RCON = SPECIES_MIXTURE(Z_INDEX)%RCON
         ELSEIF (Y_INDEX > 0) THEN
            ZZ_GET(1:N_TRACKED_SPECIES) = CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC)
            RCON = SPECIES(Y_INDEX)%RCON
            CALL GET_MASS_FRACTION(ZZ_GET,Y_INDEX,Y_SPECIES)
         ENDIF
         IF (DRY .AND. H2O_INDEX > 0) THEN
            ZZ_GET(1:N_TRACKED_SPECIES) = CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC)
            CALL GET_MASS_FRACTION(ZZ_GET,H2O_INDEX,Y_H2O)
            R_Y_H2O = SPECIES(H2O_INDEX)%RCON * Y_H2O
            IF (Y_INDEX==H2O_INDEX) Y_SPECIES=0._EB
         ENDIF
         VC = VC + CUT_CELL(ICC)%VOLUME(JCC)
         IND_SELECT_2: SELECT CASE(IND)
            CASE DEFAULT
               EXIT CCVAR_IF ! GAS_PHASE_OUTPUT_RES is unchanged
            CASE(1)   ! DENSITY
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + CUT_CELL(ICC)%RHO(JCC)*Y_SPECIES * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(5)   ! TEMPERATURE
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + (CUT_CELL(ICC)%TMP(JCC)-TMPM)    * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(9)   ! PRESSURE
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + &
                  ( PBAR(KK,PRESSURE_ZONE(II,JJ,KK)) + CUT_CELL(ICC)%RHO(JCC)*(CUT_CELL(ICC)%H(JCC)-KRES(II,JJ,KK))  - P_0(KK) ) &
                                                                                            * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(10)  ! VELOCITY
               IF(II<1 .OR. II>IBAR .OR. JJ<1 .OR. JJ>JBAR .OR. KK<1 .OR. KK>KBAR) THEN
                  GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + SQRT(2._EB*KRES(II,JJ,KK)) * CUT_CELL(ICC)%VOLUME(JCC)
               ELSE
                  CALL CC_CUTCELL_VELOCITY(NM,0._EB,ICC,JCC,VELOCITY_COMPONENT,ATOTV,RETURN_INTEGRALS=.FALSE.)
                  GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + SQRT(DOT_PRODUCT(VELOCITY_COMPONENT,VELOCITY_COMPONENT)) &
                                                                                            * CUT_CELL(ICC)%VOLUME(JCC)
               ENDIF
            CASE(11)  ! HRRPUV
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + CUT_CELL(ICC)%Q(JCC)*0.001       * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(12)  ! H
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + 0.5_EB*(CUT_CELL(ICC)%H(JCC)+CUT_CELL(ICC)%HS(JCC)) &
                                                                                            * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(14)  ! DIVERGENCE
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + CUT_CELL(ICC)%D(JCC)             * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(15)  ! MIXING TIME
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + CUT_CELL(ICC)%MIX_TIME(JCC)      * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(19)  ! RADIATION LOSS
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + CUT_CELL(ICC)%QR(JCC)*0.001      * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(22)  ! HS
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + CUT_CELL(ICC)%HS(JCC)            * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(90)  ! MASS FRACTION
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + Y_SPECIES/(1._EB-Y_H2O)          * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(94)  ! VOLUME FRACTION
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + &
                                     RCON*Y_SPECIES/CUT_CELL(ICC)%RSUM(JCC)/(1._EB-R_Y_H2O/CUT_CELL(ICC)%RSUM(JCC)) &
                                                                                            * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(138) ! HRRPUV REAC
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + CUT_CELL(ICC)%Q_REAC(JCC,REAC_INDEX)*0.001_EB &
                                                                                            * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(191) ! SCALAR UNKNOWN NUMBER
               GAS_PHASE_OUTPUT_RES = REAL(CUT_CELL(ICC)%UNKZ(JCC),EB); RETURN

            CASE(523) ! ABSOLUTE PRESSURE
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + &
                  ( PBAR(KK,PRESSURE_ZONE(II,JJ,KK)) + CUT_CELL(ICC)%RHO(JCC)*(CUT_CELL(ICC)%H(JCC)-KRES(II,JJ,KK)) ) &
                                                                                            * CUT_CELL(ICC)%VOLUME(JCC)
            CASE(550) ! CUTCELL VELOCITY DIVERGENCE
               GAS_PHASE_OUTPUT_CC = GAS_PHASE_OUTPUT_CC + CCVELDIV(II,JJ,KK)*CUT_CELL(ICC)%VOLUME(JCC)
         END SELECT IND_SELECT_2
      ENDDO CC_LOOP
      GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_CC/VC
   ENDIF CCVAR_IF

   AXIS = ABS(OUTPUT_QUANTITY(IND)%IOR)
   AXIS_IF: IF (AXIS>0) THEN
      FCVAR_IF: IF (FCVAR(II,JJ,KK,CC_IDCF,AXIS)>0) THEN
         ! face centered quantities
         GAS_PHASE_OUTPUT_CFA = 0._EB
         CFACE_AREA = 0._EB
         ICF=FCVAR(II,JJ,KK,CC_IDCF,AXIS)
         NFACE=CUT_FACE(ICF)%NFACE
         CFA_LOOP: DO JCF=1,NFACE
            CFACE_AREA = CFACE_AREA + CUT_FACE(ICF)%AREA(JCF)

            IND_SELECT_3: SELECT CASE(IND)
               CASE DEFAULT
                  EXIT FCVAR_IF ! GAS_PHASE_OUTPUT_RES is unchanged
               CASE( 2, 3, 4)  ! F_X, F_Y, F_Z
                  GAS_PHASE_OUTPUT_CFA = GAS_PHASE_OUTPUT_CFA + CUT_FACE(ICF)%FN(JCF) * CUT_FACE(ICF)%AREA(JCF)
               CASE(6)   ! U-VELOCITY
                  GAS_PHASE_OUTPUT_CFA = GAS_PHASE_OUTPUT_CFA + CUT_FACE(ICF)%VEL(JCF) * CUT_FACE(ICF)%AREA(JCF)
               CASE(7)   ! V-VELOCITY
                  GAS_PHASE_OUTPUT_CFA = GAS_PHASE_OUTPUT_CFA + CUT_FACE(ICF)%VEL(JCF) * CUT_FACE(ICF)%AREA(JCF)
               CASE(8)   ! W-VELOCITY
                  GAS_PHASE_OUTPUT_CFA = GAS_PHASE_OUTPUT_CFA + CUT_FACE(ICF)%VEL(JCF) * CUT_FACE(ICF)%AREA(JCF)
               CASE(192:194) ! F_X,F_Y,F_Z UNKNOWN NUMBER
                  GAS_PHASE_OUTPUT_RES = REAL(CUT_FACE(ICF)%UNKF(JCF),EB); RETURN
            END SELECT IND_SELECT_3

         ENDDO CFA_LOOP
         GAS_PHASE_OUTPUT_RES = GAS_PHASE_OUTPUT_CFA/CFACE_AREA
      ENDIF FCVAR_IF
   ENDIF AXIS_IF

ENDIF CC_IBM_IF

END FUNCTION GAS_PHASE_OUTPUT


!> \brief Compute solid phase device output quantities
!> \param INDX Output quantity index
!> \param Y_INDEX Index of primitive gas species
!> \param Z_INDEX Index of gas species mixture
!> \param PART_INDEX Index of Lagrangian particle class
!> \param OPT_WALL_INDEX Index of WALL boundary cell
!> \param OPT_LP_INDEX Index of Lagrangian particle
!> \param OPT_BNDF_INDEX Index of the boundary file
!> \param OPT_DEVC_INDEX Index of device
!> \param OPT_CFACE_INDEX Index of immersed boundary cell face
!> \param OPT_CUT_FACE_INDEX Index of the cut face
!> \param OPT_NODE_INDEX Index of internal heat conduction grid
!> \param OPT_PROF_INDEX Index of PROFile

REAL(EB) FUNCTION SOLID_PHASE_OUTPUT(INDX,T,NM,Y_INDEX,Z_INDEX,PART_INDEX,OPT_WALL_INDEX,OPT_LP_INDEX,OPT_BNDF_INDEX,&
                                     OPT_DEVC_INDEX,OPT_CFACE_INDEX,OPT_CUT_FACE_INDEX,OPT_NODE_INDEX,OPT_PROF_INDEX)

USE PHYSICAL_FUNCTIONS, ONLY: SURFACE_DENSITY,GET_MASS_FRACTION,GET_SENSIBLE_ENTHALPY,GET_SPECIFIC_HEAT,GET_CONDUCTIVITY,&
                              GET_VISCOSITY,HEAT_TRANSFER_COEFFICIENT
USE TURBULENCE, ONLY: TAU_WALL_IJ
INTEGER, INTENT(IN), OPTIONAL :: OPT_WALL_INDEX,OPT_LP_INDEX,OPT_CFACE_INDEX,OPT_BNDF_INDEX,OPT_DEVC_INDEX,OPT_CUT_FACE_INDEX,&
                                 OPT_NODE_INDEX,OPT_PROF_INDEX
INTEGER, INTENT(IN) :: INDX,Y_INDEX,Z_INDEX,PART_INDEX,NM
REAL(EB),INTENT(IN) :: T
REAL(EB) :: Q_CON,RHOSUM,VOLSUM,ZZ_GET(1:N_TRACKED_SPECIES),Y_SPECIES,DEPTH,ASH_DEPTH,UN,H_S,RHO_D_DYDN,U_CELL,V_CELL,W_CELL,&
            LTMP,ATMP,CTMP,H_W_EFF,X0,VOL,DN,PRESS,&
            NVEC(3),PVEC(3),TAU_IJ(3,3),VEL_CELL(3),VEL_WALL(3),MU_WALL,RHO_WALL,FVEC(3),SVEC(3),TVEC1(3),TVEC2(3),&
            PR1,PR2,Z1,Z2,RADIUS,CUT_FACE_AREA,SOLID_PHASE_OUTPUT_CTF,AAA,BBB,CCC,ALP,BET,GAM,MMM,DTMP,HTC
INTEGER :: I_DEPTH,II2,IIG,JJG,KKG,NN,IWX,SURF_INDEX,I,J,NWP,M_INDEX,ICC,IND1,IND2,IC2,ITMP,ICF,JCF,NFACE,NR,MATL_INDEX,OUTPUT_INDEX
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_PROP2_TYPE), POINTER :: B2
TYPE(BOUNDARY_RADIA_TYPE), POINTER :: BR
TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

! Assign default value to output

SOLID_PHASE_OUTPUT = OUTPUT_QUANTITY(-INDX)%AMBIENT_VALUE

IF (PRESENT(OPT_DEVC_INDEX)) DV => DEVICE(OPT_DEVC_INDEX)

IF (PRESENT(OPT_WALL_INDEX)) THEN

   IF (OPT_WALL_INDEX==0) RETURN
   IWX = OPT_WALL_INDEX
   WC=>WALL(IWX)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) RETURN
   SURF_INDEX = WC%SURF_INDEX
   IF (WC%OD_INDEX>0) ONE_D => BOUNDARY_ONE_D(WC%OD_INDEX)
   IF (WC%BC_INDEX>0) BC    => BOUNDARY_COORD(WC%BC_INDEX)
   IF (WC%B1_INDEX>0) B1    => BOUNDARY_PROP1(WC%B1_INDEX)
   IF (WC%B2_INDEX>0) B2    => BOUNDARY_PROP2(WC%B2_INDEX)
   IF (WC%BR_INDEX>0) BR    => BOUNDARY_RADIA(WC%BR_INDEX)

ELSEIF (PRESENT(OPT_LP_INDEX)) THEN

   LP => LAGRANGIAN_PARTICLE(OPT_LP_INDEX)
   IF (LP%OD_INDEX>0) ONE_D => BOUNDARY_ONE_D(LP%OD_INDEX)
   IF (LP%BC_INDEX>0) BC    => BOUNDARY_COORD(LP%BC_INDEX)
   IF (LP%B1_INDEX>0) B1    => BOUNDARY_PROP1(LP%B1_INDEX)
   IF (LP%B2_INDEX>0) B2    => BOUNDARY_PROP2(LP%B2_INDEX)
   IF (LP%BR_INDEX>0) BR    => BOUNDARY_RADIA(LP%BR_INDEX)
   SURF_INDEX = LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)%SURF_INDEX

ELSEIF (PRESENT(OPT_CFACE_INDEX)) THEN

   CFA => CFACE(OPT_CFACE_INDEX)
   SURF_INDEX = CFA%SURF_INDEX
   IF (CFA%OD_INDEX>0) ONE_D => BOUNDARY_ONE_D(CFA%OD_INDEX)
   IF (CFA%BC_INDEX>0) BC    => BOUNDARY_COORD(CFA%BC_INDEX)
   IF (CFA%B1_INDEX>0) B1    => BOUNDARY_PROP1(CFA%B1_INDEX)
   IF (CFA%B2_INDEX>0) B2    => BOUNDARY_PROP2(CFA%B2_INDEX)
   IF (CFA%BR_INDEX>0) BR    => BOUNDARY_RADIA(CFA%BR_INDEX)

ENDIF

MATL_INDEX = 0
IF (PRESENT(OPT_DEVC_INDEX)) MATL_INDEX = DEVICE(OPT_DEVC_INDEX)%MATL_INDEX
IF (PRESENT(OPT_BNDF_INDEX)) MATL_INDEX = BOUNDARY_FILE(OPT_BNDF_INDEX)%MATL_INDEX
IF (PRESENT(OPT_PROF_INDEX)) MATL_INDEX = PROFILE(OPT_PROF_INDEX)%MATL_INDEX

ICF = 0
IF (PRESENT(OPT_CUT_FACE_INDEX)) ICF = OPT_CUT_FACE_INDEX

SF => SURFACE(SURF_INDEX)

! Special cases where an in-depth quantity is needed

IF (OUTPUT_QUANTITY(-INDX)%INSIDE_SOLID) THEN
   IF (SF%THERMAL_BC_INDEX/=THERMALLY_THICK) RETURN
   IF (PRESENT(OPT_NODE_INDEX)) THEN
      I_DEPTH = OPT_NODE_INDEX
   ELSE
      I_DEPTH = DV%I_DEPTH
      IF (ONE_D%PYROLYSIS_MODEL==PYROLYSIS_PREDICTED .OR. I_DEPTH==-1) THEN
         IF (DV%DEPTH > TWENTY_EPSILON_EB) THEN
            DEPTH = DV%DEPTH
         ELSE
            DEPTH = MAX(0._EB,SUM(ONE_D%LAYER_THICKNESS)+DV%DEPTH)
         ENDIF
         II2 = SUM(ONE_D%N_LAYER_CELLS)
         IF (DEPTH>SUM(ONE_D%LAYER_THICKNESS)) THEN
            RETURN
         ELSE
            DO II2=II2,1,-1
               IF (DEPTH<=ONE_D%X(II2)) I_DEPTH = II2
            ENDDO
         ENDIF
      ENDIF
   ENDIF
ENDIF

! Find the appropriate solid phase output quantity

SOLID_PHASE_SELECT: SELECT CASE(INDX)
   CASE( 1) ! RADIATIVE HEAT FLUX
      SOLID_PHASE_OUTPUT = (B1%Q_RAD_IN-B1%Q_RAD_OUT)*0.001_EB
   CASE( 2) ! CONVECTIVE HEAT FLUX
      SOLID_PHASE_OUTPUT = B1%Q_CON_F*0.001_EB
   CASE( 3) ! NORMAL VELOCITY
      SELECT CASE(BC%IOR)
         CASE( 1) ; SOLID_PHASE_OUTPUT = -U(BC%IIG-1,BC%JJG,BC%KKG)
         CASE(-1) ; SOLID_PHASE_OUTPUT =  U(BC%IIG  ,BC%JJG,BC%KKG)
         CASE( 2) ; SOLID_PHASE_OUTPUT = -V(BC%IIG,BC%JJG-1,BC%KKG)
         CASE(-2) ; SOLID_PHASE_OUTPUT =  V(BC%IIG,BC%JJG  ,BC%KKG)
         CASE( 3) ; SOLID_PHASE_OUTPUT = -W(BC%IIG,BC%JJG,BC%KKG-1)
         CASE(-3) ; SOLID_PHASE_OUTPUT =  W(BC%IIG,BC%JJG,BC%KKG  )
      END SELECT
      IF(PRESENT(OPT_CFACE_INDEX)) THEN
         IND1 = CFA%CUT_FACE_IND1
         IND2 = CFA%CUT_FACE_IND2
         SOLID_PHASE_OUTPUT = CUT_FACE(IND1)%VEL(IND2)
      ELSEIF (ICF>0) THEN
         SOLID_PHASE_OUTPUT_CTF = 0._EB
         CUT_FACE_AREA = 0._EB
         NFACE=CUT_FACE(ICF)%NFACE
         DO JCF=1,NFACE
            CUT_FACE_AREA = CUT_FACE_AREA + CUT_FACE(ICF)%AREA(JCF)
            SOLID_PHASE_OUTPUT_CTF = SOLID_PHASE_OUTPUT_CTF &
                                   - SIGN(1._EB,REAL(BC%IOR,EB))*CUT_FACE(ICF)%VEL(JCF)*CUT_FACE(ICF)%AREA(JCF)
         ENDDO
         SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT_CTF/CUT_FACE_AREA
      ENDIF
   CASE( 4) ! GAS TEMPERATURE
      SOLID_PHASE_OUTPUT = B1%TMP_G - TMPM
   CASE( 5) ! WALL TEMPERATURE
      SOLID_PHASE_OUTPUT = B1%TMP_F - TMPM
   CASE( 6) ! INSIDE WALL TEMPERATURE
      SOLID_PHASE_OUTPUT = ONE_D%TMP(I_DEPTH) - TMPM
   CASE( 7) ! BURNING RATE
      IF (N_REACTIONS>0) THEN
         SOLID_PHASE_OUTPUT = 0._EB
         DO NR=1,N_REACTIONS
            SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT+B1%M_DOT_G_PP_ACTUAL(REACTION(NR)%FUEL_SMIX_INDEX)*B1%AREA_ADJUST
         ENDDO
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF
   CASE( 8) ! NORMALIZED HEAT RELEASE RATE
      SOLID_PHASE_OUTPUT = 0._EB
      DO NR=1,N_REACTIONS
         SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT+B1%M_DOT_G_PP_ADJUST(REACTION(NR)%FUEL_SMIX_INDEX)*&
                              REACTION(NR)%HOC_COMPLETE
      ENDDO
      SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT*0.001_EB/(SF%SURFACE_DENSITY*B1%AREA_ADJUST)
   CASE( 9) ! HRRPUA
      SOLID_PHASE_OUTPUT = 0._EB
      DO NR=1,N_REACTIONS
         SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT+B1%M_DOT_G_PP_ADJUST(REACTION(NR)%FUEL_SMIX_INDEX)*&
                              REACTION(NR)%HOC_COMPLETE
      ENDDO
      SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT*0.001_EB
   CASE(10) ! TOTAL HEAT FLUX
      SOLID_PHASE_OUTPUT = (B1%Q_RAD_IN-B1%Q_RAD_OUT+B1%Q_CON_F)*0.001_EB
   CASE(11) ! PRESSURE COEFFICIENT
      IF (PRESENT(OPT_WALL_INDEX)) THEN
         IIG = BC%IIG
         JJG = BC%JJG
         KKG = BC%KKG
         PRESS = RHO(IIG,JJG,KKG)*(H(IIG,JJG,KKG)-KRES(IIG,JJG,KKG))
      ELSEIF (PRESENT(OPT_CFACE_INDEX)) THEN
         IND1 = CFA%CUT_FACE_IND1
         IND2 = CFA%CUT_FACE_IND2
         CALL GET_PRES_CFACE(PRESS,IND1,IND2,CFA)
      ELSE
         PRESS = 0._EB
      ENDIF
      SOLID_PHASE_OUTPUT = PRESS/(0.5_EB*RHOA*PY%CHARACTERISTIC_VELOCITY**2)
   CASE(12) ! BACK WALL TEMPERATURE
      SOLID_PHASE_OUTPUT = B1%TMP_B - TMPM
   CASE(13) ! GAUGE HEAT FLUX
      IF (PY%HEAT_TRANSFER_COEFFICIENT>=0._EB) THEN
         Q_CON = PY%HEAT_TRANSFER_COEFFICIENT*(TMP(BC%IIG,BC%JJG,BC%KKG)-PY%GAUGE_TEMPERATURE)
      ELSE
         IF (SF%BLOWING) THEN
            IF (OPT_WALL_INDEX > 0) THEN
               HTC=HEAT_TRANSFER_COEFFICIENT(NM,T,TMP(BC%IIG,BC%JJG,BC%KKG)-B1%TMP_F,SF,WALL_INDEX_IN=OPT_WALL_INDEX,&
                                             SKIP_BLOWING=.TRUE.)
            ELSEIF (OPT_CFACE_INDEX > 0) THEN
               HTC=HEAT_TRANSFER_COEFFICIENT(NM,T,TMP(BC%IIG,BC%JJG,BC%KKG)-B1%TMP_F,SF,CFACE_INDEX_IN=OPT_CFACE_INDEX,&
                                             SKIP_BLOWING=.TRUE.)
            ELSE
               HTC=HEAT_TRANSFER_COEFFICIENT(NM,T,TMP(BC%IIG,BC%JJG,BC%KKG)-B1%TMP_F,SF,PARTICLE_INDEX_IN=OPT_LP_INDEX,&
                                             SKIP_BLOWING=.TRUE.)
            ENDIF
            Q_CON = HTC*(TMP(BC%IIG,BC%JJG,BC%KKG)-PY%GAUGE_TEMPERATURE)
         ELSE
            Q_CON = B1%HEAT_TRANS_COEF*(TMP(BC%IIG,BC%JJG,BC%KKG)-PY%GAUGE_TEMPERATURE)
         ENDIF
      ENDIF
      SOLID_PHASE_OUTPUT = (PY%GAUGE_EMISSIVITY*(B1%Q_RAD_IN/(B1%EMISSIVITY+1.0E-10_EB) - SIGMA*PY%GAUGE_TEMPERATURE**4) + &
                            Q_CON)*0.001_EB
   CASE(14) ! NORMALIZED HEATING RATE
      SOLID_PHASE_OUTPUT = B1%Q_CON_F*0.001_EB/SF%SURFACE_DENSITY
   CASE(15,16) ! MASS FLUX, NORMALIZED MASS LOSS RATE
      IF (Z_INDEX >=0) THEN
         SOLID_PHASE_OUTPUT = B1%M_DOT_G_PP_ACTUAL(Z_INDEX)*B1%AREA_ADJUST
      ELSEIF (Y_INDEX > 0) THEN
         SOLID_PHASE_OUTPUT = &
            DOT_PRODUCT(Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),B1%M_DOT_G_PP_ACTUAL(1:N_TRACKED_SPECIES))*B1%AREA_ADJUST
      ELSEIF (MATL_INDEX>0) THEN
         SOLID_PHASE_OUTPUT = 0._EB
         M_INDEX = 0
         DO NN=1,ONE_D%N_MATL
            IF (MATL_INDEX==ONE_D%MATL_INDEX(NN)) THEN
               M_INDEX = NN
               EXIT
            ENDIF
         ENDDO
         IF (M_INDEX>0) SOLID_PHASE_OUTPUT = ONE_D%M_DOT_S_PP(M_INDEX)*B1%AREA_ADJUST
      ELSE
         SOLID_PHASE_OUTPUT = SUM(B1%M_DOT_G_PP_ACTUAL(:))*B1%AREA_ADJUST
      ENDIF
      IF (INDX==16) SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT/(SF%SURFACE_DENSITY*B1%AREA_ADJUST)
   CASE(17) ! RADIANCE
      IF (ASSOCIATED(BR)) THEN
         SOLID_PHASE_OUTPUT = SUM(BR%IL(1:NUMBER_SPECTRAL_BANDS))*0.001_EB
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF
   CASE(20) ! INCIDENT HEAT FLUX
      SOLID_PHASE_OUTPUT = ( B1%Q_RAD_IN/(B1%EMISSIVITY+1.0E-10_EB) )*0.001_EB
   CASE(21) ! HEAT TRANSFER COEFFICENT
      SOLID_PHASE_OUTPUT = B1%HEAT_TRANS_COEF
   CASE(22) ! RADIOMETER
      SOLID_PHASE_OUTPUT = PY%GAUGE_EMISSIVITY*(B1%Q_RAD_IN/(B1%EMISSIVITY+1.0E-10_EB)-SIGMA*PY%GAUGE_TEMPERATURE**4)*0.001_EB

   CASE(23) ! ADIABATIC SURFACE TEMPERATURE (Ferrari's Method for Solving the Quartic)
      H_W_EFF = 0._EB
      LTMP = 0._EB
      IF ((PRESENT(OPT_WALL_INDEX).OR.PRESENT(OPT_CFACE_INDEX)) .AND. ASSOCIATED(B2)) THEN  ! Look for evaporating liquid
         IF (ANY(ABS(B2%LP_CPUA)>TWENTY_EPSILON_EB)) THEN
            ATMP = 0._EB
            CTMP = 0._EB
            DO NN = 1,N_LAGRANGIAN_CLASSES
               LPC => LAGRANGIAN_PARTICLE_CLASS(NN)
               IF (LPC%LIQUID_DROPLET) THEN
                  CTMP = CTMP + B2%LP_CPUA(LPC%ARRAY_INDEX)
                  ATMP = ATMP + ABS(B2%LP_CPUA(LPC%ARRAY_INDEX))
                  LTMP = LTMP + ABS(B2%LP_CPUA(LPC%ARRAY_INDEX))*B2%LP_TEMP(LPC%ARRAY_INDEX)
               ENDIF
            ENDDO
            LTMP = LTMP / (ATMP + TWENTY_EPSILON_EB)
            H_W_EFF = CTMP / (B1%TMP_F-LTMP+TWENTY_EPSILON_EB)
            H_W_EFF = MIN(10000._EB,MAX(0._EB,H_W_EFF));
         ENDIF
      ENDIF
      AAA = B1%EMISSIVITY*SIGMA
      IF (B1%HEAT_TRANS_COEF+H_W_EFF>1.E-5_EB .AND. ABS(B1%Q_RAD_IN-AAA*TMP(BC%IIG,BC%JJG,BC%KKG)**4)>5.E-3_EB) THEN
         AAA = B1%EMISSIVITY*SIGMA
         BBB = B1%HEAT_TRANS_COEF + H_W_EFF
         CCC = -B1%Q_RAD_IN - B1%HEAT_TRANS_COEF*TMP(BC%IIG,BC%JJG,BC%KKG) - H_W_EFF*LTMP
         ALP = (SR3*SQRT(MAX(0._EB,27._EB*AAA**2*BBB**4-256._EB*AAA**3*CCC**3))+9._EB*AAA*BBB**2)**ONTH
         BET = FTTOT*CCC
         GAM = EIONTH*AAA
         MMM = SQRT(MAX(0._EB,BET/ALP + ALP/GAM))
         SOLID_PHASE_OUTPUT = 0.5_EB*(-MMM+SQRT(MAX(0._EB,2._EB*BBB/(AAA*MMM+TWENTY_EPSILON_EB)-MMM**2))) - TMPM
      ELSE
         SOLID_PHASE_OUTPUT = (B1%Q_RAD_IN/(B1%EMISSIVITY*SIGMA+TWENTY_EPSILON_EB))**0.25 - TMPM
      ENDIF
   CASE(24) ! WALL THICKNESS
      IF (SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
         SOLID_PHASE_OUTPUT = SUM(ONE_D%LAYER_THICKNESS)
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF

   CASE(25,26) ! SURFACE DENSITY, NORMALIZED MASS
      IF (SF%THERMAL_BC_INDEX/=THERMALLY_THICK) RETURN
      M_INDEX = 0
      IF (MATL_INDEX>0 .AND. ALLOCATED(ONE_D%MATL_INDEX)) THEN
         DO NN=1,ONE_D%N_MATL
            IF (MATL_INDEX==ONE_D%MATL_INDEX(NN)) THEN
               M_INDEX = NN
               EXIT
            ENDIF
         ENDDO
         IF (M_INDEX==0) THEN  ! There is none of the specified MATL within the surface
            SOLID_PHASE_OUTPUT = 0._EB
            RETURN
         ENDIF
      ENDIF
      IF (M_INDEX>0) THEN
         IF (PRESENT(OPT_LP_INDEX)) THEN
            SOLID_PHASE_OUTPUT = SURFACE_DENSITY(0,SF,ONE_D,MATL_INDEX=M_INDEX)
         ELSEIF (PRESENT(OPT_WALL_INDEX)) THEN
            SOLID_PHASE_OUTPUT = SURFACE_DENSITY(0,SF,ONE_D,MATL_INDEX=M_INDEX)
         ENDIF
      ELSE
         IF (PRESENT(OPT_LP_INDEX)) THEN
            SOLID_PHASE_OUTPUT = SURFACE_DENSITY(0,SF,ONE_D)
         ELSEIF (PRESENT(OPT_WALL_INDEX)) THEN
            SOLID_PHASE_OUTPUT = SURFACE_DENSITY(0,SF,ONE_D)
         ENDIF
      ENDIF
      IF (INDX==25 .AND. SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
         RADIUS = SF%INNER_RADIUS + ONE_D%X(SUM(ONE_D%N_LAYER_CELLS)) - ONE_D%X(0)
         IF (RADIUS>TWENTY_EPSILON_EB) THEN
            SELECT CASE(SF%GEOMETRY)
               CASE(SURF_CYLINDRICAL,SURF_INNER_CYLINDRICAL) ; SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT* SF%THICKNESS/RADIUS
               CASE(SURF_SPHERICAL)                          ; SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT*(SF%THICKNESS/RADIUS)**2
            END SELECT
         ELSE
            SOLID_PHASE_OUTPUT = 0._EB
         ENDIF
      ENDIF
      SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT*B1%AREA_ADJUST
      IF (INDX==26) SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT/SF%SURFACE_DENSITY

   CASE(27) ! SOLID DENSITY
      SOLID_PHASE_OUTPUT = 0._EB
      IF (MATL_INDEX > 0) THEN
         DO NN=1,ONE_D%N_MATL
            IF (MATL_INDEX==ONE_D%MATL_INDEX(NN)) THEN
               SOLID_PHASE_OUTPUT = ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)
               RETURN
            ENDIF
         ENDDO
      ELSE
         DO NN=1,ONE_D%N_MATL
            SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT + ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)
         ENDDO
      ENDIF

   CASE(28) ! EMISSIVITY
      SOLID_PHASE_OUTPUT = B1%EMISSIVITY

   CASE(29) ! SURFACE DEPOSITION
         IF (Z_INDEX>0) SOLID_PHASE_OUTPUT = B1%AWM_AEROSOL(SPECIES_MIXTURE(Z_INDEX)%AWM_INDEX)
         IF (Y_INDEX>0) SOLID_PHASE_OUTPUT = B1%AWM_AEROSOL(SPECIES(Y_INDEX)%AWM_INDEX)

   CASE(30:32) ! MPUA, CPUA, AMPUA
      LPC => LAGRANGIAN_PARTICLE_CLASS(PART_INDEX)
      IF (.NOT.PRESENT(OPT_LP_INDEX) .AND. ASSOCIATED(B2)) THEN
         SELECT CASE(INDX)
            CASE(30) ; SOLID_PHASE_OUTPUT = B2%LP_MPUA(LPC%ARRAY_INDEX)
            CASE(31) ; SOLID_PHASE_OUTPUT = B2%LP_CPUA(LPC%ARRAY_INDEX)*0.001_EB
            CASE(32) ; SOLID_PHASE_OUTPUT = B2%A_LP_MPUA(LPC%ARRAY_INDEX)
         END SELECT
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF

   CASE(33) ! SOLID SPECIFIC HEAT
      SOLID_PHASE_OUTPUT = 0._EB
      RHOSUM = 0._EB
      MATERIAL_LOOP_CP: DO NN=1,ONE_D%N_MATL
         IF (ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)<=TWENTY_EPSILON_EB) CYCLE MATERIAL_LOOP_CP
         RHOSUM = RHOSUM + ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)
         ML  => MATERIAL(ONE_D%MATL_INDEX(NN))
         ITMP = MIN(I_MAX_TEMP,NINT(ONE_D%TMP(I_DEPTH)))
         SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT + ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)*ML%C_S(ITMP)
      ENDDO MATERIAL_LOOP_CP
      SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT / RHOSUM * 0.001_EB

   CASE(34) ! SOLID CONDUCTIVITY
      SOLID_PHASE_OUTPUT = 0._EB
      VOLSUM = 0._EB
      MATERIAL_LOOP_K: DO NN=1,ONE_D%N_MATL
         IF (ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)<=TWENTY_EPSILON_EB) CYCLE MATERIAL_LOOP_K
         ML => MATERIAL(ONE_D%MATL_INDEX(NN))
         VOLSUM = VOLSUM + ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)/ML%RHO_S
         ITMP = MIN(I_MAX_TEMP,NINT(ONE_D%TMP(I_DEPTH)))
         SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT + ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)*ML%K_S(ITMP)/ML%RHO_S
      ENDDO MATERIAL_LOOP_K
      SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT / VOLSUM

   CASE(35) ! VISCOUS WALL UNITS (distance from the wall expressed in nondimensional viscous units)
      IF ((PRESENT(OPT_WALL_INDEX).OR.PRESENT(OPT_CFACE_INDEX)) .AND. ASSOCIATED(B2)) THEN
         SOLID_PHASE_OUTPUT = B2%Y_PLUS
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF

   CASE(36) ! FRICTION VELOCITY
      IF ((PRESENT(OPT_WALL_INDEX).OR.PRESENT(OPT_CFACE_INDEX)) .AND. ASSOCIATED(B2)) THEN
         SOLID_PHASE_OUTPUT = B2%U_TAU
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF

   CASE(37) ! VELOCITY ERROR
      SOLID_PHASE_OUTPUT = B1%VEL_ERR_NEW

   CASE(38) ! WALL VISCOSITY
      IF (PRESENT(OPT_WALL_INDEX)) THEN
         SOLID_PHASE_OUTPUT = MU(BC%IIG,BC%JJG,BC%KKG)
      ELSEIF (PRESENT(OPT_CFACE_INDEX)) THEN
         SOLID_PHASE_OUTPUT = CFA%MU_G
      ENDIF

   CASE(39) ! DEPOSITION VELOCITY
      SOLID_PHASE_OUTPUT = B2%V_DEP

   CASE(41) ! WALL CELL COLOR (output VENT index for WC color)
      SOLID_PHASE_OUTPUT = REAL(WC%VENT_INDEX,EB)

   CASE(42:44) ! MPUA_Z, CPUA_Z, AMPUA_Z
      SOLID_PHASE_OUTPUT = 0._EB
      IF (ASSOCIATED(B2)) THEN
         DO NN = 1,N_LAGRANGIAN_CLASSES
            LPC => LAGRANGIAN_PARTICLE_CLASS(NN)
            IF (LPC%LIQUID_DROPLET .AND. LPC%Y_INDEX==Y_INDEX) THEN
               SELECT CASE(INDX)
                  CASE(42) ; SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT + B2%LP_MPUA(LPC%ARRAY_INDEX)
                  CASE(43) ; SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT + B2%LP_CPUA(LPC%ARRAY_INDEX)*0.001_EB
                  CASE(44) ; SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT + B2%A_LP_MPUA(LPC%ARRAY_INDEX)
               END SELECT
            ENDIF
         ENDDO
      ENDIF

   CASE(45) ! WALL CELL BOUNDARY TYPE (debug)
      SOLID_PHASE_OUTPUT = REAL(WC%BOUNDARY_TYPE,EB)

   CASE(46) ! WALL CELL THERMAL BOUNDARY TYPE (debug)
      SOLID_PHASE_OUTPUT = REAL(SF%THERMAL_BC_INDEX,EB)

   CASE(47) ! INSIDE WALL DEPTH (for use with INSIDE WALL TEMPERATURE to obtain exact TMP location)
      IF (DV%DEPTH>SUM(ONE_D%LAYER_THICKNESS)) THEN
         SOLID_PHASE_OUTPUT = DV%DEPTH
      ELSE
         SOLID_PHASE_OUTPUT = 0.5_EB*( ONE_D%X(I_DEPTH-1) + ONE_D%X(I_DEPTH) )
      ENDIF

   CASE(48) ! LAYER DIVIDE DEPTH
      IF (SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
         SOLID_PHASE_OUTPUT = ONE_D%LAYER_DIVIDE_DEPTH
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF

   CASE(51)  ! ENTHALPY FLUX WALL
      ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
      CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,B1%TMP_F)
      SOLID_PHASE_OUTPUT = (-B1%RHO_F*H_S*B1%U_NORMAL &
                            -2._EB*B1%K_G*(TMP(BC%IIG,BC%JJG,BC%KKG)-B1%TMP_F)*B1%RDN)*0.001_EB

   CASE(60)  ! MASS FLUX WALL
      IF (PRESENT(OPT_WALL_INDEX)) THEN
         SELECT CASE(BC%IOR)
            CASE(-1) ; UN = -U(BC%IIG  ,BC%JJG,BC%KKG)
            CASE( 1) ; UN =  U(BC%IIG-1,BC%JJG,BC%KKG)
            CASE(-2) ; UN = -V(BC%IIG,BC%JJG  ,BC%KKG)
            CASE( 2) ; UN =  V(BC%IIG,BC%JJG-1,BC%KKG)
            CASE(-3) ; UN = -W(BC%IIG,BC%JJG,BC%KKG  )
            CASE( 3) ; UN =  W(BC%IIG,BC%JJG,BC%KKG-1)
         END SELECT
      ELSE
         UN = -B1%U_NORMAL
      ENDIF
      IF (Z_INDEX > 0) THEN
         Y_SPECIES = B1%ZZ_F(Z_INDEX)
         RHO_D_DYDN = B1%RHO_D_DZDN_F(Z_INDEX)
      ELSEIF (Y_INDEX > 0) THEN
         ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
         CALL GET_MASS_FRACTION(ZZ_GET,Y_INDEX,Y_SPECIES)
         RHO_D_DYDN = DOT_PRODUCT(Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),B1%RHO_D_DZDN_F(1:N_TRACKED_SPECIES))
      ELSE
         Y_SPECIES = 1._EB
         RHO_D_DYDN = 0._EB
      ENDIF
      ! convention here is: inflow is positive (adds mass to domain), outflow is negative (subtracts mass)
      SOLID_PHASE_OUTPUT = Y_SPECIES*B1%RHO_F*UN - RHO_D_DYDN

   CASE(61) ! GAS DENSITY
      IF (Z_INDEX > 0) THEN
         Y_SPECIES = B1%ZZ_F(Z_INDEX)
      ELSEIF (Y_INDEX > 0) THEN
         ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
         CALL GET_MASS_FRACTION(ZZ_GET,Y_INDEX,Y_SPECIES)
      ELSE
         Y_SPECIES = 1._EB
      ENDIF
      SOLID_PHASE_OUTPUT = B1%RHO_G*Y_SPECIES

   CASE(63) ! THERMAL WALL UNITS
      IF ((PRESENT(OPT_WALL_INDEX).OR.PRESENT(OPT_CFACE_INDEX)) .AND. ASSOCIATED(B2)) THEN
         SOLID_PHASE_OUTPUT = B2%Z_STAR
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF

   CASE(64) ! TOTAL MASS FLUX WALL
      IIG = BC%IIG
      JJG = BC%JJG
      KKG = BC%KKG
      IF (STORE_SPECIES_FLUX .AND. PRESENT(OPT_WALL_INDEX)) THEN ! Case of External walls or Obsts.
         ! convention here is: inflow is positive (adds mass to domain), outflow is negative (subtracts mass)
         IF (Z_INDEX>0) THEN
            SELECT CASE(BC%IOR)
            CASE(-1); SOLID_PHASE_OUTPUT=-(ADV_FX(IIG  ,JJG  ,KKG  ,Z_INDEX)+DIF_FX(IIG  ,JJG  ,KKG  ,Z_INDEX))
            CASE( 1); SOLID_PHASE_OUTPUT= (ADV_FX(IIG-1,JJG  ,KKG  ,Z_INDEX)+DIF_FX(IIG-1,JJG  ,KKG  ,Z_INDEX))
            CASE(-2); SOLID_PHASE_OUTPUT=-(ADV_FY(IIG  ,JJG  ,KKG  ,Z_INDEX)+DIF_FY(IIG  ,JJG  ,KKG  ,Z_INDEX))
            CASE( 2); SOLID_PHASE_OUTPUT= (ADV_FY(IIG  ,JJG-1,KKG  ,Z_INDEX)+DIF_FY(IIG  ,JJG-1,KKG  ,Z_INDEX))
            CASE(-3); SOLID_PHASE_OUTPUT=-(ADV_FZ(IIG  ,JJG  ,KKG  ,Z_INDEX)+DIF_FZ(IIG  ,JJG  ,KKG  ,Z_INDEX))
            CASE( 3); SOLID_PHASE_OUTPUT= (ADV_FZ(IIG  ,JJG  ,KKG-1,Z_INDEX)+DIF_FZ(IIG  ,JJG  ,KKG-1,Z_INDEX))
            END SELECT
         ELSEIF (Y_INDEX>0) THEN
            SELECT CASE(BC%IOR)
            CASE(-1); SOLID_PHASE_OUTPUT=-1._EB*DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                      (ADV_FX(IIG  ,JJG  ,KKG  ,1:N_TRACKED_SPECIES)+DIF_FX(IIG  ,JJG  ,KKG  ,1:N_TRACKED_SPECIES)) )
            CASE( 1); SOLID_PHASE_OUTPUT= 1._EB*DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                      (ADV_FX(IIG-1,JJG  ,KKG  ,1:N_TRACKED_SPECIES)+DIF_FX(IIG-1,JJG  ,KKG  ,1:N_TRACKED_SPECIES)) )
            CASE(-2); SOLID_PHASE_OUTPUT=-1._EB*DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                      (ADV_FY(IIG  ,JJG  ,KKG  ,1:N_TRACKED_SPECIES)+DIF_FY(IIG  ,JJG  ,KKG  ,1:N_TRACKED_SPECIES)) )
            CASE( 2); SOLID_PHASE_OUTPUT= 1._EB*DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                      (ADV_FY(IIG  ,JJG-1,KKG  ,1:N_TRACKED_SPECIES)+DIF_FY(IIG  ,JJG-1,KKG  ,1:N_TRACKED_SPECIES)) )
            CASE(-3); SOLID_PHASE_OUTPUT=-1._EB*DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                      (ADV_FZ(IIG  ,JJG  ,KKG  ,1:N_TRACKED_SPECIES)+DIF_FZ(IIG  ,JJG  ,KKG  ,1:N_TRACKED_SPECIES)) )
            CASE( 3); SOLID_PHASE_OUTPUT= 1._EB*DOT_PRODUCT( Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),&
                      (ADV_FZ(IIG  ,JJG  ,KKG-1,1:N_TRACKED_SPECIES)+DIF_FZ(IIG  ,JJG  ,KKG-1,1:N_TRACKED_SPECIES)) )
            END SELECT
         ENDIF
      ELSE
         UN  = -B1%U_NORMAL
         IF (Z_INDEX > 0) THEN
            Y_SPECIES = B1%ZZ_F(Z_INDEX)
            RHO_D_DYDN = B1%RHO_D_DZDN_F(Z_INDEX)
         ELSEIF (Y_INDEX > 0) THEN
            ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
            CALL GET_MASS_FRACTION(ZZ_GET,Y_INDEX,Y_SPECIES)
            RHO_D_DYDN = DOT_PRODUCT(Z2Y(Y_INDEX,1:N_TRACKED_SPECIES),B1%RHO_D_DZDN_F(1:N_TRACKED_SPECIES))
         ELSE
            Y_SPECIES = 1._EB
            RHO_D_DYDN = 0._EB
         ENDIF
         ! convention here is: inflow is positive (adds mass to domain), outflow is negative (subtracts mass)
         SOLID_PHASE_OUTPUT = Y_SPECIES*B1%RHO_F*UN - RHO_D_DYDN
      ENDIF
   CASE(65) ! WALL PRESSURE (takes optional FORCE_DIRECTION vector)
      IF (PRESENT(OPT_WALL_INDEX)) THEN
         IIG = BC%IIG
         JJG = BC%JJG
         KKG = BC%KKG
         ! quadratic extrapolation to surface pressure
         PR1 = RHO(IIG,JJG,KKG)*(H(IIG,JJG,KKG)-KRES(IIG,JJG,KKG))
         PR2 = PR1
         SELECT CASE(BC%IOR)
            CASE( 1)
               NVEC=(/ 1._EB,0._EB,0._EB/)
               Z1 = 0.5_EB*DX(IIG)
               Z2 = DX(IIG)+0.5_EB*DX(IIG+1)
               IC2 = CELL_INDEX(IIG+1,JJG,KKG)
               IF (.NOT.CELL(IC2)%SOLID) PR2 = RHO(IIG,JJG,KKG)*(H(IIG+1,JJG,KKG)-KRES(IIG+1,JJG,KKG))
            CASE(-1)
               NVEC=(/-1._EB,0._EB,0._EB/)
               Z1 = 0.5_EB*DX(IIG)
               Z2 = DX(IIG)+0.5_EB*DX(IIG-1)
               IC2 = CELL_INDEX(IIG-1,JJG,KKG)
               IF (.NOT.CELL(IC2)%SOLID) PR2 = RHO(IIG,JJG,KKG)*(H(IIG-1,JJG,KKG)-KRES(IIG-1,JJG,KKG))
            CASE( 2)
               NVEC=(/0._EB, 1._EB,0._EB/)
               Z1 = 0.5_EB*DY(JJG)
               Z2 = DY(JJG)+0.5_EB*DY(JJG+1)
               IC2 = CELL_INDEX(IIG,JJG+1,KKG)
               IF (.NOT.CELL(IC2)%SOLID) PR2 = RHO(IIG,JJG,KKG)*(H(IIG,JJG+1,KKG)-KRES(IIG,JJG+1,KKG))
            CASE(-2)
               NVEC=(/0._EB,-1._EB,0._EB/)
               Z1 = 0.5_EB*DY(JJG)
               Z2 = DY(JJG)+0.5_EB*DY(JJG-1)
               IC2 = CELL_INDEX(IIG,JJG-1,KKG)
               IF (.NOT.CELL(IC2)%SOLID) PR2 = RHO(IIG,JJG,KKG)*(H(IIG,JJG-1,KKG)-KRES(IIG,JJG-1,KKG))
            CASE( 3)
               NVEC=(/0._EB,0._EB, 1._EB/)
               Z1 = 0.5_EB*DZ(KKG)
               Z2 = DZ(KKG)+0.5_EB*DZ(KKG+1)
               IC2 = CELL_INDEX(IIG,JJG,KKG+1)
               IF (.NOT.CELL(IC2)%SOLID) PR2 = RHO(IIG,JJG,KKG)*(H(IIG,JJG,KKG+1)-KRES(IIG,JJG,KKG+1))
            CASE(-3)
               NVEC=(/0._EB,0._EB,-1._EB/)
               Z1 = 0.5_EB*DZ(KKG)
               Z2 = DZ(KKG)+0.5_EB*DZ(KKG-1)
               IC2 = CELL_INDEX(IIG,JJG,KKG-1)
               IF (.NOT.CELL(IC2)%SOLID) PR2 = RHO(IIG,JJG,KKG)*(H(IIG,JJG,KKG-1)-KRES(IIG,JJG,KKG-1))
         END SELECT

         PVEC = ( PR1 - (PR2-PR1)*Z1**2 / (Z2**2-Z1**2) ) * NVEC ! surface normal pressure force
      ELSEIF (PRESENT(OPT_CFACE_INDEX)) THEN
         NVEC = BC%NVEC
         ! find cut-cell adjacent to CFACE
         IND1 = CFA%CUT_FACE_IND1
         IND2 = CFA%CUT_FACE_IND2
         CALL GET_PRES_CFACE(PRESS,IND1,IND2,CFA)
         PVEC = PRESS * NVEC ! surface normal pressure force
      ENDIF

      SOLID_PHASE_OUTPUT = DOT_PRODUCT(PVEC,NVEC)

      IF(FROM_BNDF) RETURN

      IF (ASSOCIATED(DV)) THEN
         IF (NORM2(DV%DFVEC)>TWENTY_EPSILON_EB) THEN
            SOLID_PHASE_OUTPUT = -DOT_PRODUCT(PVEC,DV%DFVEC)
         ENDIF
      ENDIF

   CASE(66) ! VISCOUS STRESS WALL (takes optional FORCE_DIRECTION vector)
      IF (PRESENT(OPT_WALL_INDEX)) THEN
         SELECT CASE(BC%IOR)
            ! note: VEL_T does not follow a right hand rule, see user guide
            CASE( 1); NVEC=(/ 1._EB,0._EB,0._EB/); TVEC1=(/ 0._EB,1._EB,0._EB/); TVEC2=(/ 0._EB,0._EB,1._EB/)
            CASE(-1); NVEC=(/-1._EB,0._EB,0._EB/); TVEC1=(/ 0._EB,1._EB,0._EB/); TVEC2=(/ 0._EB,0._EB,1._EB/)
            CASE( 2); NVEC=(/0._EB, 1._EB,0._EB/); TVEC1=(/ 1._EB,0._EB,0._EB/); TVEC2=(/ 0._EB,0._EB,1._EB/)
            CASE(-2); NVEC=(/0._EB,-1._EB,0._EB/); TVEC1=(/ 1._EB,0._EB,0._EB/); TVEC2=(/ 0._EB,0._EB,1._EB/)
            CASE( 3); NVEC=(/0._EB,0._EB, 1._EB/); TVEC1=(/ 1._EB,0._EB,0._EB/); TVEC2=(/ 0._EB,1._EB,0._EB/)
            CASE(-3); NVEC=(/0._EB,0._EB,-1._EB/); TVEC1=(/ 1._EB,0._EB,0._EB/); TVEC2=(/ 0._EB,1._EB,0._EB/)
         END SELECT
         IIG = BC%IIG
         JJG = BC%JJG
         KKG = BC%KKG
         U_CELL = 0.5_EB*(U(IIG-1,JJG,KKG)+U(IIG,JJG,KKG))
         V_CELL = 0.5_EB*(V(IIG,JJG-1,KKG)+V(IIG,JJG,KKG))
         W_CELL = 0.5_EB*(W(IIG,JJG,KKG-1)+W(IIG,JJG,KKG))
         MU_WALL = MU_DNS(IIG,JJG,KKG)
      ELSEIF (PRESENT(OPT_CFACE_INDEX)) THEN
         NVEC = BC%NVEC
         ! right now VEL_T not defined for CFACEs
         TVEC1=(/ 0._EB,0._EB,0._EB/)
         TVEC2=(/ 0._EB,0._EB,0._EB/)
         ! find cut-cell adjacent to CFACE
         IND1 = CFA%CUT_FACE_IND1
         IND2 = CFA%CUT_FACE_IND2
         CALL GET_UVWGAS_CFACE(U_CELL,V_CELL,W_CELL,IND1,IND2,U,V,W,PREDFCT=1._EB)
         CALL GET_MUDNS_CFACE(MU_WALL,IND1,IND2)
         ICC = CUT_FACE(IND1)%CELL_LIST(2,LOW_IND,IND2)
         IIG = CUT_CELL(ICC)%IJK(1)
         JJG = CUT_CELL(ICC)%IJK(2)
         KKG = CUT_CELL(ICC)%IJK(3)
      ENDIF

      IF (PRESENT(OPT_WALL_INDEX) .OR. PRESENT(OPT_CFACE_INDEX)) THEN
         DN  = 1._EB/B1%RDN
         ! velocity vector in the centroid of the gas (cut) cell
         VEL_CELL = (/U_CELL,V_CELL,W_CELL/)
         ! velocity vector of the surface
         IF (SF%VELOCITY_BC_INDEX == FREE_SLIP_BC) THEN
            ! U_NORMAL velocity in Normal direction, same tangential velocities as VEL_CELL:
            VEL_WALL = -B1%U_NORMAL*NVEC + ( VEL_CELL - DOT_PRODUCT(VEL_CELL,NVEC)*NVEC )
         ELSE
            VEL_WALL = -B1%U_NORMAL*NVEC + SF%VEL_T(1)*TVEC1 + SF%VEL_T(2)*TVEC2
         ENDIF
         RHO_WALL = B1%RHO_F
         CALL TAU_WALL_IJ(TAU_IJ,SVEC,VEL_CELL,VEL_WALL,NVEC,DN,D(IIG,JJG,KKG),MU_WALL,RHO_WALL,SF%ROUGHNESS)
         DO I=1,3
            FVEC(I) = DOT_PRODUCT(TAU_IJ(I,:),NVEC(:))
         ENDDO
         SOLID_PHASE_OUTPUT = DOT_PRODUCT(FVEC,SVEC)
         IF (FROM_BNDF) RETURN
         IF (ASSOCIATED(DV)) THEN
            IF (NORM2(DV%DFVEC)>TWENTY_EPSILON_EB) SOLID_PHASE_OUTPUT = DOT_PRODUCT(FVEC,DV%DFVEC)
         ENDIF
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF

   CASE(67) ! WALL PRESSURE TEST (takes optional FORCE_DIRECTION vector)
      IF (PRESENT(OPT_WALL_INDEX)) THEN
         SELECT CASE(BC%IOR)
            CASE( 1); NVEC=(/ 1._EB,0._EB,0._EB/)
            CASE(-1); NVEC=(/-1._EB,0._EB,0._EB/)
            CASE( 2); NVEC=(/0._EB, 1._EB,0._EB/)
            CASE(-2); NVEC=(/0._EB,-1._EB,0._EB/)
            CASE( 3); NVEC=(/0._EB,0._EB, 1._EB/)
            CASE(-3); NVEC=(/0._EB,0._EB,-1._EB/)
         END SELECT
         IIG = BC%IIG
         JJG = BC%JJG
         KKG = BC%KKG
         PVEC = RHO(IIG,JJG,KKG)*H(IIG,JJG,KKG) * NVEC ! surface normal pressure force
      ELSEIF (PRESENT(OPT_CFACE_INDEX)) THEN
         NVEC = BC%NVEC
         ! find cut-cell adjacent to CFACE
         IND1 = CFA%CUT_FACE_IND1
         IND2 = CFA%CUT_FACE_IND2
         CALL GET_PRES_CFACE_TEST(PRESS,IND1,IND2,CFA)
         PVEC = PRESS * NVEC ! surface normal pressure force
      ENDIF

      SOLID_PHASE_OUTPUT = DOT_PRODUCT(PVEC,NVEC)

      IF(FROM_BNDF) RETURN

      IF (ASSOCIATED(DV)) THEN
         IF (NORM2(DV%DFVEC)>TWENTY_EPSILON_EB) THEN
            SOLID_PHASE_OUTPUT = -DOT_PRODUCT(PVEC,DV%DFVEC)
         ENDIF
      ENDIF

   CASE(68) ! LEVEL SET
      IF (ASSOCIATED(B2)) THEN
         SOLID_PHASE_OUTPUT = B2%PHI_LS
      ELSE
         SOLID_PHASE_OUTPUT = 0._EB
      ENDIF

   CASE(69) ! WALL ENTHALPY
      SOLID_PHASE_OUTPUT = 0._EB
      IF (SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
         IF (ONE_D%PYROLYSIS_MODEL==PYROLYSIS_PREDICTED .OR. SF%HT_DIM>1) THEN
            NWP = SUM(ONE_D%N_LAYER_CELLS(1:ONE_D%N_LAYERS))
            X0 = SUM(ONE_D%LAYER_THICKNESS)
         ELSE
            NWP = SF%N_CELLS_INI
            X0 = ONE_D%X(NWP)
         ENDIF
         DO I=1,NWP
            SELECT CASE (SF%GEOMETRY)
               CASE DEFAULT
                  VOL = B1%AREA*(ONE_D%X(I)-ONE_D%X(I-1))
               CASE (SURF_CYLINDRICAL)
                  VOL = SF%LENGTH*PI*((SF%INNER_RADIUS+X0-ONE_D%X(I-1))**2-(SF%INNER_RADIUS+X0-ONE_D%X(I))**2)
               CASE (SURF_INNER_CYLINDRICAL)
                  VOL = SF%LENGTH*PI*((SF%INNER_RADIUS+ONE_D%X(I))**2-(SF%INNER_RADIUS+ONE_D%X(I-1))**2)
               CASE (SURF_SPHERICAL)
                  VOL = FOTHPI*((X0-ONE_D%X(I-1))**3-(X0-ONE_D%X(I))**3)
            END SELECT
            H_MATL_LOOP: DO J=1,ONE_D%N_MATL
               IF (ONE_D%MATL_COMP(J)%RHO(I)<=TWENTY_EPSILON_EB) CYCLE H_MATL_LOOP
               ML  => MATERIAL(ONE_D%MATL_INDEX(J))
               ITMP = INT(ONE_D%TMP(I))
               SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT + ONE_D%MATL_COMP(J)%RHO(I)*VOL * &
                              (ML%H(ITMP)+(ONE_D%TMP(I)-REAL(ITMP,EB))*(ML%H(MIN(I_MAX_TEMP,ITMP+1))-ML%H(ITMP)))
            ENDDO H_MATL_LOOP
         ENDDO
         IF (PRESENT(OPT_LP_INDEX)) SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT*LP%PWT
         SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT * 0.001_EB
      ENDIF

   CASE(70) ! SUBSTEPS
      SOLID_PHASE_OUTPUT = REAL(B1%N_SUBSTEPS,EB)

   CASE(71) ! EFFECTIVE HEAT TRANSFER COEFFICIENT
      DTMP = TMP(BC%IIG,BC%JJG,BC%KKG)-0.5_EB*(B1%TMP_F_OLD+B1%TMP_F)
      IF (ABS(DTMP)>TWENTY_EPSILON_EB .AND. ABS(B1%Q_CON_F)>TWENTY_EPSILON_EB) THEN
         SOLID_PHASE_OUTPUT = B1%Q_CON_F/DTMP
      ELSE
         SOLID_PHASE_OUTPUT = B1%HEAT_TRANS_COEF
      ENDIF
   CASE(72) ! SCALING HEAT FLUX
      SOLID_PHASE_OUTPUT = B1%Q_IN_SMOOTH*0.001_EB

   CASE(73) ! VEGETATION FUEL TYPE
      SOLID_PHASE_OUTPUT = SF%VEG_LSET_FUEL_INDEX

   CASE(74) ! SOLID MASS FRACTION
      SOLID_PHASE_OUTPUT = 0._EB
      X0 = 0._EB
      DO NN=1,ONE_D%N_MATL
         X0 = X0 + ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)
         IF (MATL_INDEX==ONE_D%MATL_INDEX(NN)) SOLID_PHASE_OUTPUT = ONE_D%MATL_COMP(NN)%RHO(I_DEPTH)
      ENDDO
      SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT / (X0 + TWENTY_EPSILON_EB)

   CASE(75) ! SOLID ENTHALPY
      SOLID_PHASE_OUTPUT = 0._EB
      SH_MATL_LOOP: DO J=1,ONE_D%N_MATL
         IF (ONE_D%MATL_COMP(J)%RHO(I_DEPTH)<=TWENTY_EPSILON_EB) CYCLE SH_MATL_LOOP
         ITMP = INT(ONE_D%TMP(I_DEPTH))
         ML  => MATERIAL(ONE_D%MATL_INDEX(J))
         SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT + ONE_D%MATL_COMP(J)%RHO(I_DEPTH) * &
                              (ML%H(ITMP)+(ONE_D%TMP(I_DEPTH)-REAL(ITMP,EB))*(ML%H(MIN(I_MAX_TEMP,ITMP+1))-ML%H(ITMP)))
      ENDDO SH_MATL_LOOP
      SOLID_PHASE_OUTPUT = SOLID_PHASE_OUTPUT * 0.001_EB

   CASE(76) ! CONVECTIVE HEAT FLUX GAUGE
      IF (PY%HEAT_TRANSFER_COEFFICIENT>=0._EB) THEN
         Q_CON = PY%HEAT_TRANSFER_COEFFICIENT*(TMP(BC%IIG,BC%JJG,BC%KKG)-PY%GAUGE_TEMPERATURE)
      ELSE
         Q_CON = B1%Q_CON_F + B1%HEAT_TRANS_COEF*(B1%TMP_F-PY%GAUGE_TEMPERATURE)
      ENDIF
      SOLID_PHASE_OUTPUT = Q_CON*0.001_EB

   CASE(77) ! CONVECTIVE HEAT TRANSFER REGIME
      SOLID_PHASE_OUTPUT = 0._EB
      IF (ASSOCIATED(B2)) SOLID_PHASE_OUTPUT = B2%HEAT_TRANSFER_REGIME
   CASE(78) ! SURFACE OXYGEN MASS FRACTION
      SOLID_PHASE_OUTPUT = 0._EB
      IF (ASSOCIATED(B2)) SOLID_PHASE_OUTPUT = B2%Y_O2_F
   CASE(79) ! SURFACE OXYGEN ITERATIONS
      SOLID_PHASE_OUTPUT = 0._EB
      IF (ASSOCIATED(B2)) SOLID_PHASE_OUTPUT = B2%Y_O2_ITER
   CASE(80) ! OXIDATIVE HRRPUA
      SOLID_PHASE_OUTPUT = B1%Q_DOT_O2_PP*0.001_EB
   CASE(81) ! SOLID OXYGEN MASS FRACTION
      SOLID_PHASE_OUTPUT = 0._EB
      IF (ASSOCIATED(B2) .AND. MATL_INDEX>0) THEN
         ML => MATERIAL(MATL_INDEX)
         ! for the moment this assumes there is only one char reaction
         IF (ML%N_O2(1)>0._EB) THEN
            DEPTH = 0.5_EB*(ONE_D%X(I_DEPTH-1)+ONE_D%X(I_DEPTH))
            ASH_DEPTH = 0._EB
            IF (TEST_NEW_CHAR_MODEL) ASH_DEPTH = ONE_D%X(B2%I_ASH_DEPTH-1)
            SOLID_PHASE_OUTPUT = B2%Y_O2_F*EXP(-MAX(0._EB,DEPTH-ASH_DEPTH)/(TWENTY_EPSILON_EB+ML%GAS_DIFFUSION_DEPTH(1)))
         ENDIF
      ENDIF
   CASE(82) ! BLOWING CORRECTION
      SOLID_PHASE_OUTPUT = 0._EB
      IF (ASSOCIATED(B2)) SOLID_PHASE_OUTPUT = B2%BLOWING_CORRECTION
   CASE(90:92) ! FIRE ARRIVAL TIME, FIRE RESIDENCE TIME, LS SPREAD RATE
      IF (PRESENT(OPT_WALL_INDEX)) THEN
         OUTPUT_INDEX = OPT_WALL_INDEX
      ELSEIF (PRESENT(OPT_CFACE_INDEX)) THEN
         OUTPUT_INDEX = OPT_CFACE_INDEX-INTERNAL_CFACE_CELLS_LB+N_INTERNAL_WALL_CELLS+N_EXTERNAL_WALL_CELLS
      ENDIF
      SELECT CASE(INDX)
         CASE(90); SOLID_PHASE_OUTPUT = FIRE_ARRIVAL_TIME(OUTPUT_INDEX)
         CASE(91); SOLID_PHASE_OUTPUT = FIRE_RESIDENCE_TIME(OUTPUT_INDEX)
         CASE(92); SOLID_PHASE_OUTPUT = LS_SPREAD_RATE(OUTPUT_INDEX)
      END SELECT

   CASE(100) ! CONDENSATION HEAT FLUX
      SOLID_PHASE_OUTPUT = B1%Q_CONDENSE * 0.001_EB

   CASE(101) ! TANGENTIAL VELOCITY
      IF (ASSOCIATED(B1)) SOLID_PHASE_OUTPUT = B1%U_TANG

END SELECT SOLID_PHASE_SELECT

END FUNCTION SOLID_PHASE_OUTPUT


END MODULE GET_DATA
