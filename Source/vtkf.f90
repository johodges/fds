!> \brief Routines for handling vtk output

MODULE VTK_FDS_INTERFACE

USE PRECISION_PARAMETERS
USE MESH_VARIABLES
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE OUTPUT_DATA
USE OUTPUT_CLOCKS
USE DEVICE_VARIABLES
USE PROPERTY_DATA
USE DUMP, ONLY : DUMP_MESH_SPREADSHEET_OUTPUTS,GAS_PHASE_OUTPUT,SOLID_PHASE_OUTPUT,PARTICLE_OUTPUT,&
                 GET_SMOKE3D_QQ,GET_GEOMSIZES,GET_GEOMINFO,GET_GEOMVALS,GETSLICEDIR,&
                 COMPUTE_PARTICLE_FLUXES,GET_SLICE_QUANTITY,GET_SLICE_CORNER_WEIGHTS,DRY
#ifdef WITH_HDF5
USE HDF5
#endif
USE MPI_F08
USE COMPLEX_GEOMETRY, ONLY : WRITE_GEOM,WRITE_GEOM_ALL,CC_FGSC,CC_IDCF,CC_IDCC,CC_UNKZ,CC_UNKF,CC_FTYPE_RCGAS,&
                             CC_FTYPE_CFGAS,CC_FTYPE_CFINB,CC_SOLID,CC_CGSC,CC_CUTCFE,TRIANGULATE,&
                             CC_VGSC,CC_GASPHASE,MAKE_UNIQUE_VERT_ARRAY,AVERAGE_FACE_VALUES

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

! Pointers used by the VTKHDF output drivers below (mirrors the set carried by DUMP)

TYPE (LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (PROPERTY_TYPE), POINTER :: PY
TYPE (SLICE_TYPE), POINTER :: SL
TYPE (BOUNDARY_FILE_TYPE), POINTER :: BF

! HDF5 identifiers and per-mesh sizing for the VTKHDF output files

#ifdef WITH_HDF5
! Total number of obst and geom patches in each mesh
! 2*NMESHES arrays. odds include OBST patches, evens include GEOM patches
INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: NCELLS_VTK, NPOINTS_VTK, NCONNECTIONS_VTK

INTEGER(HID_T) :: HDF_SM3D_FILE_ID, HDF_SM3D_PLIST_ID, HDF_SM3D_CRP_LIST
INTEGER(HID_T) :: HDF_SM3D_G1,HDF_SM3D_G2,HDF_SM3D_G3,HDF_SM3D_G4,HDF_SM3D_G5,HDF_SM3D_G6,HDF_SM3D_G7 ! Group identifier
INTEGER(HID_T) :: HDF_SM3D_COUNTER=0

INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: HDF_SLCF_FILE_ID, HDF_SLCF_PLIST_ID, HDF_SLCF_CRP_LIST
INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: HDF_SLCF_G1,HDF_SLCF_G2,HDF_SLCF_G3,HDF_SLCF_G4
INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: HDF_SLCF_G5,HDF_SLCF_G6,HDF_SLCF_G7 ! Group identifier
INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: HDF_SLCF_COUNTER
! Per-mesh point and cell counts for each unique slice plane
INTEGER(IB32), ALLOCATABLE, DIMENSION(:,:) :: HDF_SLCF_G1_NCELLS, HDF_SLCF_G1_NPOINTS

INTEGER(HID_T) :: HDF_BNDF_FILE_ID, HDF_BNDF_PLIST_ID, HDF_BNDF_CRP_LIST
INTEGER(HID_T) :: HDF_BNDF_G1,HDF_BNDF_G2,HDF_BNDF_G3,HDF_BNDF_G4,HDF_BNDF_G5,HDF_BNDF_G6,HDF_BNDF_G7 ! Group identifier
INTEGER(HID_T) :: HDF_BNDF_COUNTER=0

INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: HDF_PART_FILE_ID, HDF_PART_PLIST_ID, HDF_PART_CRP_LIST
INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: HDF_PART_G1,HDF_PART_G2,HDF_PART_G3,HDF_PART_G4
INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: HDF_PART_G5,HDF_PART_G6,HDF_PART_G7
INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: HDF_PART_COUNTER

! One collective data-transfer property list is shared by every VTKHDF write.  It is
! created on first use and lives for the run.  Creating one per dataset access (the
! previous behaviour) leaked a property list on every write of every quantity of every
! mesh at every output time.

INTEGER(HID_T) :: VTK_DXPL_ID=-1_HID_T


! Chunk sizes are clamped into this range (in elements).  Chunks of one element, which
! is what the per-time-step scalar datasets used to get, cost a B-tree entry and a
! filter-pipeline pass each; chunks of a whole mesh blow past the chunk cache.

INTEGER(HSIZE_T), PARAMETER :: VTK_CHUNK_MIN=1024_HSIZE_T, VTK_CHUNK_MAX=1048576_HSIZE_T

! Chunk length (in points) for the particle arrays.  Their length changes at every output
! time, and they are created at the first one, when a class typically has no particles at
! all -- deriving the chunk from that gave every particle dataset the clamp floor above,
! so a run ended up with ten thousand chunks of 1024 points, each its own B-tree entry and
! its own filter pipeline pass.

INTEGER, PARAMETER :: VTK_PART_CHUNK=4096

! Metadata block size (bytes) for the VTKHDF files.  HDF5 allocates metadata in blocks
! of this size, and whatever is left over at the end of the last block stays in the file
! as slack, so an oversized value costs every file a fixed lump of dead space: at 1 MB
! the Verification/VTK cases carried more slack than data.

INTEGER(HSIZE_T), PARAMETER :: VTK_META_BLOCK_SIZE=65536_HSIZE_T

! Open dataset handles, keyed by the group they live in and their name.  Every append
! used to reach its dataset by name at every output time, which costs an H5Lexists and
! an H5Dopen -- both collective, since the file access property list asks for collective
! metadata operations -- and an H5Dclose that evicts the dataset's chunk from the chunk
! cache, forcing a filter pipeline round trip on the next append.  Holding the handles
! open for the life of the file removes all three.  VTK_DSET_CACHE_PURGE closes them.

INTEGER, PARAMETER :: VTK_DSET_CACHE_MAX=2048
INTEGER, PARAMETER :: VTK_DSET_NAME_LEN=64
INTEGER(HID_T) :: VTK_DSET_CACHE_GROUP(VTK_DSET_CACHE_MAX)=-1_HID_T
INTEGER(HID_T) :: VTK_DSET_CACHE_ID(VTK_DSET_CACHE_MAX)=-1_HID_T
CHARACTER(VTK_DSET_NAME_LEN) :: VTK_DSET_CACHE_NAME(VTK_DSET_CACHE_MAX)=''
INTEGER :: VTK_DSET_CACHE_N=0

#endif

PUBLIC WRITE_PARAVIEW_STATE_FILE,EXCHANGE_NSLICE_INFO,EXCHANGE_NPATCH_INFO
#ifdef WITH_HDF5
! Everything the rest of FDS needs from the VTKHDF writers.  The HDF5 plumbing below
! (PARALLEL_INIT_*, PARALLEL_WRITE_*, the OPEN/CLOSE/EXTEND helpers) is private now that
! all of the VTKHDF output lives in this module.

PUBLIC WRITE_VTKHDF_GEOM_FILE,DUMP_VTK_MESH_OUTPUTS_SERIES,EXCHANGE_NOBST_INFO,&
       CLOSE_VTKHDF_BNDF,CLOSE_VTKHDF_SMOKE3D,CLOSE_VTKHDF_SLICE,CLOSE_VTKHDF_PART
#endif

CONTAINS



#ifdef WITH_HDF5

!> \brief Build the VTK hexahedral cells spanning an entire mesh

!>

!> \param NM Mesh number

!> \param NC Number of cells (out)

!> \param NP Number of points (out)

!> \param VERTICES Point coordinates, allocated here (out)

!> \param CONNECT Point indices of each cell, allocated here (out)

!> \param OFFSETS Where each cell starts in CONNECT, allocated here (out)

!> \param VTKC_TYPE VTK cell type of each cell, allocated here (out)

!>

!> Connectivity is local to the mesh, as it is for every VTKHDF output, so the

!> mesh becomes one partition of the assembled grid.


SUBROUTINE BUILD_VTK_GAS_PHASE_GEOMETRY2(NM, &
                                        NC, NP, VERTICES, CONNECT, OFFSETS, VTKC_TYPE)

INTEGER :: NX, NY, NZ, NC, NP, I, J, K, IFACT, JFACT, KFACT
REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: VERTICES
INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS
INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: VTKC_TYPE
INTEGER, INTENT(IN) :: NM

NX = SIZE(MESHES(NM)%X)
NY = SIZE(MESHES(NM)%Y)
NZ = SIZE(MESHES(NM)%Z)

NP = NX*NY*NZ
NC = (NX-1)*(NY-1)*(NZ-1)

! Fill point data
ALLOCATE(VERTICES(3,NP))
IFACT = 1
DO K = 0, NZ-1
   DO J = 0, NY-1
      DO I = 0, NX-1
         VERTICES(1, IFACT)=REAL(MESHES(NM)%X(I),FB)
         VERTICES(2, IFACT)=REAL(MESHES(NM)%Y(J),FB)
         VERTICES(3,IFACT)=REAL(MESHES(NM)%Z(K),FB)
         IFACT = IFACT + 1
      ENDDO
   ENDDO
ENDDO

! Fill cell data
ALLOCATE(CONNECT(NC*8))
DO I = 1, NX-1
   IFACT = (I-1)
   DO J = 1, NY-1
      JFACT = (J-1)*(NX-1)
      DO K = 1, NZ-1
         KFACT = (K - 1)*(NY-1)*(NX-1)
         CONNECT((IFACT+JFACT+KFACT)*8+1) = (K-1)*(NY*NX) + (J-1)*NX + I-1
         CONNECT((IFACT+JFACT+KFACT)*8+2) = (K-1)*(NY*NX) + (J-1)*NX + I
         CONNECT((IFACT+JFACT+KFACT)*8+3) = (K-1)*(NY*NX) + (J)*NX + I-1
         CONNECT((IFACT+JFACT+KFACT)*8+4) = (K-1)*(NY*NX) + (J)*NX + I
         CONNECT((IFACT+JFACT+KFACT)*8+5) = (K)*(NY*NX) + (J-1)*NX + I-1
         CONNECT((IFACT+JFACT+KFACT)*8+6) = (K)*(NY*NX) + (J-1)*NX + I
         CONNECT((IFACT+JFACT+KFACT)*8+7) = (K)*(NY*NX) + (J)*NX + I-1
         CONNECT((IFACT+JFACT+KFACT)*8+8) = (K)*(NY*NX) + (J)*NX + I
      ENDDO
   ENDDO
ENDDO

ALLOCATE(OFFSETS(NC+1))
ALLOCATE(VTKC_TYPE(NC))

OFFSETS(1) = 0
DO I=1,NC
   OFFSETS(I+1) = (I)*8_IB32
   VTKC_TYPE(I) = 11_IB8
ENDDO

ENDSUBROUTINE BUILD_VTK_GAS_PHASE_GEOMETRY2



!> \brief Build the VTK cells covering one slice within one mesh



!>



!> \param NM Mesh number



!> \param SL Slice to build the geometry for



!> \param NTSL Terrain slice counter, used to index K_AGL_SLICE



!> \param NC Number of cells (out)



!> \param NP Number of points (out)



!> \param VERTICES Point coordinates, allocated here (out)



!> \param CONNECT Point indices of each cell, allocated here (out)



!> \param OFFSETS Where each cell starts in CONNECT, allocated here (out)



!> \param VTKC_TYPE VTK cell type of each cell, allocated here (out)



!>



!> A slice with a degenerate direction is a sheet of quadrilaterals; one with none



!> is a block of hexahedra.  A terrain slice follows the ground, so its points take



!> their height from K_AGL_SLICE rather than from the slice plane.




SUBROUTINE BUILD_VTK_SLICE_GEOMETRY2(NM, SL, NTSL, &
                                        NC, NP, VERTICES, CONNECT, OFFSETS, VTKC_TYPE)

INTEGER :: NX, NY, NZ, NC, NP, I, J, K, IFACT, JFACT, KFACT, KTS
INTEGER :: I1,I2,J1,J2,K1,K2,L1,L2,N1,N2
REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: VERTICES
INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS
INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: VTKC_TYPE
INTEGER, INTENT(IN) :: NM, NTSL
TYPE(SLICE_TYPE), POINTER, INTENT(IN) :: SL

I1  = SL%I1
I2  = SL%I2
J1  = SL%J1
J2  = SL%J2
K1  = SL%K1
K2  = SL%K2

NX = I2 + 1 - I1
NY = J2 + 1 - J1
NZ = K2 + 1 - K1

NP = NX*NY*NZ
NC = MAX((NX-1),1)*MAX((NY-1),1)*MAX((NZ-1),1)

! Fill point data
ALLOCATE(VERTICES(3,NP))
IF (I2-I1==0 .OR. J2-J1==0 .OR. K2-K1==0) THEN
   IFACT = 1
   IF (I2-I1==0) THEN
      DO K = K1, K2
         DO J = J1, J2
            VERTICES(1,IFACT)=REAL(MESHES(NM)%X(I1),FB)
            VERTICES(2,IFACT)=REAL(MESHES(NM)%Y(J),FB)
            VERTICES(3,IFACT)=REAL(MESHES(NM)%Z(K),FB)
            IFACT = IFACT + 1
         ENDDO
      ENDDO
   ELSEIF (J2-J1==0) THEN
      DO K = K1, K2
         DO I = I1, I2
            VERTICES(1,IFACT)=REAL(MESHES(NM)%X(I),FB)
            VERTICES(2,IFACT)=REAL(MESHES(NM)%Y(J1),FB)
            VERTICES(3,IFACT)=REAL(MESHES(NM)%Z(K),FB)
            IFACT = IFACT + 1
         ENDDO
      ENDDO
   ELSEIF (K2-K1==0) THEN
      DO J = J1, J2
         DO I = I1, I2
            VERTICES(1,IFACT)=REAL(MESHES(NM)%X(I),FB)
            VERTICES(2,IFACT)=REAL(MESHES(NM)%Y(J),FB)
            IF (SL%TERRAIN_SLICE) THEN
               KTS = MESHES(NM)%K_AGL_SLICE(I,J,NTSL)
               VERTICES(3,IFACT)=REAL(MESHES(NM)%Z(KTS),FB)
            ELSE
               VERTICES(3,IFACT)=REAL(MESHES(NM)%Z(K1),FB)
            ENDIF
            IFACT = IFACT + 1
         ENDDO
      ENDDO
   ENDIF
ELSE
   IFACT = 1
   DO K = K1, K2
      DO J = J1, J2
         DO I = I1, I2
            VERTICES(1,IFACT)=REAL(MESHES(NM)%X(I),FB)
            VERTICES(2,IFACT)=REAL(MESHES(NM)%Y(J),FB)
            VERTICES(3,IFACT)=REAL(MESHES(NM)%Z(K),FB)
            IFACT = IFACT + 1
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Fill cell data
ALLOCATE(OFFSETS(NC+1))
ALLOCATE(VTKC_TYPE(NC))
IF (I2-I1==0 .OR. J2-J1==0 .OR. K2-K1==0) THEN
   !2-D slice
   ALLOCATE(CONNECT(NC*4))
   
   IF (I2-I1==0) THEN
      L1=J1 ; L2=J2
      N1=K1 ; N2=K2;
   ELSEIF (J2-J1==0) THEN
      L1=I1 ; L2=I2
      N1=K1 ; N2=K2;
   ELSEIF (K2-K1==0) THEN
      L1=I1 ; L2=I2
      N1=J1 ; N2=J2;
   ENDIF
   
   IFACT = 0
   DO I = 1, (L2-L1)
      DO J = 1, (N2-N1)
         CONNECT((IFACT)*4+1) = (L2-L1+1)*(J-1)+(I-1)
         CONNECT((IFACT)*4+2) = (L2-L1+1)*(J-1)+(I)
         CONNECT((IFACT)*4+3) = (L2-L1+1)*(J)+(I-1)
         CONNECT((IFACT)*4+4) = (L2-L1+1)*(J)+(I)
         IFACT = IFACT+1
      ENDDO
   ENDDO
   OFFSETS(1) = 0
   DO I=1,NC
      OFFSETS(I+1) = (I)*4_IB32
      VTKC_TYPE(I) = 8_IB8
   ENDDO
   
ELSE
   !3-D slice
   ALLOCATE(CONNECT(NC*8))
   DO I = 1, NX-1
      IFACT = (I-1)
      DO J = 1, NY-1
         JFACT = (J-1)*(NX-1)
         DO K = 1, NZ-1
            KFACT = (K - 1)*(NY-1)*(NX-1)
            CONNECT((IFACT+JFACT+KFACT)*8+1) = (K-1)*(NY*NX) + (J-1)*NX + I-1
            CONNECT((IFACT+JFACT+KFACT)*8+2) = (K-1)*(NY*NX) + (J-1)*NX + I
            CONNECT((IFACT+JFACT+KFACT)*8+3) = (K-1)*(NY*NX) + (J)*NX + I-1
            CONNECT((IFACT+JFACT+KFACT)*8+4) = (K-1)*(NY*NX) + (J)*NX + I
            CONNECT((IFACT+JFACT+KFACT)*8+5) = (K)*(NY*NX) + (J-1)*NX + I-1
            CONNECT((IFACT+JFACT+KFACT)*8+6) = (K)*(NY*NX) + (J-1)*NX + I
            CONNECT((IFACT+JFACT+KFACT)*8+7) = (K)*(NY*NX) + (J)*NX + I-1
            CONNECT((IFACT+JFACT+KFACT)*8+8) = (K)*(NY*NX) + (J)*NX + I
         ENDDO
      ENDDO
   ENDDO
   
   OFFSETS(1) = 0
   DO I=1,NC
      OFFSETS(I+1) = (I)*8_IB32
      VTKC_TYPE(I) = 11_IB8
   ENDDO
   
ENDIF
         
ENDSUBROUTINE BUILD_VTK_SLICE_GEOMETRY2


!> \brief Build the VTK quadrilaterals covering one boundary patch


!>


!> \param NM Mesh number


!> \param PA Patch to build the geometry for


!> \param NCELLS Number of cells (out)


!> \param NPOINTS Number of points (out)


!> \param X_PTS Point x coordinates, allocated here (out)


!> \param Y_PTS Point y coordinates, allocated here (out)


!> \param Z_PTS Point z coordinates, allocated here (out)


!> \param CONNECT Point indices of each cell, allocated here (out)


!> \param OFFSETS Where each cell starts in CONNECT, allocated here (out)


!> \param VTKC_TYPE VTK cell type of each cell, allocated here (out)


!>


!> The patch lies in the plane normal to PA%IOR, so which pair of mesh directions


!> the quadrilaterals span depends on that orientation.



SUBROUTINE BUILD_VTK_SOLID_PHASE_GEOMETRY(NM, PA, &
                                          NCELLS, NPOINTS, X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)

INTEGER :: NX,NY,NZ,NCELLS,NPOINTS,I,J,K,IFACT,L1,L2,N1,N2
REAL(FB), ALLOCATABLE, DIMENSION(:) :: X_PTS, Y_PTS, Z_PTS
INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS
INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: VTKC_TYPE
INTEGER, INTENT(IN) :: NM
TYPE(PATCH_TYPE), POINTER, INTENT(IN) :: PA

! Fill point data
SELECT CASE(ABS(PA%IOR))
   CASE(1) ; L1=PA%JG1 ; L2=PA%JG2 ; N1=PA%KG1 ; N2=PA%KG2 ; NX=1; NY=(L2-L1); NZ=(N2-N1)
   CASE(2) ; L1=PA%IG1 ; L2=PA%IG2 ; N1=PA%KG1 ; N2=PA%KG2 ; NX=(L2-L1); NY=1; NZ=(N2-N1)
   CASE(3) ; L1=PA%IG1 ; L2=PA%IG2 ; N1=PA%JG1 ; N2=PA%JG2 ; NX=(L2-L1); NY=(N2-N1); NZ=1
END SELECT
NPOINTS = (L2-L1+2)*(N2-N1+2)
NCELLS = (L2-L1+1)*(N2-N1+1)
ALLOCATE(X_PTS(NPOINTS))
ALLOCATE(Y_PTS(NPOINTS))
ALLOCATE(Z_PTS(NPOINTS))
ALLOCATE(CONNECT(NCELLS*4))
ALLOCATE(OFFSETS(NCELLS))
ALLOCATE(VTKC_TYPE(NCELLS))
IFACT = 1
SELECT CASE(ABS(PA%IOR))
   CASE(1)
      DO K = PA%K1,PA%K2
         DO J = PA%J1,PA%J2
            X_PTS(IFACT)=REAL(MESHES(NM)%X(PA%I1),FB)
            Y_PTS(IFACT)=REAL(MESHES(NM)%Y(J),FB)
            Z_PTS(IFACT)=REAL(MESHES(NM)%Z(K),FB)
            IFACT = IFACT + 1
         ENDDO
      ENDDO
   CASE(2)
      DO K = PA%K1,PA%K2
         DO I = PA%I1,PA%I2
            X_PTS(IFACT)=REAL(MESHES(NM)%X(I),FB)
            Y_PTS(IFACT)=REAL(MESHES(NM)%Y(PA%J1),FB)
            Z_PTS(IFACT)=REAL(MESHES(NM)%Z(K),FB)
            IFACT = IFACT + 1
         ENDDO
      ENDDO
   CASE(3)
      DO J = PA%J1,PA%J2
         DO I = PA%I1,PA%I2
            X_PTS(IFACT)=REAL(MESHES(NM)%X(I),FB)
            Y_PTS(IFACT)=REAL(MESHES(NM)%Y(J),FB)
            Z_PTS(IFACT)=REAL(MESHES(NM)%Z(PA%K1),FB)
            IFACT = IFACT + 1
         ENDDO
      ENDDO
END SELECT

! Fill cell data

IFACT = 0
DO I = 1, (L2-L1+1)
   DO J = 1, (N2-N1+1)
      CONNECT((IFACT)*4+1) = (L2-L1+2)*(J-1)+(I-1)
      CONNECT((IFACT)*4+2) = (L2-L1+2)*(J-1)+(I)
      CONNECT((IFACT)*4+3) = (L2-L1+2)*(J)+(I-1)
      CONNECT((IFACT)*4+4) = (L2-L1+2)*(J)+(I)
      IFACT = IFACT+1
   ENDDO
ENDDO

DO I=1,NCELLS
   OFFSETS(I) = (I)*4_IB32
   VTKC_TYPE(I) = 8_IB8
ENDDO

ENDSUBROUTINE BUILD_VTK_SOLID_PHASE_GEOMETRY

!> \brief Report whether one face of an obstruction belongs in the VTKHDF geometry
!>
!> The geometry file is built from the obstruction list rather than from the boundary
!> PATCHes.  A PATCH is only created where an obstruction face abuts a gas cell of the
!> same mesh, so a face lying in a mesh interface plane has none and used to be left
!> out, opening a seam in the rendered solid.  EXPOSED_FACE_INDEX marks the faces that
!> abut gas or a mesh boundary and is what the .smv file hands Smokeview, so keying off
!> it here makes the two renderings agree.  Faces permanently covered by another
!> obstruction stay unset and remain hidden.
!>
!> \param OB Obstruction
!> \param IOR Orientation index of the face

LOGICAL FUNCTION VTK_OBST_FACE_DRAWN(OB,IOR)

TYPE(OBSTRUCTION_TYPE), INTENT(IN) :: OB
INTEGER, INTENT(IN) :: IOR
INTEGER :: FI

VTK_OBST_FACE_DRAWN = .FALSE.
IF (OB%HIDDEN) RETURN
FI = ABS(IOR)*2 ; IF (IOR<0) FI = FI-1
IF (OB%EXPOSED_FACE_INDEX(FI)/=1) RETURN

! A mesh boundary can clip an obstruction to zero cells thick, leaving a copy whose
! faces coincide with ones the neighboring mesh draws in full.  Leave those out.  An
! obstruction the user made zero cells thick is flagged THIN and is still drawn.

IF (.NOT.OB%THIN) THEN
   SELECT CASE(ABS(IOR))
      CASE(1) ; IF (OB%I1==OB%I2) RETURN
      CASE(2) ; IF (OB%J1==OB%J2) RETURN
      CASE(3) ; IF (OB%K1==OB%K2) RETURN
   END SELECT
ENDIF

! Skip faces of zero area, which span no cells

SELECT CASE(ABS(IOR))
   CASE(1) ; IF (OB%J2<=OB%J1 .OR. OB%K2<=OB%K1) RETURN
   CASE(2) ; IF (OB%I2<=OB%I1 .OR. OB%K2<=OB%K1) RETURN
   CASE(3) ; IF (OB%I2<=OB%I1 .OR. OB%J2<=OB%J1) RETURN
END SELECT

VTK_OBST_FACE_DRAWN = .TRUE.

END FUNCTION VTK_OBST_FACE_DRAWN


!> \brief Describe one face of an obstruction as a PATCH
!>
!> BUILD_VTK_SOLID_PHASE_GEOMETRY works from a PATCH, so an obstruction face is handed
!> to it as one.  The node indices are collapsed onto the plane of the face and the gas
!> cell ranges span the two directions the face extends in.
!>
!> \param OB Obstruction
!> \param IOR Orientation index of the face
!> \param PA PATCH descriptor to fill

SUBROUTINE SET_VTK_OBST_FACE_PATCH(OB,IOR,PA)

TYPE(OBSTRUCTION_TYPE), INTENT(IN) :: OB
INTEGER, INTENT(IN) :: IOR
TYPE(PATCH_TYPE), INTENT(OUT) :: PA

PA%IOR = IOR
PA%I1 = OB%I1 ; PA%I2 = OB%I2 ; PA%IG1 = OB%I1+1 ; PA%IG2 = OB%I2
PA%J1 = OB%J1 ; PA%J2 = OB%J2 ; PA%JG1 = OB%J1+1 ; PA%JG2 = OB%J2
PA%K1 = OB%K1 ; PA%K2 = OB%K2 ; PA%KG1 = OB%K1+1 ; PA%KG2 = OB%K2

SELECT CASE(IOR)
   CASE(-1) ; PA%I2 = OB%I1
   CASE( 1) ; PA%I1 = OB%I2
   CASE(-2) ; PA%J2 = OB%J1
   CASE( 2) ; PA%J1 = OB%J2
   CASE(-3) ; PA%K2 = OB%K1
   CASE( 3) ; PA%K1 = OB%K2
END SELECT

END SUBROUTINE SET_VTK_OBST_FACE_PATCH



!> \brief Convert an FDS triangulated surface into VTK triangles


!>


!> \param VERTS Vertex coordinates, three per vertex


!> \param FACES Vertex indices, three per face


!> \param NCELLS Number of cells (out)


!> \param NPOINTS Number of points (out)


!> \param X_PTS Point x coordinates, allocated here (out)


!> \param Y_PTS Point y coordinates, allocated here (out)


!> \param Z_PTS Point z coordinates, allocated here (out)


!> \param CONNECT Point indices of each cell, allocated here (out)


!> \param OFFSETS Where each cell starts in CONNECT, allocated here (out)


!> \param VTKC_TYPE VTK cell type of each cell, allocated here (out)



SUBROUTINE BUILD_VTK_GEOM_GEOMETRY(VERTS, FACES, NCELLS, NPOINTS,&
                                   X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)

REAL(FB), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: VERTS
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: FACES
REAL(FB), ALLOCATABLE, DIMENSION(:) :: X_PTS, Y_PTS, Z_PTS
INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS
INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: VTKC_TYPE
INTEGER :: I
INTEGER :: NPOINTS, NCELLS, NFACES

NPOINTS = SIZE(VERTS)/3
NCELLS = SIZE(FACES)/3
NFACES = SIZE(FACES)

ALLOCATE(X_PTS(NPOINTS))
ALLOCATE(Y_PTS(NPOINTS))
ALLOCATE(Z_PTS(NPOINTS))
ALLOCATE(CONNECT(NCELLS*3))
ALLOCATE(OFFSETS(NCELLS))
ALLOCATE(VTKC_TYPE(NCELLS))

DO I=1,NPOINTS
   X_PTS(I) = VERTS((I-1)*3+1)
   Y_PTS(I) = VERTS((I-1)*3+2)
   Z_PTS(I) = VERTS((I-1)*3+3)
ENDDO

DO I=1,NCELLS
   OFFSETS(I) = (I)*3_IB32
   VTKC_TYPE(I) = 5_IB8
   CONNECT(3*(I-1)+1) = FACES(3*(I-1)+1)-1
   CONNECT(3*(I-1)+2) = FACES(3*(I-1)+2)-1
   CONNECT(3*(I-1)+3) = FACES(3*(I-1)+3)-1
ENDDO

ENDSUBROUTINE BUILD_VTK_GEOM_GEOMETRY




!> \brief Release the arrays returned by the geometry builders




!>




!> \param X_PTS Point x coordinates




!> \param Y_PTS Point y coordinates




!> \param Z_PTS Point z coordinates




!> \param OFFSETS Where each cell starts in CONNECT




!> \param VTKC_TYPE VTK cell type of each cell




!> \param CONNECT Point indices of each cell





SUBROUTINE DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)

REAL(FB), DIMENSION(:), ALLOCATABLE :: X_PTS, Y_PTS, Z_PTS
INTEGER(IB32), DIMENSION(:), ALLOCATABLE :: CONNECT, OFFSETS
INTEGER(IB8), DIMENSION(:), ALLOCATABLE :: VTKC_TYPE

DEALLOCATE(VTKC_TYPE)
DEALLOCATE(OFFSETS)
DEALLOCATE(CONNECT)
DEALLOCATE(X_PTS)
DEALLOCATE(Y_PTS)
DEALLOCATE(Z_PTS)

ENDSUBROUTINE DEALLOCATE_VTK_GAS_PHASE_GEOMETRY







!> \brief Build the file-access property list used for every VTKHDF file
!>
!> Collective metadata I/O is the important part.  By default every rank independently
!> reads the same superblock, group and object-header blocks when a file or dataset is
!> opened, and independently writes the same metadata when one is created, which turns
!> into a read/write storm on the same file offsets as the rank count grows.  With
!> collective metadata operations rank 0 does the I/O and broadcasts.  The larger metadata
!> block size trades a little file space for far fewer, larger metadata transfers.

SUBROUTINE GET_VTK_FAPL(PLIST_ID)

INTEGER(HID_T), INTENT(OUT) :: PLIST_ID
INTEGER :: ERROR

CALL H5PCREATE_F(H5P_FILE_ACCESS_F, PLIST_ID, ERROR)

! On a single process there is nothing to coordinate, and the MPI-IO driver is not free:
! it turns every flush into an MPI_File_sync and routes every write through its collective
! machinery.  HDF5's default driver does the same work without any of that.

IF (N_MPI_PROCESSES>1) THEN
   CALL H5PSET_FAPL_MPIO_F(PLIST_ID, MPI_COMM_WORLD, MPI_INFO_NULL, ERROR)
   CALL H5PSET_ALL_COLL_METADATA_OPS_F(PLIST_ID, .TRUE., ERROR)
   CALL H5PSET_COLL_METADATA_WRITE_F(PLIST_ID, .TRUE., ERROR)
ENDIF
CALL H5PSET_META_BLOCK_SIZE_F(PLIST_ID, VTK_META_BLOCK_SIZE, ERROR)

END SUBROUTINE GET_VTK_FAPL


!> \brief Create a single time VTKHDF file and open its groups


!>


!> \param FILENAME Name of the file to create


!> \param FILE_ID HDF5 file identifier (out)


!> \param PLIST_ID Data transfer property list shared by writes to this file (out)


!> \param GROUP_ID1 VTKHDF group (out)


!> \param GROUP_ID2 VTKHDF/CellData group (out)


!> \param GROUP_ID3 VTKHDF/FieldData group (out)


!> \param GROUP_ID4 VTKHDF/PointData group (out)


!>


!> Used for output that does not vary over time, such as the geometry file.



SUBROUTINE CREATE_OPEN_VTKHDF(FILENAME,&
                       FILE_ID, PLIST_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4)
   CHARACTER(*), INTENT(IN) :: FILENAME
   INTEGER(HID_T), INTENT(OUT) :: FILE_ID, PLIST_ID       ! Identifiers
   INTEGER(HID_T), INTENT(OUT) :: GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4 ! Group identifier
   INTEGER  :: ERROR                    !< IO Error status.
   INTEGER(KIND=MPI_INTEGER_KIND) :: MPI_SIZE, MPI_RANK, MPIERROR
   INTEGER(HSIZE_T) :: VER(2)
   INTEGER(HSIZE_T), DIMENSION(1) :: DATA_DIMS    ! Attribute dimension
   
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_SIZE, MPIERROR)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, MPIERROR)
   
   ! Setup file access property list with parallel I/O access.
   CALL GET_VTK_FAPL(PLIST_ID)
   
   ! Create the file collectively.
   CALL H5FCREATE_F(FILENAME, H5F_ACC_TRUNC_F, FILE_ID, ERROR, ACCESS_PRP = PLIST_ID)
   CALL H5PCLOSE_F(PLIST_ID, ERROR)
   
   ! Create VTKHDF groups
   CALL H5GCREATE_F(FILE_ID, "VTKHDF", GROUP_ID1, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/CellData", GROUP_ID2, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/FieldData", GROUP_ID3, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/PointData", GROUP_ID4, ERROR)
   
   ! Add VTKHDF header information
   VER = 2_HSIZE_T
   DATA_DIMS = 1_HSIZE_T
   CALL ADD_VERSION(GROUP_ID1,VER(1),1,"Version",VER)
   CALL ADD_ATTRIBUTE_CHAR(GROUP_ID1,DATA_DIMS,"Type","UnstructuredGrid",16_SIZE_T)

END SUBROUTINE CREATE_OPEN_VTKHDF


!> \brief Create a time series VTKHDF file and open its groups


!>


!> \param FILENAME Name of the file to create


!> \param FILE_ID HDF5 file identifier (out)


!> \param PLIST_ID Data transfer property list shared by writes to this file (out)


!> \param G1 VTKHDF group (out)


!> \param G2 VTKHDF/CellData group (out)


!> \param G3 VTKHDF/FieldData group (out)


!> \param G4 VTKHDF/PointData group (out)


!> \param G5 VTKHDF/Steps group (out)


!> \param G6 VTKHDF/Steps/CellDataOffsets group (out)


!> \param G7 VTKHDF/Steps/PointDataOffsets group (out)


!>


!> The Steps groups are what make the file temporal: each output time appends one


!> entry to the arrays under them saying where that step's data begins.



SUBROUTINE CREATE_OPEN_VTKHDF_SERIES(FILENAME,FILE_ID,PLIST_ID,G1,G2,G3,G4,G5,G6,G7)
   CHARACTER(*), INTENT(IN) :: FILENAME
   INTEGER(HID_T), INTENT(OUT) :: FILE_ID        ! Identifiers
   INTEGER(HID_T), INTENT(INOUT) :: PLIST_ID       ! Identifiers
   INTEGER(HID_T), INTENT(OUT) :: G1,G2,G3,G4,G5,G6,G7 ! Group identifier
   INTEGER  :: ERROR                    !< IO Error status.
   INTEGER(KIND=MPI_INTEGER_KIND) :: MPI_SIZE, MPI_RANK, MPIERROR
   INTEGER(HSIZE_T) :: VER(2)
   INTEGER(HSIZE_T), DIMENSION(1) :: DATA_DIMS    ! Attribute dimension
   
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_SIZE, MPIERROR)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, MPIERROR)
   
   ! Setup file access property list with parallel I/O access.
   CALL GET_VTK_FAPL(PLIST_ID)
   
   ! Create the file collectively.
   CALL H5FCREATE_F(FILENAME, H5F_ACC_TRUNC_F, FILE_ID, ERROR, ACCESS_PRP = PLIST_ID)
   CALL H5PCLOSE_F(PLIST_ID, ERROR)
   
   ! Create VTKHDF groups
   CALL H5GCREATE_F(FILE_ID, "VTKHDF", G1, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/CellData", G2, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/FieldData", G3, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/PointData", G4, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/Steps", G5, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/Steps/CellDataOffsets", G6, ERROR)
   CALL H5GCREATE_F(FILE_ID, "VTKHDF/Steps/PointDataOffsets", G7, ERROR)
   
   ! Add VTKHDF header information
   VER = 2_HSIZE_T
   DATA_DIMS = 1_HSIZE_T
   CALL ADD_VERSION(G1,VER(1),1,"Version",VER)
   CALL ADD_ATTRIBUTE_CHAR(G1,DATA_DIMS,"Type","UnstructuredGrid",16_SIZE_T)
   
   DATA_DIMS = 1_HSIZE_T
   CALL ADD_ATTRIBUTE_INT(G5,DATA_DIMS,"NSteps",INT(NFRAMES,HID_T))

END SUBROUTINE CREATE_OPEN_VTKHDF_SERIES


!> \brief Reopen an existing time series VTKHDF file and its groups


!>


!> \param FILENAME Name of the file to open


!> \param FILE_ID HDF5 file identifier (out)


!> \param PLIST_ID Data transfer property list shared by writes to this file (out)


!> \param G1 VTKHDF group (out)


!> \param G2 VTKHDF/CellData group (out)


!> \param G3 VTKHDF/FieldData group (out)


!> \param G4 VTKHDF/PointData group (out)


!> \param G5 VTKHDF/Steps group (out)


!> \param G6 VTKHDF/Steps/CellDataOffsets group (out)


!> \param G7 VTKHDF/Steps/PointDataOffsets group (out)


!>


!> Only reached when VTK_KEEPOPEN is false, which reopens each file at every


!> output time rather than holding it open for the run.



SUBROUTINE OPEN_VTKHDF_SERIES(FILENAME,FILE_ID,PLIST_ID,G1,G2,G3,G4,G5,G6,G7)
   CHARACTER(*), INTENT(IN) :: FILENAME
   INTEGER(HID_T), INTENT(OUT) :: FILE_ID       ! Identifiers
   INTEGER(HID_T), INTENT(INOUT) :: PLIST_ID       ! Identifiers
   INTEGER(HID_T), INTENT(OUT) :: G1,G2,G3,G4,G5,G6,G7 ! Group identifier
   INTEGER  :: ERROR                    !< IO Error status.
   INTEGER(KIND=MPI_INTEGER_KIND) :: MPI_SIZE, MPI_RANK, MPIERROR
   
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_SIZE, MPIERROR)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, MPIERROR)
   
   ! Setup file access property list with parallel I/O access.
   CALL GET_VTK_FAPL(PLIST_ID)
   
   ! Create the file collectively.
   CALL H5FOPEN_F(FILENAME, H5F_ACC_RDWR_F, FILE_ID, ERROR, ACCESS_PRP = PLIST_ID)
   CALL H5PCLOSE_F(PLIST_ID, ERROR)
   
   ! Create VTKHDF groups
   CALL H5GOPEN_F(FILE_ID, "VTKHDF", G1, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/CellData", G2, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/FieldData", G3, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/PointData", G4, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/Steps", G5, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/Steps/CellDataOffsets", G6, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/Steps/PointDataOffsets", G7, ERROR)
   
END SUBROUTINE OPEN_VTKHDF_SERIES



!> \brief Reopen an existing single time VTKHDF file and its groups



!>



!> \param FILENAME Name of the file to open



!> \param FILE_ID HDF5 file identifier (out)



!> \param PLIST_ID Data transfer property list shared by writes to this file (out)



!> \param GROUP_ID1 VTKHDF group (out)



!> \param GROUP_ID2 VTKHDF/CellData group (out)



!> \param GROUP_ID3 VTKHDF/FieldData group (out)



!> \param GROUP_ID4 VTKHDF/PointData group (out)




SUBROUTINE OPEN_VTKHDF(FILENAME,&
                       FILE_ID, PLIST_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4)
   CHARACTER(*), INTENT(IN) :: FILENAME
   INTEGER(HID_T), INTENT(OUT) :: FILE_ID, PLIST_ID       ! Identifiers
   INTEGER(HID_T), INTENT(OUT) :: GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4 ! Group identifier
   INTEGER  :: ERROR                    !< IO Error status.
   INTEGER(KIND=MPI_INTEGER_KIND) :: MPI_SIZE, MPI_RANK, MPIERROR
   
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_SIZE, MPIERROR)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, MPIERROR)
   
   ! Setup file access property list with parallel I/O access.
   CALL GET_VTK_FAPL(PLIST_ID)
   
   ! Create the file collectively.
   CALL H5FOPEN_F(FILENAME, H5F_ACC_RDWR_F, FILE_ID, ERROR, ACCESS_PRP = PLIST_ID)
   CALL H5PCLOSE_F(PLIST_ID, ERROR)
   
   ! Create VTKHDF groups
   CALL H5GOPEN_F(FILE_ID, "VTKHDF", GROUP_ID1, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/CellData", GROUP_ID2, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/FieldData", GROUP_ID3, ERROR)
   CALL H5GOPEN_F(FILE_ID, "VTKHDF/PointData", GROUP_ID4, ERROR)

END SUBROUTINE OPEN_VTKHDF


!> \brief Close a single time VTKHDF file and its groups


!>


!> \param FILE_ID HDF5 file identifier


!> \param GROUP_ID1 VTKHDF group


!> \param GROUP_ID2 VTKHDF/CellData group


!> \param GROUP_ID3 VTKHDF/FieldData group


!> \param GROUP_ID4 VTKHDF/PointData group



SUBROUTINE CLOSE_VTKHDF(FILE_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4)
   INTEGER(HID_T), INTENT(IN) :: FILE_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4
   INTEGER  :: ERROR                    !< IO Error status.
   CALL VTK_DSET_CACHE_PURGE()
   CALL H5GCLOSE_F(GROUP_ID4, ERROR)
   CALL H5GCLOSE_F(GROUP_ID3, ERROR)
   CALL H5GCLOSE_F(GROUP_ID2, ERROR)
   CALL H5GCLOSE_F(GROUP_ID1, ERROR)
   CALL H5FCLOSE_F(FILE_ID, ERROR)
END SUBROUTINE CLOSE_VTKHDF

!> \brief Close a time series VTKHDF file and its groups

!>

!> \param FILE_ID HDF5 file identifier

!> \param G1 VTKHDF group

!> \param G2 VTKHDF/CellData group

!> \param G3 VTKHDF/FieldData group

!> \param G4 VTKHDF/PointData group

!> \param G5 VTKHDF/Steps group

!> \param G6 VTKHDF/Steps/CellDataOffsets group

!> \param G7 VTKHDF/Steps/PointDataOffsets group

!>

!> The cached dataset handles are closed first; HDF5 will not write a file out

!> completely while objects inside it are still open.


SUBROUTINE CLOSE_VTKHDF_SERIES(FILE_ID,G1,G2,G3,G4,G5,G6,G7)
   INTEGER(HID_T), INTENT(IN) :: FILE_ID, G1,G2,G3,G4,G5,G6,G7
   INTEGER  :: ERROR                    !< IO Error status.
   CALL VTK_DSET_CACHE_PURGE()
   CALL H5FFLUSH_F(FILE_ID,H5F_SCOPE_LOCAL_F,ERROR)
   CALL H5GCLOSE_F(G7, ERROR)
   CALL H5GCLOSE_F(G6, ERROR)
   CALL H5GCLOSE_F(G5, ERROR)
   CALL H5GCLOSE_F(G4, ERROR)
   CALL H5GCLOSE_F(G3, ERROR)
   CALL H5GCLOSE_F(G2, ERROR)
   CALL H5GCLOSE_F(G1, ERROR)
   CALL H5FCLOSE_F(FILE_ID, ERROR)
END SUBROUTINE CLOSE_VTKHDF_SERIES

!> \brief Lengthen a slice quantity's dataset to make room for one more output time

!>

!> \param DATANAME Name of the quantity being written

!> \param II Index of the unique slice plane, which selects the file

!> \param BASE_OFFSET Where this output time's data begins in the dataset (out)

!> \param NCELLS Number of cells each mesh contributes (out)

!> \param NPOINTS Number of points each mesh contributes (out)

!>

!> Every rank needs the same picture of how much each mesh contributes, so the

!> per mesh counts are exchanged with MPI_ALLREDUCE rather than gathered on one

!> rank and scattered back.


SUBROUTINE EXTEND_SLICE_VTKHDF(DATANAME, II, BASE_OFFSET, NCELLS, NPOINTS)
   CHARACTER(*), INTENT(IN) :: DATANAME
   INTEGER, INTENT(IN) :: II
   !INTEGER(HID_T) :: CRP_LIST       ! Identifiers
   INTEGER(HID_T) :: DSET_ID    ! Dataset identifiers
   INTEGER :: NPOINTS_TOTAL, NCELLS_TOTAL
   INTEGER(IB32), DIMENSION(1:NMESHES), INTENT(OUT) :: NPOINTS, NCELLS
   INTEGER :: NCELLS_MAX, NPOINTS_MAX
   INTEGER :: NM,ERROR
   INTEGER(HSIZE_T), DIMENSION(1) :: NPOINTS_MAX_ARRAY(1), NPOINTS_TOTAL_ARRAY(1)
   INTEGER(HSIZE_T), DIMENSION(1) :: SIZE1, MDIM, CDIM, DDIM
   INTEGER(HSIZE_T), INTENT(OUT) :: BASE_OFFSET
   INTEGER(HID_T) :: DATASPACE
   !INTEGER(HSIZE_T), DIMENSION(1) :: START1(1)

   ! Read number of cells
   !CALL H5DOPEN_F(HDF_SLCF_G1(II), "NumberOfCells", DSET_ID, ERROR)
   !START1 = NMESHES
   !CALL H5DREAD_F(DSET_ID, H5T_STD_I32LE, NCELLS, START1, ERROR)
   !CALL VTK_DSET_RELEASE(DSET_ID, ERROR)
   NCELLS = HDF_SLCF_G1_NCELLS(II,1:NMESHES)

   ! Read number of points
   !CALL H5DOPEN_F(HDF_SLCF_G1(II), "NumberOfPoints", DSET_ID, ERROR)
   !CALL H5DREAD_F(DSET_ID, H5T_STD_I32LE, NPOINTS, START1, ERROR)
   !CALL VTK_DSET_RELEASE(DSET_ID, ERROR)
   NPOINTS = HDF_SLCF_G1_NPOINTS(II,1:NMESHES)

   NPOINTS_TOTAL=0
   NCELLS_TOTAL=0
   NCELLS_MAX = 0
   NPOINTS_MAX = 0
   
   MESH_LOOP_HDF_COUNT: DO NM=1,NMESHES
      NPOINTS_TOTAL = NPOINTS_TOTAL + NPOINTS(NM) !+ NVERTS*3
      NCELLS_TOTAL = NCELLS_TOTAL + NCELLS(NM) !+ NFACES*3
      NCELLS_MAX = MAX(NCELLS_MAX,NCELLS(NM))
      NPOINTS_MAX = MAX(NPOINTS_MAX,NPOINTS(NM))
   ENDDO MESH_LOOP_HDF_COUNT
   
   NPOINTS_MAX_ARRAY(1) = INT(NPOINTS_MAX, HSIZE_T)
   NPOINTS_TOTAL_ARRAY(1) = INT(NPOINTS_TOTAL, HSIZE_T)
   MDIM = (/H5S_UNLIMITED_F/)
   DDIM = (/0_HSIZE_T/)
   CALL PARALLEL_INIT_F32(HDF_SLCF_G4(II), TRIM(DATANAME), HDF_SLCF_CRP_LIST(II), 1,&
      NPOINTS_MAX_ARRAY, DDIM, DSET_ID, HDF_SLCF_PLIST_ID(II), MDIM)
   
   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE, CDIM, MDIM, ERROR)
   SIZE1 = CDIM(1)+NPOINTS_TOTAL_ARRAY(1)
   BASE_OFFSET = CDIM(1)
   CALL H5DSET_EXTENT_F(DSET_ID, SIZE1, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID, ERROR)
   
END SUBROUTINE EXTEND_SLICE_VTKHDF


!> \brief Write every staged boundary-data block in a single collective operation

SUBROUTINE WRITE_VTKHDF_BNDF_DATA_MULTI(DATASET,GROUP_ID,N_BLK,BLK_START,BLK_COUNT,DATA)
   CHARACTER(*), INTENT(IN) :: DATASET
   INTEGER(HID_T), INTENT(IN) :: GROUP_ID
   INTEGER, INTENT(IN) :: N_BLK
   INTEGER, DIMENSION(:), INTENT(IN) :: BLK_START, BLK_COUNT
   REAL(FB), DIMENSION(*), INTENT(IN) :: DATA
   INTEGER(HID_T) :: DSET_ID, MEMSPACE, DATASPACE
   INTEGER(HSIZE_T), DIMENSION(1) :: START1, EXTENT1
   INTEGER(HSIZE_T) :: NLOC
   INTEGER :: I, ERROR

   START1 = 0_HSIZE_T
   EXTENT1 = 0_HSIZE_T
   CALL PARALLEL_INIT_F32(GROUP_ID, DATASET, HDF_BNDF_CRP_LIST, 1,&
      START1, EXTENT1, DSET_ID, HDF_BNDF_PLIST_ID)

   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SSELECT_NONE_F(DATASPACE, ERROR)

   NLOC = 0_HSIZE_T
   DO I=1,N_BLK
      IF (BLK_COUNT(I)<=0) CYCLE
      START1  = INT(BLK_START(I),HSIZE_T)
      EXTENT1 = INT(BLK_COUNT(I),HSIZE_T)
      CALL H5SSELECT_HYPERSLAB_F(DATASPACE, H5S_SELECT_OR_F, START1, EXTENT1, ERROR)
      NLOC = NLOC + EXTENT1(1)
   ENDDO

   EXTENT1 = NLOC
   CALL H5SCREATE_SIMPLE_F(1, EXTENT1, MEMSPACE, ERROR)
   IF (NLOC==0_HSIZE_T) CALL H5SSELECT_NONE_F(MEMSPACE, ERROR)
   CALL H5DWRITE_F(DSET_ID, H5T_IEEE_F32LE, DATA, EXTENT1, ERROR, &
                   MEM_SPACE_ID = MEMSPACE, FILE_SPACE_ID = DATASPACE, XFER_PRP = HDF_BNDF_PLIST_ID)

   CALL H5SCLOSE_F(MEMSPACE, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID, ERROR)

END SUBROUTINE WRITE_VTKHDF_BNDF_DATA_MULTI


!> \brief Lengthen the boundary file's datasets to make room for one more output time


!>


!> \param BASE_OFFSET_CELLS Where this output time's cell data begins (out)


!> \param BASE_OFFSET_PTS Where this output time's point data begins (out)



SUBROUTINE EXTEND_BNDF_VTKHDF(BASE_OFFSET_CELLS, BASE_OFFSET_PTS)
   TYPE (BOUNDARY_FILE_TYPE), POINTER :: BF
   INTEGER(HID_T) :: DSET_ID    ! Dataset identifiers
   INTEGER :: NF,ERROR
   INTEGER(HSIZE_T), DIMENSION(1) :: SIZE1, MDIM, CDIM
   INTEGER(HSIZE_T), INTENT(OUT) :: BASE_OFFSET_CELLS, BASE_OFFSET_PTS
   INTEGER(HID_T) :: DATASPACE
   INTEGER(HSIZE_T), DIMENSION(1) :: START1(1), EXTENT1(1)

EXTEND_LOOP: DO NF=1,N_BNDF
   BF => BOUNDARY_FILE(NF)
   MDIM = (/H5S_UNLIMITED_F/)
   START1=0
   EXTENT1=0
   IF (BF%CELL_CENTERED) THEN
      CALL PARALLEL_INIT_F32(HDF_BNDF_G2, BF%SMOKEVIEW_LABEL(1:30), HDF_BNDF_CRP_LIST, 1,&
         START1, EXTENT1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM) ! Open cell centered to read
   ELSE
      CALL PARALLEL_INIT_F32(HDF_BNDF_G4, BF%SMOKEVIEW_LABEL(1:30), HDF_BNDF_CRP_LIST, 1,&
         START1, EXTENT1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM) ! Open nodal to read
   ENDIF
   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE, CDIM, MDIM, ERROR)
   IF (BF%CELL_CENTERED) THEN
      SIZE1 = CDIM(1)+SUM(NCELLS_VTK)
      BASE_OFFSET_CELLS = CDIM(1)
   ELSE
      SIZE1 = CDIM(1)+SUM(NPOINTS_VTK)
      BASE_OFFSET_PTS = CDIM(1)
   ENDIF
   CALL H5DSET_EXTENT_F(DSET_ID, SIZE1, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID, ERROR)
ENDDO EXTEND_LOOP
   
END SUBROUTINE EXTEND_BNDF_VTKHDF



!> \brief Lengthen a Smoke3D quantity's dataset to make room for one more output time



!>



!> \param DATANAME Name of the quantity being written



!> \param BASE_OFFSET Where this output time's data begins in the dataset (out)




SUBROUTINE EXTEND_SMOKE3D_VTKHDF(DATANAME, BASE_OFFSET)
   CHARACTER(*), INTENT(IN) :: DATANAME
   !INTEGER(HID_T) :: CRP_LIST       ! Identifiers
   INTEGER(HID_T) :: DSET_ID    ! Dataset identifiers
   INTEGER :: NPOINTS_TOTAL, NCELLS_TOTAL, NOFFSETS_TOTAL
   INTEGER :: NPOINTS, NCELLS
   INTEGER :: NCONN_TOTAL, NPIECES
   INTEGER :: NOFFSETS_ACCUM, NCONN_ACCUM, NPIECES_ACCUM, NCELLS_ACCUM, NPOINTS_ACCUM !, NFACES
   INTEGER :: NCONN_MAX, NOFFSETS_MAX, NCELLS_MAX, NPOINTS_MAX
   INTEGER :: NM,ERROR
   INTEGER, DIMENSION(1:N_MPI_PROCESSES) :: MESHES_PER_PROCESS
   !INTEGER :: N_WRITTEN
   INTEGER(HSIZE_T), DIMENSION(1) :: NPOINTS_MAX_ARRAY(1), NPOINTS_TOTAL_ARRAY(1)
   INTEGER(HSIZE_T), DIMENSION(1) :: SIZE1, MDIM, CDIM, DDIM
   INTEGER(HSIZE_T), INTENT(OUT) :: BASE_OFFSET
   INTEGER(HID_T) :: DATASPACE


   DO NM=1,N_MPI_PROCESSES
      MESHES_PER_PROCESS(NM) = 0
   ENDDO
   
   NPOINTS_TOTAL=0
   NCELLS_TOTAL=0
   NCELLS_ACCUM=0
   NPOINTS_ACCUM=0
   NOFFSETS_TOTAL=0
   NOFFSETS_ACCUM = 0
   NCONN_TOTAL = 0
   NPIECES = 0
   NCONN_ACCUM = 0
   NPIECES_ACCUM = 0
   NCELLS_MAX = 0
   NPOINTS_MAX = 0
   NCONN_MAX = 0
   NOFFSETS_MAX = 0
   
   MESH_LOOP_HDF_COUNT: DO NM=1,NMESHES
      NCELLS = MESHES(NM)%NC
      NPOINTS = MESHES(NM)%NP
      NPOINTS_TOTAL = NPOINTS_TOTAL + NPOINTS !+ NVERTS*3
      NCELLS_TOTAL = NCELLS_TOTAL + NCELLS !+ NFACES*3
      NOFFSETS_TOTAL = NOFFSETS_TOTAL + NCELLS + 1
      NCONN_TOTAL = NCONN_TOTAL + NCELLS*8
      NPIECES = NPIECES + 1
      NCELLS_MAX = MAX(NCELLS_MAX,NCELLS)
      NPOINTS_MAX = MAX(NPOINTS_MAX,NPOINTS)
      NOFFSETS_MAX = MAX(NOFFSETS_MAX,NCELLS + 1)
      NCONN_MAX = MAX(NCONN_MAX,NCELLS*8)
      MESHES_PER_PROCESS(PROCESS(NM)+1) = MESHES_PER_PROCESS(PROCESS(NM)+1) + 1
   ENDDO MESH_LOOP_HDF_COUNT
   
   ! Data Quantity
   NPOINTS_MAX_ARRAY(1) = INT(NPOINTS_MAX, HSIZE_T)
   NPOINTS_TOTAL_ARRAY(1) = INT(NPOINTS_TOTAL, HSIZE_T)
   MDIM = (/H5S_UNLIMITED_F/)
   DDIM = (/0_HSIZE_T/)
   CALL PARALLEL_INIT_U8(HDF_SM3D_G4, TRIM(DATANAME), HDF_SM3D_CRP_LIST, 1,&
      NPOINTS_MAX_ARRAY, DDIM, DSET_ID, HDF_SM3D_PLIST_ID, MDIM) ! Values = Timesteps
   
   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE, CDIM, MDIM, ERROR)
   SIZE1 = CDIM(1)+NPOINTS_TOTAL_ARRAY(1)
   BASE_OFFSET = CDIM(1)
   CALL H5DSET_EXTENT_F(DSET_ID, SIZE1, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID, ERROR)
   
END SUBROUTINE EXTEND_SMOKE3D_VTKHDF



!> \brief Write every locally owned mesh's Smoke3D data in a single collective operation

SUBROUTINE ADD_DATA_TO_SMOKE3D_VTKHDF_MULTI(DATANAME,DATA,BASE_OFFSET,N_LOCAL,WRITE_NM)
   CHARACTER(*), INTENT(IN) :: DATANAME
   INTEGER(IB8), DIMENSION(*), INTENT(IN) :: DATA
   INTEGER(HSIZE_T), INTENT(IN) :: BASE_OFFSET
   INTEGER, INTENT(IN) :: N_LOCAL
   INTEGER, DIMENSION(:), INTENT(IN) :: WRITE_NM
   INTEGER(HID_T) :: DSET_ID, MEMSPACE, DATASPACE
   INTEGER(HSIZE_T) :: OFFSET_ACCUM(NMESHES), NPOINTS_MAX, NPOINTS_TOTAL, NLOC
   INTEGER(HSIZE_T), DIMENSION(1) :: START1, EXTENT1, MDIM
   INTEGER :: NM, I, ERROR

   NPOINTS_TOTAL = 0_HSIZE_T
   NPOINTS_MAX = 0_HSIZE_T
   DO NM=1,NMESHES
      OFFSET_ACCUM(NM) = NPOINTS_TOTAL
      NPOINTS_TOTAL = NPOINTS_TOTAL + INT(MESHES(NM)%NP,HSIZE_T)
      NPOINTS_MAX = MAX(NPOINTS_MAX,INT(MESHES(NM)%NP,HSIZE_T))
   ENDDO

   MDIM = (/H5S_UNLIMITED_F/)
   START1 = NPOINTS_MAX
   EXTENT1 = 0_HSIZE_T
   CALL PARALLEL_INIT_U8(HDF_SM3D_G4, TRIM(DATANAME), HDF_SM3D_CRP_LIST, 1,&
      START1, EXTENT1, DSET_ID, HDF_SM3D_PLIST_ID, MDIM)

   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SSELECT_NONE_F(DATASPACE, ERROR)

   NLOC = 0_HSIZE_T
   DO I=1,N_LOCAL
      NM = WRITE_NM(I)
      START1  = OFFSET_ACCUM(NM) + BASE_OFFSET
      EXTENT1 = INT(MESHES(NM)%NP,HSIZE_T)
      IF (EXTENT1(1)==0_HSIZE_T) CYCLE
      CALL H5SSELECT_HYPERSLAB_F(DATASPACE, H5S_SELECT_OR_F, START1, EXTENT1, ERROR)
      NLOC = NLOC + EXTENT1(1)
   ENDDO

   EXTENT1 = NLOC
   CALL H5SCREATE_SIMPLE_F(1, EXTENT1, MEMSPACE, ERROR)
   IF (NLOC==0_HSIZE_T) CALL H5SSELECT_NONE_F(MEMSPACE, ERROR)
   CALL H5DWRITE_F(DSET_ID, H5T_STD_U8LE, DATA, EXTENT1, ERROR, &
                   MEM_SPACE_ID = MEMSPACE, FILE_SPACE_ID = DATASPACE, XFER_PRP = HDF_SM3D_PLIST_ID)

   CALL H5SCLOSE_F(MEMSPACE, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID, ERROR)

END SUBROUTINE ADD_DATA_TO_SMOKE3D_VTKHDF_MULTI






!> \brief Create the Smoke3D VTKHDF file and write the grid it will refer to






!>






!> The grid does not change over the run, so it is written once here and every






!> output time's Steps entry points back at it.







SUBROUTINE INITIALIZE_VTKHDF_SMOKE3D()
   CHARACTER(FN_LENGTH)   :: FILENAME
   INTEGER(HID_T) :: CRP_LIST       ! Identifiers
   INTEGER(HID_T) :: DSET_ID_CON, DSET_ID_NCELLS, DSET_ID_NCON    ! Dataset identifiers
   INTEGER(HID_T) :: DSET_ID_NPTS, DSET_ID_OFF, DSET_ID_PTS, DSET_ID_TYP       ! Dataset identifiers
   INTEGER :: NPOINTS_TOTAL, NCELLS_TOTAL, NOFFSETS_TOTAL
   INTEGER :: NPOINTS, NCELLS
   INTEGER :: NCONN_TOTAL, NPIECES
   INTEGER :: NOFFSETS_ACCUM, NCONN_ACCUM, NPIECES_ACCUM, NCELLS_ACCUM, NPOINTS_ACCUM !, NFACES
   INTEGER :: NCONN_MAX, NOFFSETS_MAX, NCELLS_MAX, NPOINTS_MAX
   INTEGER :: NM, ERROR
   TYPE (MESH_TYPE), POINTER :: M
   REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: VERTICES
   INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS
   INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: VTKC_TYPE
   INTEGER, DIMENSION(1:N_MPI_PROCESSES) :: MESHES_PER_PROCESS
   INTEGER :: N_WRITTEN
   INTEGER(HSIZE_T), DIMENSION(1) :: START1(1), EXTENT1(1), START2(2), EXTENT2(2)
   INTEGER(IB32), DIMENSION(1) :: IDATA_OUT1(1)
   REAL(FB), DIMENSION(1) :: RDATA_OUT3(3)
   WRITE(FILENAME,'(A,A,A)') "",TRIM(VTK_DIR)//TRIM(CHID),'_SM3D.vtkhdf'
   DO NM=1,N_MPI_PROCESSES
      MESHES_PER_PROCESS(NM) = 0
   ENDDO
   
   NPOINTS_TOTAL=0
   NCELLS_TOTAL=0
   NCELLS_ACCUM=0
   NPOINTS_ACCUM=0
   NOFFSETS_TOTAL=0
   NOFFSETS_ACCUM = 0
   NCONN_TOTAL = 0
   NPIECES = 0
   NCONN_ACCUM = 0
   NPIECES_ACCUM = 0
   NCELLS_MAX = 0
   NPOINTS_MAX = 0
   NCONN_MAX = 0
   NOFFSETS_MAX = 0
   N_WRITTEN=0
   
   MESH_LOOP_HDF_COUNT: DO NM=1,NMESHES
      NCELLS = MESHES(NM)%NC
      NPOINTS = MESHES(NM)%NP
      NPOINTS_TOTAL = NPOINTS_TOTAL + NPOINTS !+ NVERTS*3
      NCELLS_TOTAL = NCELLS_TOTAL + NCELLS !+ NFACES*3
      NOFFSETS_TOTAL = NOFFSETS_TOTAL + NCELLS + 1
      NCONN_TOTAL = NCONN_TOTAL + NCELLS*8
      NPIECES = NPIECES + 1
      NCELLS_MAX = MAX(NCELLS_MAX,NCELLS)
      NPOINTS_MAX = MAX(NPOINTS_MAX,NPOINTS)
      NOFFSETS_MAX = MAX(NOFFSETS_MAX,NCELLS + 1)
      NCONN_MAX = MAX(NCONN_MAX,NCELLS*8)
      MESHES_PER_PROCESS(PROCESS(NM)+1) = MESHES_PER_PROCESS(PROCESS(NM)+1) + 1
   ENDDO MESH_LOOP_HDF_COUNT
   
   ! Initialize HDF5 datafiles
   CALL CREATE_OPEN_VTKHDF_SERIES(FILENAME,HDF_SM3D_FILE_ID,HDF_SM3D_PLIST_ID,&
      HDF_SM3D_G1,HDF_SM3D_G2,HDF_SM3D_G3,HDF_SM3D_G4,HDF_SM3D_G5,HDF_SM3D_G6,HDF_SM3D_G7)
   
   START1 = NCONN_MAX
   EXTENT1 = NCONN_TOTAL
   CALL PARALLEL_INIT_I32(HDF_SM3D_G1, TRIM("Connectivity"), CRP_LIST, 1, START1, EXTENT1,&
      DSET_ID_CON, HDF_SM3D_PLIST_ID) ! Connectivity
   START1 = NPIECES
   EXTENT1 = NPIECES
   CALL PARALLEL_INIT_I32(HDF_SM3D_G1, TRIM("NumberOfCells"), CRP_LIST, 1, START1, EXTENT1,&
      DSET_ID_NCELLS, HDF_SM3D_PLIST_ID) ! NumberOfCells
   CALL PARALLEL_INIT_I32(HDF_SM3D_G1, TRIM("NumberOfPoints"), CRP_LIST, 1, START1, EXTENT1,&
      DSET_ID_NPTS, HDF_SM3D_PLIST_ID) ! NumberOfPoints
   CALL PARALLEL_INIT_I32(HDF_SM3D_G1, TRIM("NumberOfConnectivityIds"), CRP_LIST, 1, START1, EXTENT1,&
      DSET_ID_NCON, HDF_SM3D_PLIST_ID) ! NumberOfConnectivityIds
   START1 = NOFFSETS_MAX
   EXTENT1 = NOFFSETS_TOTAL
   CALL PARALLEL_INIT_I32(HDF_SM3D_G1, TRIM("Offsets"), CRP_LIST, 1, START1, EXTENT1,&
      DSET_ID_OFF, HDF_SM3D_PLIST_ID) ! Offsets
   START2(1:2) = INT((/3, NPOINTS_MAX/),HSIZE_T)
   EXTENT2(1:2) = INT((/3, NPOINTS_TOTAL/),HSIZE_T)
   CALL PARALLEL_INIT_F32(HDF_SM3D_G1, TRIM("Points"), CRP_LIST, 2, START2,EXTENT2,&
      DSET_ID_PTS, HDF_SM3D_PLIST_ID) ! Points
   START1 = NCELLS_MAX
   EXTENT1 = NCELLS_TOTAL
   CALL PARALLEL_INIT_U8(HDF_SM3D_G1, TRIM("Types"), CRP_LIST, 1, START1, EXTENT1,&
      DSET_ID_TYP, HDF_SM3D_PLIST_ID) ! Types

   MESH_LOOP_HDF: DO NM=1,NMESHES
      NCELLS = MESHES(NM)%NC
      NPOINTS = MESHES(NM)%NP
      IF (PROCESS(NM)/=MY_RANK) THEN
         NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS
         NCELLS_ACCUM = NCELLS_ACCUM + NCELLS
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS + 1
         NCONN_ACCUM = NCONN_ACCUM + NCELLS*8
         NPIECES_ACCUM = NPIECES_ACCUM+1
         CYCLE MESH_LOOP_HDF
      ENDIF
      N_WRITTEN = N_WRITTEN+1
      M => MESHES(NM)
      CALL BUILD_VTK_GAS_PHASE_GEOMETRY2(NM, NCELLS, NPOINTS, VERTICES, CONNECT, OFFSETS, VTKC_TYPE)
      
      ! Write connectivity data to file
      START1 = NCELLS_ACCUM*8
      EXTENT1 = NCELLS*8
      CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, HDF_SM3D_PLIST_ID, START1, EXTENT1, CONNECT)
      
      ! Write number of cells data to file
      START1 = NM-1
      EXTENT1 = 1
      IDATA_OUT1 = NCELLS
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, HDF_SM3D_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      ! Write number of points data to file
      IDATA_OUT1 = NPOINTS
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, HDF_SM3D_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      ! Write NumberOfConnectivityIds data to file
      IDATA_OUT1 = NCELLS*8
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, HDF_SM3D_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      ! Write offsets data to file
      START1 = NOFFSETS_ACCUM
      EXTENT1 = NCELLS+1
      CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, HDF_SM3D_PLIST_ID, START1, EXTENT1, OFFSETS)
      
      ! Write point data to file
      START2(1:2) = INT((/0, NPOINTS_ACCUM/),HSIZE_T)
      EXTENT2(1:2) = INT((/3, NPOINTS/),HSIZE_T)
      CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, HDF_SM3D_PLIST_ID, START2, EXTENT2, VERTICES)
      
      ! Write types data to file
      START1 = NCELLS_ACCUM
      EXTENT1 = NCELLS
      CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, HDF_SM3D_PLIST_ID, START1, EXTENT1, VTKC_TYPE)
      
      DEALLOCATE(VERTICES)
      DEALLOCATE(OFFSETS)
      DEALLOCATE(CONNECT)
      DEALLOCATE(VTKC_TYPE)
      NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS
      NCELLS_ACCUM = NCELLS_ACCUM + NCELLS
      NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS + 1
      NCONN_ACCUM = NCONN_ACCUM + NCELLS*8
      NPIECES_ACCUM = NPIECES_ACCUM+1
   ENDDO MESH_LOOP_HDF
   
   
   ! Write fake data in processes that have less meshes than max process
   DO NM=1,MAXVAL(MESHES_PER_PROCESS)
      IF (N_WRITTEN < MAXVAL(MESHES_PER_PROCESS)) THEN
         NCELLS = 0
         NPOINTS = 0
         ALLOCATE(CONNECT(NCELLS))
         ALLOCATE(OFFSETS(NCELLS))
         ALLOCATE(VERTICES(3,NCELLS))
         ALLOCATE(VTKC_TYPE(NCELLS))
         ! Write connectivity data to file
         START1 = NCELLS_ACCUM*8
         EXTENT1 = NCELLS*8
         CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, HDF_SM3D_PLIST_ID, START1, EXTENT1, CONNECT)
         
         ! Write number of cells data to file
         START1 = 0
         EXTENT1 = 0
         IDATA_OUT1 = 0
         CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, HDF_SM3D_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
         
         ! Write number of points data to file
         CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, HDF_SM3D_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
         
         ! Write NumberOfConnectivityIds data to file
         CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, HDF_SM3D_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
         
         ! Write offsets data to file
         CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, HDF_SM3D_PLIST_ID, START1, EXTENT1, OFFSETS)
         
         ! Write point data to file
         START2(1:2) = INT((/0, 0/),HSIZE_T)
         EXTENT2(1:2) = INT((/3, NPOINTS/),HSIZE_T)
         RDATA_OUT3 = 0._FB
         CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, HDF_SM3D_PLIST_ID, START2, EXTENT2, RDATA_OUT3)
         
         ! Write types data to file
         CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, HDF_SM3D_PLIST_ID, START1, EXTENT1, VTKC_TYPE)
         
         DEALLOCATE(VERTICES)
         DEALLOCATE(OFFSETS)
         DEALLOCATE(CONNECT)
         DEALLOCATE(VTKC_TYPE)
         N_WRITTEN = N_WRITTEN + 1
      ENDIF
   ENDDO
   
   ! Close HDF5 datasets
   CALL VTK_DSET_RELEASE(DSET_ID_CON, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NCELLS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NPTS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NCON, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_OFF, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_PTS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_TYP, ERROR)
   
END SUBROUTINE INITIALIZE_VTKHDF_SMOKE3D

!> \brief Flush every open VTKHDF file to disk
!>
!> With MPI-IO, HDF5 defers the superblock, group and object-header writes.  A run that is
!> killed without closing its files therefore leaves them with no readable header at all --
!> the raw data blocks are present but nothing can open the file.  Flushing pushes the
!> metadata out so that the file on disk stays loadable while the run continues.

!> \brief Make the VTKHDF files written at this output time readable on disk
!>
!> \param WROTE_SLCF Slice files were written since the last flush
!> \param WROTE_SM3D The Smoke3D file was written since the last flush
!> \param WROTE_BNDF The boundary file was written since the last flush
!> \param WROTE_PART Particle files were written since the last flush
!>
!> Each output has its own clock, and this routine runs on every time step, so on most
!> time steps there is nothing new in any of these files: storage flushed 211 times for
!> 11 output times, couch 617 times for 51.  A flush of an unchanged file is not free --
!> it still costs a metadata cache walk and an MPI_File_sync per file -- so only flush
!> what was written.

SUBROUTINE FLUSH_VTKHDF(WROTE_SLCF,WROTE_SM3D,WROTE_BNDF,WROTE_PART)

LOGICAL, INTENT(IN) :: WROTE_SLCF,WROTE_SM3D,WROTE_BNDF,WROTE_PART
INTEGER :: II,N,ERROR

IF (WROTE_SLCF .AND. ALLOCATED(HDF_SLCF_FILE_ID)) THEN
   DO II=1,MESHES(1)%N_UNIQUE_SLCF
      CALL H5FFLUSH_F(HDF_SLCF_FILE_ID(II),H5F_SCOPE_GLOBAL_F,ERROR)
   ENDDO
ENDIF
IF (WROTE_SM3D .AND. N_SMOKE3D>0) CALL H5FFLUSH_F(HDF_SM3D_FILE_ID,H5F_SCOPE_GLOBAL_F,ERROR)
IF (WROTE_BNDF) CALL H5FFLUSH_F(HDF_BNDF_FILE_ID,H5F_SCOPE_GLOBAL_F,ERROR)
IF (WROTE_PART .AND. ALLOCATED(HDF_PART_FILE_ID)) THEN
   DO N=1,N_LAGRANGIAN_CLASSES
      CALL H5FFLUSH_F(HDF_PART_FILE_ID(N),H5F_SCOPE_GLOBAL_F,ERROR)
   ENDDO
ENDIF

END SUBROUTINE FLUSH_VTKHDF


!> \brief Close the boundary VTKHDF file



SUBROUTINE CLOSE_VTKHDF_BNDF()
   CALL CLOSE_VTKHDF_SERIES(HDF_BNDF_FILE_ID,&
      HDF_BNDF_G1,HDF_BNDF_G2,HDF_BNDF_G3,HDF_BNDF_G4,HDF_BNDF_G5,HDF_BNDF_G6,HDF_BNDF_G7)
END SUBROUTINE CLOSE_VTKHDF_BNDF

!> \brief Close the Smoke3D VTKHDF file


SUBROUTINE CLOSE_VTKHDF_SMOKE3D()
   IF (N_SMOKE3D > 0) THEN
      CALL CLOSE_VTKHDF_SERIES(HDF_SM3D_FILE_ID,&
         HDF_SM3D_G1,HDF_SM3D_G2,HDF_SM3D_G3,HDF_SM3D_G4,HDF_SM3D_G5,HDF_SM3D_G6,HDF_SM3D_G7)
   ENDIF
END SUBROUTINE CLOSE_VTKHDF_SMOKE3D

!> \brief Close the slice VTKHDF files


SUBROUTINE CLOSE_VTKHDF_SLICE()
INTEGER :: IQ
   DO IQ=1,MESHES(1)%N_UNIQUE_SLCF
      CALL CLOSE_VTKHDF_SERIES(HDF_SLCF_FILE_ID(IQ),&
         HDF_SLCF_G1(IQ),HDF_SLCF_G2(IQ),HDF_SLCF_G3(IQ),HDF_SLCF_G4(IQ),&
         HDF_SLCF_G5(IQ),HDF_SLCF_G6(IQ),HDF_SLCF_G7(IQ))
   ENDDO
   DEALLOCATE(HDF_SLCF_FILE_ID)
   DEALLOCATE(HDF_SLCF_PLIST_ID)
   DEALLOCATE(HDF_SLCF_CRP_LIST)
   DEALLOCATE(HDF_SLCF_G1)
   DEALLOCATE(HDF_SLCF_G2)
   DEALLOCATE(HDF_SLCF_G3)
   DEALLOCATE(HDF_SLCF_G4)
   DEALLOCATE(HDF_SLCF_G5)
   DEALLOCATE(HDF_SLCF_G6)
   DEALLOCATE(HDF_SLCF_G7)
END SUBROUTINE CLOSE_VTKHDF_SLICE

!> \brief Reopen the boundary VTKHDF file, when the files are not held open


SUBROUTINE OPEN_VTKHDF_BNDF()
   CHARACTER(FN_LENGTH) :: FILENAME
   WRITE(FILENAME,'(A,A,A)') "",TRIM(VTK_DIR)//TRIM(CHID),'_BNDF.vtkhdf'
   ! Initialize HDF5 datafiles
   CALL OPEN_VTKHDF_SERIES(FILENAME,HDF_BNDF_FILE_ID,HDF_BNDF_PLIST_ID,&
      HDF_BNDF_G1,HDF_BNDF_G2,HDF_BNDF_G3,HDF_BNDF_G4,HDF_BNDF_G5,HDF_BNDF_G6,HDF_BNDF_G7)
END SUBROUTINE OPEN_VTKHDF_BNDF

!> \brief Reopen the Smoke3D VTKHDF file, when the files are not held open


SUBROUTINE OPEN_VTKHDF_SMOKE3D()
   CHARACTER(FN_LENGTH) :: FILENAME
   WRITE(FILENAME,'(A,A,A)') "",TRIM(VTK_DIR)//TRIM(CHID),'_SM3D.vtkhdf'
   CALL OPEN_VTKHDF_SERIES(FILENAME,HDF_SM3D_FILE_ID,HDF_SM3D_PLIST_ID,&
      HDF_SM3D_G1,HDF_SM3D_G2,HDF_SM3D_G3,HDF_SM3D_G4,HDF_SM3D_G5,HDF_SM3D_G6,HDF_SM3D_G7)
END SUBROUTINE OPEN_VTKHDF_SMOKE3D

!> \brief Reopen the slice VTKHDF files, when the files are not held open


SUBROUTINE OPEN_VTKHDF_SLICE()
CHARACTER(FN_LENGTH) :: FILENAME,SLCFNAME
INTEGER :: IQ
ALLOCATE(HDF_SLCF_FILE_ID(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_PLIST_ID(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_CRP_LIST(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G1(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G2(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G3(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G4(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G5(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G6(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G7(MESHES(1)%N_UNIQUE_SLCF))
   UNIQUE_LOOPF: DO IQ=1,MESHES(1)%N_UNIQUE_SLCF
      SLCFNAME = MESHES(1)%UNIQUE_SLICE_NAMES(IQ)
      WRITE(FILENAME,'(A,A,A,A)') TRIM(VTK_DIR)//TRIM(CHID),'_',&
         TRIM(SLCFNAME),'.vtkhdf'
      CALL OPEN_VTKHDF_SERIES(FILENAME,HDF_SLCF_FILE_ID(IQ),HDF_SLCF_PLIST_ID(IQ),&
         HDF_SLCF_G1(IQ),HDF_SLCF_G2(IQ),HDF_SLCF_G3(IQ),HDF_SLCF_G4(IQ),&
         HDF_SLCF_G5(IQ),HDF_SLCF_G6(IQ),HDF_SLCF_G7(IQ))
   ENDDO UNIQUE_LOOPF
END SUBROUTINE OPEN_VTKHDF_SLICE

!> \brief Create one slice plane's VTKHDF file and write the grid it will refer to

!>

!> \param FILENAME Name of the file to create

!> \param SLCFNAME Name of the slice plane, shared by every quantity on it

!> \param II Index of the unique slice plane

!> \param NTSL Terrain slice counter, used to index K_AGL_SLICE

!>

!> Slices that lie on the same plane share a file, one quantity per point data

!> array, so the grid is written once no matter how many quantities there are.


SUBROUTINE INITIALIZE_VTKHDF_SLICE(FILENAME,SLCFNAME,II,NTSL)
   CHARACTER(*), INTENT(IN) :: FILENAME,SLCFNAME
   INTEGER, INTENT(IN) :: II, NTSL
   INTEGER(HID_T) :: DSET_ID_NCELLS, DSET_ID_NCON, DSET_ID_NPTS    ! Dataset identifiers
   INTEGER :: NPOINTS, NCELLS, IQ, NQT
   INTEGER :: NM, ERROR, NX, NY, NZ, NCONNECTIONS
   TYPE(SLICE_TYPE), POINTER :: SL
   INTEGER, DIMENSION(1:N_MPI_PROCESSES) :: MESHES_PER_PROCESS
   INTEGER :: N_WRITTEN
   INTEGER(HSIZE_T), DIMENSION(1) :: START1(1), EXTENT1(1)
   INTEGER(IB32), DIMENSION(1) :: IDATA_OUT_32(1)
   LOGICAL :: SL3D
   INTEGER(IB32), DIMENSION(1:NMESHES) :: NCELLS_ALL, NPOINTS_ALL
   INTEGER :: NCONNECTIONS_ALL
   INTEGER :: NINFO(2)
   TYPE (MPI_STATUS) :: MPISTATUS
   INTEGER :: IPROC, IERR
   
   DO NM=1,N_MPI_PROCESSES
      MESHES_PER_PROCESS(NM) = 0
   ENDDO
   DO NM=1,NMESHES
      MESHES_PER_PROCESS(PROCESS(NM)+1) = MESHES_PER_PROCESS(PROCESS(NM)+1) + 1
   ENDDO
   
   ! Initialize HDF5 datafiles
   CALL CREATE_OPEN_VTKHDF_SERIES(FILENAME,HDF_SLCF_FILE_ID(II),HDF_SLCF_PLIST_ID(II),&
      HDF_SLCF_G1(II),HDF_SLCF_G2(II),HDF_SLCF_G3(II),HDF_SLCF_G4(II),&
      HDF_SLCF_G5(II),HDF_SLCF_G6(II),HDF_SLCF_G7(II))
   
   ! Initialize file
   START1 = NMESHES
   EXTENT1 = NMESHES
   CALL PARALLEL_INIT_I32(HDF_SLCF_G1(II), TRIM("NumberOfCells"), HDF_SLCF_CRP_LIST(II), 1,&
      START1, EXTENT1, DSET_ID_NCELLS, HDF_SLCF_PLIST_ID(II)) ! NumberOfCells
   CALL PARALLEL_INIT_I32(HDF_SLCF_G1(II), TRIM("NumberOfPoints"), HDF_SLCF_CRP_LIST(II), 1,&
      START1, EXTENT1, DSET_ID_NPTS, HDF_SLCF_PLIST_ID(II)) ! NumberOfPoints
   CALL PARALLEL_INIT_I32(HDF_SLCF_G1(II), TRIM("NumberOfConnectivityIds"), HDF_SLCF_CRP_LIST(II), 1,&
      START1, EXTENT1, DSET_ID_NCON, HDF_SLCF_PLIST_ID(II)) ! NumberOfConnectivityIds
   
   NQT = MESHES(1)%N_UNIQUE_SLCF
   ! Fill metadata
   N_WRITTEN=0
   MESH_LOOP: DO NM=1,NMESHES
      CALL POINT_TO_MESH(NM)
      IF (PROCESS(NM)/=MY_RANK) CYCLE MESH_LOOP
      NCELLS = 0
      NPOINTS = 0
      IF (MESHES(NM)%EMPTY_UNIQUE_SLICE(II)) THEN
         NCELLS = 0
         NPOINTS = 0
         NCONNECTIONS = 0
      ELSE
         QUANTITY_LOOP: DO IQ=1,MESHES(1)%N_SLCF_VTK
            SL => SLICE(IQ)
            IF (TRIM(SL%SLCF_NAME).NE.TRIM(SLCFNAME)) CYCLE QUANTITY_LOOP
            NX = SL%I2 + 1 - SL%I1
            NY = SL%J2 + 1 - SL%J1
            NZ = SL%K2 + 1 - SL%K1
            IF (SL%I2-SL%I1==0 .OR. SL%J2-SL%J1==0 .OR. SL%K2-SL%K1==0) THEN
               NCONNECTIONS=4 ! 2-D slice
            ELSEIF (MESHES(1)%UNIQUE_SLCF_AGL(II)>0) THEN
               NCONNECTIONS=4 ! 2-D slice
               NZ=1
            ELSE
               NCONNECTIONS=8 ! 3-D slice
            ENDIF
            NPOINTS = NX*NY*NZ
            NCELLS = MAX((NX-1),1)*MAX((NY-1),1)*MAX((NZ-1),1)
            EXIT
         ENDDO QUANTITY_LOOP
      ENDIF

      ! Write number of cells data to file
      START1 = NM-1
      EXTENT1 = 1
      IDATA_OUT_32 = NCELLS
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, HDF_SLCF_PLIST_ID(II), START1, EXTENT1, IDATA_OUT_32)
      HDF_SLCF_G1_NCELLS(II,NM) = INT(NCELLS,IB32)
      
      ! Write number of points data to file
      START1 = NM-1
      EXTENT1 = 1
      IDATA_OUT_32 = NPOINTS
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, HDF_SLCF_PLIST_ID(II), START1, EXTENT1, IDATA_OUT_32)
      HDF_SLCF_G1_NPOINTS(II,NM) = INT(NPOINTS,IB32)
      
      ! Write NumberOfConnectivityIds data to file
      START1 = NM-1
      EXTENT1 = 1
      IDATA_OUT_32 = NCELLS*NCONNECTIONS
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, HDF_SLCF_PLIST_ID(II), START1, EXTENT1, IDATA_OUT_32)
      
      N_WRITTEN = N_WRITTEN + 1
   ENDDO MESH_LOOP
   
   ! Write fake data in processes that have less meshes than max process
   DO NM=1,MAXVAL(MESHES_PER_PROCESS)
      IF (N_WRITTEN < MAXVAL(MESHES_PER_PROCESS)) THEN
         NCELLS = 0
         NPOINTS = 0
         START1 = NM-1
         EXTENT1 = 0
         IDATA_OUT_32 = NCELLS
         CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, HDF_SLCF_PLIST_ID(II), START1, EXTENT1, IDATA_OUT_32)
      
         ! Write number of points data to file
         START1 = NM-1
         EXTENT1 = 0
         IDATA_OUT_32 = NPOINTS
         CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, HDF_SLCF_PLIST_ID(II), START1, EXTENT1, IDATA_OUT_32)
      
         ! Write NumberOfConnectivityIds data to file
         START1 = NM-1
         EXTENT1 = 0
         IDATA_OUT_32 = NCELLS*NCONNECTIONS
         CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, HDF_SLCF_PLIST_ID(II), START1, EXTENT1, IDATA_OUT_32)
         N_WRITTEN = N_WRITTEN + 1
      ENDIF
   ENDDO
   
   ! Write geometry data to file
   N_WRITTEN = 0
   MESH_LOOP2: DO NM=1,NMESHES
      CALL POINT_TO_MESH(NM)
      IF (PROCESS(NM)/=MY_RANK) CYCLE MESH_LOOP2
      SL3D = MESHES(1)%UNIQUE_SLICE_IS_SL3D(II)
      CALL WRITE_VTKHDF_SLICE_CELL_FILE_NOOPEN(SLCFNAME,SL3D,NM,NCELLS_ALL,NPOINTS_ALL,NCONNECTIONS_ALL,NTSL,&
         HDF_SLCF_PLIST_ID(II),HDF_SLCF_CRP_LIST(II),HDF_SLCF_G1(II))
      N_WRITTEN = N_WRITTEN + 1
   ENDDO MESH_LOOP2
   
   ! Write fake data in processes that have less meshes than max process
   DO NM=1,MAXVAL(MESHES_PER_PROCESS)
      IF (N_WRITTEN < MAXVAL(MESHES_PER_PROCESS)) THEN
         CALL WRITE_VTKHDF_SLICE_CELL_FILE_NOOPEN(SLCFNAME,SL3D,NM,NCELLS_ALL,NPOINTS_ALL,NCONNECTIONS_ALL,NTSL,&
            HDF_SLCF_PLIST_ID(II),HDF_SLCF_CRP_LIST(II),HDF_SLCF_G1(II))
         N_WRITTEN = N_WRITTEN + 1
      ENDIF
   ENDDO
   
   ! Close VTKHDF interface
   CALL VTK_DSET_RELEASE(DSET_ID_NCELLS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NPTS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NCON, ERROR)
   !CALL H5PCLOSE_F(HDF_SLCF_PLIST_ID(ii), ERROR)
   
   ! Exchange HDF_SLCF_G1_NCELLS(II,1:NMESHES)
   EXCHANGE_SLCF_NCELLS: DO NM=1,NMESHES
      IF (PROCESS(NM)/=MY_RANK) THEN
         CALL MPI_RECV(NINFO,2,MPI_INTEGER,PROCESS(NM),PROCESS(NM),MPI_COMM_WORLD,MPISTATUS,IERR)
         HDF_SLCF_G1_NCELLS(II,NM) = INT(NINFO(1),IB32)
         HDF_SLCF_G1_NPOINTS(II,NM) = INT(NINFO(2),IB32)
         CYCLE EXCHANGE_SLCF_NCELLS
      ELSE
         NINFO = INT((/HDF_SLCF_G1_NCELLS(II,NM),HDF_SLCF_G1_NPOINTS(II,NM)/))
         DO IPROC=0,N_MPI_PROCESSES-1
            IF (MY_RANK/=IPROC) THEN
               CALL MPI_SEND(NINFO,2,MPI_INTEGER,IPROC,MY_RANK,MPI_COMM_WORLD,IERR)
            ENDIF
         ENDDO
      ENDIF
   ENDDO EXCHANGE_SLCF_NCELLS
   
END SUBROUTINE INITIALIZE_VTKHDF_SLICE


!> \brief Write the obstruction and terrain geometry to a single time VTKHDF file


!>


!> The geometry is written once, before the time loop, and is what ParaView draws


!> as the solid part of the scene.



SUBROUTINE WRITE_VTKHDF_GEOM_FILE()
   USE MPI_F08
   USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
   INTEGER(HID_T) :: FILE_ID, PLIST_ID, CRP_LIST       ! Identifiers
   INTEGER(HID_T) :: GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4 ! Group identifier
   INTEGER(HID_T) :: DSET_ID_NCELLS, DSET_ID_NCON, DSET_ID_NPTS, DSET_ID    ! Dataset identifiers
   INTEGER(HID_T) :: DSET_ID_CON, DSET_ID_OFF, DSET_ID_PTS, DSET_ID_TYP
   INTEGER :: IP, II, LAST_OFFSET_VALUE
   INTEGER :: NFACES, NFACES_CUTCELLS, NVERTS, NVERTS_CUTCELLS
   INTEGER :: NM, NMNM, NM1, NM2, ERROR, PA_NCELLS, PA_NPOINTS
   CHARACTER(200) :: FILENAME
   INTEGER, DIMENSION(1:N_MPI_PROCESSES) :: MESHES_PER_PROCESS
   INTEGER :: N_WRITTEN, I, IERR, IPROC
   TYPE (MPI_STATUS) :: MPISTATUS
   TYPE(PATCH_TYPE), POINTER :: PA
   TYPE(PATCH_TYPE), TARGET :: OB_PATCH
   INTEGER :: IOR, NOB
   INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCATIONS,FACES,SURFIND,GEOMIND
   REAL(FB), ALLOCATABLE, DIMENSION(:) :: VERTS
   REAL(FB), ALLOCATABLE, DIMENSION(:) :: X_PTS, Y_PTS, Z_PTS
   INTEGER :: NCELLS_ACCUM, NCELLS_MAX, NCELLS_TOTAL, NCELLS_START
   INTEGER :: NCONN_ACCUM, NCONN_MAX, NCONN_TOTAL, NCONN_START
   INTEGER :: NOFFSETS_ACCUM, NOFFSETS_MAX, NOFFSETS_TOTAL, NOFFSETS_START
   INTEGER :: NPOINTS_ACCUM, NPOINTS_MAX, NPOINTS_TOTAL, NPOINTS_START
   INTEGER(IB32), DIMENSION(1:2*NMESHES) :: NCELLS, NPOINTS
   INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS, ALL_CONNECT, ALL_OFFSETS
   INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: VTKC_TYPE, ALL_VTKC_TYPE
   REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: VERTICES, ALL_VERTICES, COLORS, ALL_COLORS
   TYPE(MESH_TYPE), POINTER :: M !=>NULL()
   TYPE(OBSTRUCTION_TYPE), POINTER :: OB !=>NULL()
   REAL(FB) :: COLOR(3)
   INTEGER :: NINFO(6)
   INTEGER(HSIZE_T), DIMENSION(1) :: CDIM1(1), DDIM1(1), CDIM2(2), DDIM2(2)
   INTEGER(IB32), DIMENSION(1) :: DATA_OUT(1)
   INTEGER(HSIZE_T), DIMENSION(1) :: START1(1), EXTENT1(1), START2(2), EXTENT2(2)
   REAL(FB), DIMENSION(1) :: DATA_OUT2(2)
   INTEGER(IB8), DIMENSION(1) :: DATA_OUT_U8(1)
   REAL(EB) :: TNOW
   TNOW = CURRENT_TIME()

   ALLOCATE(NCELLS_VTK(1:2*NMESHES))
   ALLOCATE(NPOINTS_VTK(1:2*NMESHES))
   ALLOCATE(NCONNECTIONS_VTK(1:2*NMESHES))

   DO NM=1,N_MPI_PROCESSES
      MESHES_PER_PROCESS(NM) = 0
   ENDDO
   DO NM=1,NMESHES
      MESHES_PER_PROCESS(PROCESS(NM)+1) = MESHES_PER_PROCESS(PROCESS(NM)+1) + 1
      NM1 = 2*NM-1
      NM2 = 2*NM
      NCELLS_VTK(NM1) = 0
      NCELLS_VTK(NM2) = 0
      NPOINTS_VTK(NM1) = 0
      NPOINTS_VTK(NM2) = 0
      NCONNECTIONS_VTK(NM1) = 0
      NCONNECTIONS_VTK(NM2) = 0
   ENDDO
   
   ! Initialize file
   WRITE(FILENAME,'(A,A,A)') "",TRIM(CHID),'_GEOM.vtkhdf'
   CALL CREATE_OPEN_VTKHDF(FILENAME, FILE_ID, PLIST_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4)

   DDIM1 = 2*NMESHES
   CDIM1 = 2*NMESHES
   ! 2*NMESHES arrays. odds include OBST patches, evens include GEOM patches
   CALL PARALLEL_INIT_I32(GROUP_ID1, TRIM("NumberOfCells"), CRP_LIST, 1,&
      DDIM1, CDIM1, DSET_ID_NCELLS, PLIST_ID) ! NumberOfCells
   CALL PARALLEL_INIT_I32(GROUP_ID1, TRIM("NumberOfPoints"), CRP_LIST, 1,&
      DDIM1, CDIM1, DSET_ID_NPTS, PLIST_ID) ! NumberOfPoints
   CALL PARALLEL_INIT_I32(GROUP_ID1, TRIM("NumberOfConnectivityIds"), CRP_LIST, 1,&
      DDIM1, CDIM1, DSET_ID_NCON, PLIST_ID) ! NumberOfConnectivityIds
      
   ! Fill metadata
   N_WRITTEN=0
   MESH_LOOP: DO NM=1,NMESHES
      CALL POINT_TO_MESH(NM)
      M => MESHES(NM)
      NM1 = 2*NM-1
      NM2 = 2*NM
      IF (PROCESS(NM)/=MY_RANK) THEN
         CYCLE MESH_LOOP
      ENDIF
      NCELLS_ACCUM = 0
      NPOINTS_ACCUM = 0
      NCONN_ACCUM = 0
      ! Count OBST face info
      OBST_LOOP1: DO NOB=1,M%N_OBST
         OB => M%OBSTRUCTION(NOB)
         FACE_LOOP1: DO IOR=-3,3
            IF (IOR==0) CYCLE FACE_LOOP1
            IF (.NOT.VTK_OBST_FACE_DRAWN(OB,IOR)) CYCLE FACE_LOOP1
            CALL SET_VTK_OBST_FACE_PATCH(OB,IOR,OB_PATCH)
            PA => OB_PATCH
            ! Initialize piece
            CALL BUILD_VTK_SOLID_PHASE_GEOMETRY(NM, PA, NCELLS(NM1), NPOINTS(NM1),&
               X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)
            NCONN_ACCUM = NCONN_ACCUM + NCELLS(NM1)*4
            NCELLS_ACCUM = NCELLS_ACCUM + NCELLS(NM1)
            NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS(NM1)
            CALL DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)
         ENDDO FACE_LOOP1
      ENDDO OBST_LOOP1
      
      NCELLS_VTK(NM1) = NCELLS_ACCUM
      NPOINTS_VTK(NM1) = NPOINTS_ACCUM
      NCONNECTIONS_VTK(NM1) = NCONN_ACCUM
      DDIM1 = NM1-1
      CDIM1 = 1
      ! Write number of cells data to file
      DATA_OUT = NCELLS_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
      
      ! Write number of points data to file
      DATA_OUT = NPOINTS_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
      
      ! Write NumberOfConnectivityIds data to file
      DATA_OUT = NCONN_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
      
      ! Count GEOM patch info
      NCELLS_ACCUM = 0
      NPOINTS_ACCUM = 0
      NCONN_ACCUM = 0
      IF (MESHES(NM)%N_INTERNAL_CFACE_CELLS>0) THEN
         CALL GET_GEOMSIZES('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS)
         IF (NVERTS>0 .AND. NFACES>0) THEN
            ALLOCATE(VERTS(3*NVERTS))
            ALLOCATE(FACES(3*NFACES))
            ALLOCATE(LOCATIONS(NFACES))
            ALLOCATE(SURFIND(NFACES))
            ALLOCATE(GEOMIND(NFACES))
            CALL GET_GEOMINFO('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS,&
                              VERTS,FACES,LOCATIONS,SURFIND=SURFIND,GEOMIND=GEOMIND)
            CALL BUILD_VTK_GEOM_GEOMETRY(VERTS, FACES, NCELLS(NM2), NPOINTS(NM2), X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)
            DEALLOCATE(VERTS)
            DEALLOCATE(FACES)
            DEALLOCATE(LOCATIONS)
            DEALLOCATE(SURFIND)
            DEALLOCATE(GEOMIND)
            
            CALL DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)
            
            NCONN_ACCUM = NCONN_ACCUM + NCELLS(NM2)*3
            NCELLS_ACCUM = NCELLS_ACCUM + NCELLS(NM2)
            NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS(NM2)
         ENDIF
      ENDIF
      
      NCELLS_VTK(NM2) = NCELLS_ACCUM
      NPOINTS_VTK(NM2) = NPOINTS_ACCUM
      NCONNECTIONS_VTK(NM2) = NCONN_ACCUM
      DDIM1 = NM2-1
      CDIM1 = 1
      ! Write number of cells data to file
      DATA_OUT = (/NCELLS_ACCUM/)
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
      
      ! Write number of points data to file
      DATA_OUT = (/NPOINTS_ACCUM/)
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
      
      ! Write NumberOfConnectivityIds data to file
      DATA_OUT = (/NCONN_ACCUM/)
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
      
      N_WRITTEN = N_WRITTEN + 1
   ENDDO MESH_LOOP

   ! Write fake data in processes that have less meshes than max process
   DO NM=1,MAXVAL(MESHES_PER_PROCESS)
      IF (N_WRITTEN < MAXVAL(MESHES_PER_PROCESS)) THEN
         DDIM1 = NM1-1
         CDIM1 = 0
         DATA_OUT = 0
         DO I=1,2 ! Loop through 2x to write empty twice per mesh
            CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
         
            ! Write number of points data to file
            CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
         
            ! Write NumberOfConnectivityIds data to file
            CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, PLIST_ID, DDIM1, CDIM1, DATA_OUT)
         ENDDO
         N_WRITTEN = N_WRITTEN + 1
      ENDIF
   ENDDO

   ! Close VTKHDF interface
   !CALL MPI_BARRIER(MPI_COMM_WORLD, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NCELLS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NPTS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NCON, ERROR)
   
   CALL CLOSE_VTKHDF(FILE_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4)

   ! Exchange boundary NCELLS, NPOINTS, and NCONNECTIONS
   MESH_LOOP2: DO NM=1,NMESHES
      NM1 = 2*NM-1
      NM2 = 2*NM
      IF (PROCESS(NM)/=MY_RANK) THEN
         CALL MPI_RECV(NINFO,6,MPI_INTEGER,PROCESS(NM),PROCESS(NM),MPI_COMM_WORLD,MPISTATUS,IERR)
         NCELLS_VTK(NM1) = INT(NINFO(1),IB32)
         NPOINTS_VTK(NM1) = INT(NINFO(2),IB32)
         NCONNECTIONS_VTK(NM1) = INT(NINFO(3),IB32)
         NCELLS_VTK(NM2) = INT(NINFO(4),IB32)
         NPOINTS_VTK(NM2) = INT(NINFO(5),IB32)
         NCONNECTIONS_VTK(NM2) = INT(NINFO(6),IB32)
         CYCLE MESH_LOOP2
      ELSE
         NINFO = INT((/NCELLS_VTK(NM1),NPOINTS_VTK(NM1),NCONNECTIONS_VTK(NM1),&
                      NCELLS_VTK(NM2),NPOINTS_VTK(NM2),NCONNECTIONS_VTK(NM2)/))
         DO IPROC=0,N_MPI_PROCESSES-1
            IF (MY_RANK/=IPROC) THEN
               CALL MPI_SEND(NINFO,6,MPI_INTEGER,IPROC,MY_RANK,MPI_COMM_WORLD,IERR)
            ENDIF
         ENDDO
      ENDIF
   ENDDO MESH_LOOP2
   
   ! Reopen VTKHDF INTERFACE
   CALL OPEN_VTKHDF(FILENAME, FILE_ID, PLIST_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4)
   
   NPOINTS_TOTAL=0
   NCELLS_TOTAL=0
   NCELLS_ACCUM=0
   NPOINTS_ACCUM=0
   NOFFSETS_TOTAL=0
   NOFFSETS_ACCUM = 0
   NCONN_TOTAL = 0
   NCONN_ACCUM = 0
   NCELLS_MAX = 1
   NPOINTS_MAX = 1
   NCONN_MAX = 1
   NOFFSETS_MAX = 1
   NPOINTS_START = 0
   
   MESH_LOOP_HDF_COUNT: DO NMNM=1,2*NMESHES
      NPOINTS_TOTAL = NPOINTS_TOTAL + NPOINTS_VTK(NMNM)
      NCELLS_TOTAL = NCELLS_TOTAL + NCELLS_VTK(NMNM)
      NOFFSETS_TOTAL = NOFFSETS_TOTAL + NCELLS_VTK(NMNM) + 1
      NCONN_TOTAL = NCONN_TOTAL + NCONNECTIONS_VTK(NMNM)
      NCELLS_MAX = MAX(NCELLS_MAX,NCELLS_VTK(NMNM))
      NPOINTS_MAX = MAX(NPOINTS_MAX,NPOINTS_VTK(NMNM))
      NOFFSETS_MAX = MAX(NOFFSETS_MAX,NCELLS_VTK(NMNM) + 1)
      NCONN_MAX = MAX(NCONN_MAX,NCONNECTIONS_VTK(NMNM))
   ENDDO MESH_LOOP_HDF_COUNT
   
   DDIM1 = NCONN_MAX
   CDIM1 = NCONN_TOTAL
   CALL PARALLEL_INIT_I32(GROUP_ID1, TRIM("Connectivity"), CRP_LIST, 1, DDIM1,CDIM1, DSET_ID_CON, PLIST_ID) ! Connectivity
   DDIM1 = NOFFSETS_MAX
   CDIM1 = NOFFSETS_TOTAL
   CALL PARALLEL_INIT_I32(GROUP_ID1, TRIM("Offsets"), CRP_LIST, 1, DDIM1,CDIM1, DSET_ID_OFF, PLIST_ID) ! Offsets
   DDIM2(1:2) = (/3_HSIZE_T, INT(NPOINTS_MAX,HSIZE_T)/)
   CDIM2(1:2) = (/3_HSIZE_T, INT(NPOINTS_TOTAL,HSIZE_T)/)
   CALL PARALLEL_INIT_F32(GROUP_ID1, TRIM("Points"), CRP_LIST, 2, DDIM2,CDIM2, DSET_ID_PTS, PLIST_ID) ! Points
   DDIM1 = NCELLS_MAX
   CDIM1 = NCELLS_TOTAL
   CALL PARALLEL_INIT_U8(GROUP_ID1, TRIM("Types"), CRP_LIST, 1, DDIM1,CDIM1, DSET_ID_TYP, PLIST_ID) ! Types
   
   CALL POINT_TO_MESH(1)
   DDIM2(1:2) = (/3_HSIZE_T, INT(NCELLS_MAX,HSIZE_T)/)
   CDIM2(1:2) = (/3_HSIZE_T, INT(NCELLS_TOTAL,HSIZE_T)/)
   CALL PARALLEL_INIT_F32(GROUP_ID2, "Color", CRP_LIST, 2, DDIM2,CDIM2, DSET_ID, PLIST_ID) ! Data
   
   ! Fill metadata
   N_WRITTEN=0
   MESH_LOOP_HDF: DO NM=1,NMESHES
      NM1 = 2*NM-1
      NM2 = 2*NM
      NPOINTS_START = NPOINTS_ACCUM
      NCELLS_START = NCELLS_ACCUM
      NCONN_START = NCONN_ACCUM
      NOFFSETS_START = NOFFSETS_ACCUM
      LAST_OFFSET_VALUE = 0
      IF (PROCESS(NM)/=MY_RANK) THEN
         NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS_VTK(NM1) + NPOINTS_VTK(NM2)
         NCELLS_ACCUM = NCELLS_ACCUM + NCELLS_VTK(NM1) + NCELLS_VTK(NM2)
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS_VTK(NM1)+1
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS_VTK(NM2)+1
         NCONN_ACCUM = NCONN_ACCUM + NCONNECTIONS_VTK(NM1) + NCONNECTIONS_VTK(NM2)
         CYCLE MESH_LOOP_HDF
      ENDIF
      CALL POINT_TO_MESH(NM)
      ! Build OBST boundary geometry
      IF (NCELLS_VTK(NM1)>0) THEN
         ALLOCATE(ALL_VERTICES(3,NPOINTS_VTK(NM1)))
         ALLOCATE(ALL_CONNECT(NCONNECTIONS_VTK(NM1)))
         ALLOCATE(ALL_OFFSETS(NCELLS_VTK(NM1)+1))
         ALLOCATE(ALL_VTKC_TYPE(NCELLS_VTK(NM1)))
         ALLOCATE(ALL_COLORS(3,NCELLS_VTK(NM1)))
         ALL_OFFSETS(NOFFSETS_ACCUM-NOFFSETS_START+1) = 0
         M => MESHES(NM)
         OBST_LOOP2: DO NOB=1,M%N_OBST
            OB => M%OBSTRUCTION(NOB)
            FACE_LOOP2: DO IOR=-3,3
               IF (IOR==0) CYCLE FACE_LOOP2
               IF (.NOT.VTK_OBST_FACE_DRAWN(OB,IOR)) CYCLE FACE_LOOP2
               CALL SET_VTK_OBST_FACE_PATCH(OB,IOR,OB_PATCH)
               PA => OB_PATCH
               IF (OB%RGB(1)==-1) THEN
                  COLOR = REAL(SURFACE(OB%SURF_INDEX(IOR))%RGB,FB)/255._FB
               ELSE
                  COLOR = REAL(OB%RGB,FB)/255._FB
               ENDIF

               ! Initialize piece
               CALL BUILD_VTK_SOLID_PHASE_GEOMETRY(NM, PA, PA_NCELLS, PA_NPOINTS,&
                  X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)
               ALLOCATE(VERTICES(3,PA_NPOINTS))
               DO II=1,PA_NPOINTS
                  VERTICES(1:3,II) = (/X_PTS(II),Y_PTS(II),Z_PTS(II)/)
               ENDDO
               ALLOCATE(COLORS(3,PA_NCELLS))
               DO II=1,PA_NCELLS
                  COLORS(1:3,II) = COLOR
               ENDDO
               ALL_VERTICES(:,NPOINTS_ACCUM-NPOINTS_START+1:NPOINTS_ACCUM-NPOINTS_START+PA_NPOINTS) = VERTICES
               ALL_CONNECT(NCONN_ACCUM-NCONN_START+1:NCONN_ACCUM-NCONN_START+PA_NCELLS*4) = CONNECT + NPOINTS_ACCUM-NPOINTS_START
               OFFSETS = OFFSETS + LAST_OFFSET_VALUE
               LAST_OFFSET_VALUE = OFFSETS(SIZE(OFFSETS))
               ALL_OFFSETS(NOFFSETS_ACCUM-NOFFSETS_START+2:NOFFSETS_ACCUM-NOFFSETS_START+PA_NCELLS+1) = OFFSETS
               ALL_VTKC_TYPE(NCELLS_ACCUM-NCELLS_START+1:NCELLS_ACCUM-NCELLS_START+PA_NCELLS) = VTKC_TYPE
               ALL_COLORS(:,NCELLS_ACCUM-NCELLS_START+1:NCELLS_ACCUM-NCELLS_START+PA_NCELLS) = COLORS
               NCONN_ACCUM = NCONN_ACCUM + PA_NCELLS*4
               NCELLS_ACCUM = NCELLS_ACCUM + PA_NCELLS
               NPOINTS_ACCUM = NPOINTS_ACCUM + PA_NPOINTS
               NOFFSETS_ACCUM = NOFFSETS_ACCUM + PA_NCELLS
               CALL DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)
               DEALLOCATE(VERTICES)
               DEALLOCATE(COLORS)
            ENDDO FACE_LOOP2
         ENDDO OBST_LOOP2
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + 1
         
         ! Write connectivity data to file
         START1 = INT(NCONN_START,HSIZE_T)
         EXTENT1 = INT(NCONNECTIONS_VTK(NM1),HSIZE_T)
         CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, PLIST_ID, START1, EXTENT1, ALL_CONNECT)
            
         ! Write offsets data to file
         START1 = INT(NOFFSETS_START,HSIZE_T)
         EXTENT1 = INT(NCELLS_VTK(NM1)+1,HSIZE_T)
         CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, PLIST_ID, START1, EXTENT1, ALL_OFFSETS)
         
         ! Write point data to file
         START2(1:2) = INT((/0,NPOINTS_START/),HSIZE_T)
         EXTENT2(1:2) = INT((/3, NPOINTS_VTK(NM1)/),HSIZE_T)
         CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, PLIST_ID, START2, EXTENT2, ALL_VERTICES)
         
         ! Write types data to file
         START1 = INT(NCELLS_START,HSIZE_T)
         EXTENT1 = INT(NCELLS_VTK(NM1),HSIZE_T)
         CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, PLIST_ID, START1, EXTENT1, ALL_VTKC_TYPE)
         
         ! Write color data to file
         START2(1:2) = INT((/0,NCELLS_START/),HSIZE_T)
         EXTENT2(1:2) = INT((/3, NCELLS_VTK(NM1)/),HSIZE_T)
         CALL PARALLEL_WRITE_F32(2, DSET_ID, PLIST_ID, START2, EXTENT2, ALL_COLORS)
         
         DEALLOCATE(ALL_VERTICES)
         DEALLOCATE(ALL_CONNECT)
         DEALLOCATE(ALL_OFFSETS)
         DEALLOCATE(ALL_VTKC_TYPE)
         DEALLOCATE(ALL_COLORS)
         N_WRITTEN = N_WRITTEN + 1
      ELSE
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + 1
      ENDIF
      ! Build GEOM patch geometry
      IF ((MESHES(NM)%N_INTERNAL_CFACE_CELLS>0).AND.(.TRUE.)) THEN
         NPOINTS_START = NPOINTS_ACCUM
         NCELLS_START = NCELLS_ACCUM
         NCONN_START = NCONN_ACCUM
         NOFFSETS_START = NOFFSETS_ACCUM
         CALL GET_GEOMSIZES('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS)
         IF (NVERTS>0 .AND. NFACES>0) THEN
            ALLOCATE(ALL_OFFSETS(NCELLS_VTK(NM2)+1))
            ALLOCATE(VERTS(3*NVERTS))
            ALLOCATE(FACES(3*NFACES))
            ALLOCATE(LOCATIONS(NFACES))
            ALLOCATE(SURFIND(NFACES))
            ALLOCATE(GEOMIND(NFACES))
            CALL GET_GEOMINFO('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS,&
                              VERTS,FACES,LOCATIONS,SURFIND=SURFIND,GEOMIND=GEOMIND)
            CALL BUILD_VTK_GEOM_GEOMETRY(VERTS, FACES, PA_NCELLS, PA_NPOINTS, X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)
            ALLOCATE(VERTICES(3,PA_NPOINTS))
            DO II=1,PA_NPOINTS
               VERTICES(1:3,II) = (/X_PTS(II),Y_PTS(II),Z_PTS(II)/)
            ENDDO
            ALLOCATE(COLORS(3,PA_NCELLS))
            DO II=1,PA_NCELLS
               COLOR = REAL(SURFACE(SURFIND(II))%RGB,FB)/255._FB
               COLORS(1:3,II) = COLOR
            ENDDO
            
            ALL_OFFSETS(1) = 0
            ALL_OFFSETS(2:NCELLS_VTK(NM2)+1) = OFFSETS
            NCONN_ACCUM = NCONN_ACCUM + PA_NCELLS*3
            NCELLS_ACCUM = NCELLS_ACCUM + PA_NCELLS
            NPOINTS_ACCUM = NPOINTS_ACCUM + PA_NPOINTS
            NOFFSETS_ACCUM = NOFFSETS_ACCUM + PA_NCELLS + 1
            
            ! Write connectivity data to file
            START1 = INT(NCONN_START,HSIZE_T)
            EXTENT1 = INT(NCONNECTIONS_VTK(NM2),HSIZE_T)
            CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, PLIST_ID, START1, EXTENT1, CONNECT)
               
            ! Write offsets data to file
            START1 = INT(NOFFSETS_START+1,HSIZE_T)
            EXTENT1 = INT(NCELLS_VTK(NM2)+1,HSIZE_T)
            CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, PLIST_ID, START1, EXTENT1, ALL_OFFSETS)
            
            ! Write point data to file
            START2(1:2) = INT((/0,NPOINTS_START/),HSIZE_T)
            EXTENT2(1:2) = INT((/3, NPOINTS_VTK(NM2)/),HSIZE_T)
            CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, PLIST_ID, START2, EXTENT2, VERTICES)
            
            ! Write types data to file
            START1 = INT(NCELLS_START,HSIZE_T)
            EXTENT1 = INT(NCELLS_VTK(NM2),HSIZE_T)
            CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, PLIST_ID, START1, EXTENT1, VTKC_TYPE)
            
            ! Write color data to file
            START2(1:2) = INT((/0,NCELLS_START/),HSIZE_T)
            EXTENT2(1:2) = INT((/3, NCELLS_VTK(NM2)/),HSIZE_T)
            CALL PARALLEL_WRITE_F32(2, DSET_ID, PLIST_ID, START2, EXTENT2, COLORS)
            
            DEALLOCATE(VERTS)
            DEALLOCATE(FACES)
            DEALLOCATE(LOCATIONS)
            DEALLOCATE(SURFIND)
            DEALLOCATE(GEOMIND)
            DEALLOCATE(VERTICES)
            DEALLOCATE(ALL_OFFSETS)
            DEALLOCATE(COLORS)
            CALL DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)
            N_WRITTEN = N_WRITTEN + 1
         ENDIF
      ELSE
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + 1
      ENDIF
   ENDDO MESH_LOOP_HDF
   
   ! Write fake data in processes that have less meshes than max process
   DO NM=1,MAXVAL(MESHES_PER_PROCESS)*2
      IF (N_WRITTEN < MAXVAL(MESHES_PER_PROCESS)*2) THEN

         ! Write connectivity data to file
         START1 = 0
         EXTENT1 = 0
         DATA_OUT = 0
         CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, PLIST_ID, START1, EXTENT1, DATA_OUT)
            
         ! Write offsets data to file
         CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, PLIST_ID, START1, EXTENT1, DATA_OUT)
         
         ! Write point data to file
         START2 = 0
         EXTENT2 = 0
         DATA_OUT2 = 0
         CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, PLIST_ID, START2, EXTENT2, DATA_OUT2)
         
         ! Write types data to file
         DATA_OUT_U8 = 0
         CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, PLIST_ID, START1, EXTENT1, DATA_OUT_U8)
         
         ! Write color data to file
         CALL PARALLEL_WRITE_F32(2, DSET_ID, PLIST_ID, START2, EXTENT2, DATA_OUT2)
         
         N_WRITTEN = N_WRITTEN + 1
      ENDIF
   ENDDO
   
   ! Close VTKHDF interface
   CALL VTK_DSET_RELEASE(DSET_ID, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_CON, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_OFF, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_PTS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_TYP, ERROR)
   CALL CLOSE_VTKHDF(FILE_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4)
   T_USED(7) = T_USED(7) + CURRENT_TIME() - TNOW
END SUBROUTINE WRITE_VTKHDF_GEOM_FILE



!> \brief Create the boundary VTKHDF file and write the patch geometry it will refer to



!>



!> \param FILENAME Name of the file to create




SUBROUTINE INITIALIZE_VTKHDF_BNDF(FILENAME)
   CHARACTER(*), INTENT(IN) :: FILENAME
   INTEGER(HID_T) :: DSET_ID_NCELLS, DSET_ID_NCON, DSET_ID_NPTS    ! Dataset identifiers
   INTEGER(HID_T) :: DSET_ID, DSET_ID_CON, DSET_ID_OFF, DSET_ID_PTS, DSET_ID_TYP       ! Dataset identifiers
   INTEGER :: PA_NPOINTS, PA_NCELLS, II, IP
   INTEGER :: NPOINTS, NCELLS, NPOINTS_ACCUM, NCELLS_ACCUM, NCONNECTIONS_ACCUM
   INTEGER :: NFACES, NFACES_CUTCELLS, NVERTS, NVERTS_CUTCELLS
   INTEGER :: N, NM, NM1, NM2, NMNM, ERROR, NCONNECTIONS
   INTEGER :: NCELLS_MAX, NCELLS_TOTAL, NCELLS_START
   INTEGER :: NCONN_ACCUM, NCONN_MAX, NCONN_TOTAL, NCONN_START
   INTEGER :: NOFFSETS_ACCUM, NOFFSETS_MAX, NOFFSETS_TOTAL, NOFFSETS_START
   INTEGER :: NPOINTS_MAX, NPOINTS_TOTAL, NPOINTS_START
   INTEGER, DIMENSION(1:N_MPI_PROCESSES) :: MESHES_PER_PROCESS
   INTEGER :: N_WRITTEN, I, IPROC, IERR
   TYPE(PATCH_TYPE), POINTER :: PA
   INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCATIONS,FACES,SURFIND,GEOMIND
   REAL(FB), ALLOCATABLE, DIMENSION(:) :: VERTS
   REAL(FB), ALLOCATABLE, DIMENSION(:) :: X_PTS, Y_PTS, Z_PTS
   REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: VERTICES, ALL_VERTICES
   INTEGER(IB32) :: LAST_OFFSET_VALUE
   INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS, ALL_CONNECT, ALL_OFFSETS
   INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: VTKC_TYPE, ALL_VTKC_TYPE
   INTEGER(HSIZE_T), DIMENSION(1) :: START1(1), EXTENT1(1), START2(2), EXTENT2(2)
   INTEGER(IB32), DIMENSION(1) :: IDATA_OUT1(1)
   INTEGER(IB8), DIMENSION(1) :: DATA_OUT1_U8(1)
   REAL(FB), DIMENSION(1) :: FDATA_OUT3(3)
   INTEGER :: NINFO(6)
   TYPE (MPI_STATUS) :: MPISTATUS
   TYPE (BOUNDARY_FILE_TYPE), POINTER :: BF=>NULL()
   INTEGER(HSIZE_T), DIMENSION(1) :: MDIM

   DO NM=1,N_MPI_PROCESSES
      MESHES_PER_PROCESS(NM) = 0
   ENDDO
   DO NM=1,NMESHES
      MESHES_PER_PROCESS(PROCESS(NM)+1) = MESHES_PER_PROCESS(PROCESS(NM)+1) + 1
   ENDDO

   ! Initialize HDF5 datafiles
   CALL CREATE_OPEN_VTKHDF_SERIES(FILENAME,HDF_BNDF_FILE_ID,HDF_BNDF_PLIST_ID,&
      HDF_BNDF_G1,HDF_BNDF_G2,HDF_BNDF_G3,HDF_BNDF_G4,HDF_BNDF_G5,HDF_BNDF_G6,HDF_BNDF_G7)
   
   ! 2*NMESHES arrays. odds include OBST patches, evens include GEOM patches
   START1 = 2*NMESHES
   EXTENT1 = 2*NMESHES
   CALL PARALLEL_INIT_I32(HDF_BNDF_G1, TRIM("NumberOfCells"), HDF_BNDF_CRP_LIST, 1,&
      START1, EXTENT1, DSET_ID_NCELLS, HDF_BNDF_PLIST_ID) ! NumberOfCells
   CALL PARALLEL_INIT_I32(HDF_BNDF_G1, TRIM("NumberOfPoints"), HDF_BNDF_CRP_LIST, 1,&
      START1, EXTENT1, DSET_ID_NPTS, HDF_BNDF_PLIST_ID) ! NumberOfPoints
   CALL PARALLEL_INIT_I32(HDF_BNDF_G1, TRIM("NumberOfConnectivityIds"), HDF_BNDF_CRP_LIST, 1,&
      START1, EXTENT1, DSET_ID_NCON, HDF_BNDF_PLIST_ID) ! NumberOfConnectivityIds
   
   ! Fill metadata
   N_WRITTEN=0
   MESH_LOOP: DO NM=1,NMESHES
      CALL POINT_TO_MESH(NM)
      IF (PROCESS(NM)/=MY_RANK) CYCLE MESH_LOOP
      NCELLS = 0
      NPOINTS = 0
      NCELLS_ACCUM = 0
      NPOINTS_ACCUM = 0
      NCONNECTIONS_ACCUM = 0
      ! Count OBST patch info
      IF (MESHES(NM)%N_PATCH>0) THEN
         PATCH_LOOP1: DO IP=1,N_PATCH
            PA => PATCH(IP)
            IF (PA%OBST_INDEX<=0) CYCLE PATCH_LOOP1
            ! Initialize piece
            CALL BUILD_VTK_SOLID_PHASE_GEOMETRY(NM, PA, NCELLS, NPOINTS,&
               X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)
            NCONNECTIONS_ACCUM = NCONNECTIONS_ACCUM + NCELLS*4
            NCELLS_ACCUM = NCELLS_ACCUM + NCELLS
            NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS
            CALL DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)
         ENDDO PATCH_LOOP1
      ENDIF
      
      ! Write number of cells data to file
      START1 = 2*(NM-1)
      EXTENT1 = 1
      IDATA_OUT1 = NCELLS_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      ! Write number of points data to file
      START1 = 2*(NM-1)
      EXTENT1 = 1
      IDATA_OUT1 = NPOINTS_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      ! Write NumberOfConnectivityIds data to file
      START1 = 2*(NM-1)
      EXTENT1 = 1
      IDATA_OUT1 = NCONNECTIONS_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      ! Count GEOM patch info
      NCELLS = 0
      NPOINTS = 0
      NCELLS_ACCUM = 0
      NPOINTS_ACCUM = 0
      NCONNECTIONS_ACCUM = 0
      IF (MESHES(NM)%N_INTERNAL_CFACE_CELLS>0) THEN
         CALL GET_GEOMSIZES('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS)
         IF (NVERTS>0 .AND. NFACES>0) THEN
            ALLOCATE(VERTS(3*NVERTS))
            ALLOCATE(FACES(3*NFACES))
            ALLOCATE(LOCATIONS(NFACES))
            ALLOCATE(SURFIND(NFACES))
            ALLOCATE(GEOMIND(NFACES))
            CALL GET_GEOMINFO('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS,&
                              VERTS,FACES,LOCATIONS,SURFIND=SURFIND,GEOMIND=GEOMIND)
            CALL BUILD_VTK_GEOM_GEOMETRY(VERTS, FACES, NCELLS, NPOINTS, X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)
            DEALLOCATE(VERTS)
            DEALLOCATE(FACES)
            DEALLOCATE(LOCATIONS)
            DEALLOCATE(SURFIND)
            DEALLOCATE(GEOMIND)
            
            CALL DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)
            
            NCONNECTIONS_ACCUM = NCONNECTIONS_ACCUM + NCELLS*3
            NCELLS_ACCUM = NCELLS_ACCUM + NCELLS
            NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS
         ENDIF
      ENDIF
      
      ! Write number of cells data to file
      START1 = 2*(NM-1)+1
      EXTENT1 = 1
      IDATA_OUT1 = NCELLS_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      ! Write number of points data to file
      START1 = 2*(NM-1)+1
      EXTENT1 = 1
      IDATA_OUT1 = NPOINTS_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      ! Write NumberOfConnectivityIds data to file
      START1 = 2*(NM-1)+1
      EXTENT1 = 1
      IDATA_OUT1 = NCONNECTIONS_ACCUM
      CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
      N_WRITTEN = N_WRITTEN + 1
   ENDDO MESH_LOOP
   
   ! Write fake data in processes that have less meshes than max process
   DO NM=1,MAXVAL(MESHES_PER_PROCESS)
      IF (N_WRITTEN < MAXVAL(MESHES_PER_PROCESS)) THEN
         NCELLS = 0
         NPOINTS = 0
         NCONNECTIONS = 0
         START1 = NM-1
         EXTENT1 = 0
         IDATA_OUT1 = 0
         DO I=1,2 ! Loop through 2x to write empty twice per mesh
            CALL PARALLEL_WRITE_I32(1, DSET_ID_NCELLS, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
         
            ! Write number of points data to file
            CALL PARALLEL_WRITE_I32(1, DSET_ID_NPTS, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
         
            ! Write NumberOfConnectivityIds data to file
            CALL PARALLEL_WRITE_I32(1, DSET_ID_NCON, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
         ENDDO
            
         N_WRITTEN = N_WRITTEN + 1
      ENDIF
   ENDDO
   
   ! Close VTKHDF interface
   CALL VTK_DSET_RELEASE(DSET_ID_NCELLS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NPTS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NCON, ERROR)
   !CALL H5PCLOSE_F(PLIST_ID, ERROR)
   
   ! Exchange boundary NCELLS, NPOINTS, and NCONNECTIONS
   MESH_LOOP2: DO NM=1,NMESHES
      NM1 = 2*NM-1
      NM2 = 2*NM
      IF (PROCESS(NM)/=MY_RANK) THEN
         CALL MPI_RECV(NINFO,6,MPI_INTEGER,PROCESS(NM),PROCESS(NM),MPI_COMM_WORLD,MPISTATUS,IERR)
         NCELLS_VTK(NM1) = INT(NINFO(1),IB32)
         NPOINTS_VTK(NM1) = INT(NINFO(2),IB32)
         NCONNECTIONS_VTK(NM1) = INT(NINFO(3),IB32)
         NCELLS_VTK(NM2) = INT(NINFO(4),IB32)
         NPOINTS_VTK(NM2) = INT(NINFO(5),IB32)
         NCONNECTIONS_VTK(NM2) = INT(NINFO(6),IB32)
         CYCLE MESH_LOOP2
      ELSE
         NINFO = INT((/NCELLS_VTK(NM1),NPOINTS_VTK(NM1),NCONNECTIONS_VTK(NM1),&
                      NCELLS_VTK(NM2),NPOINTS_VTK(NM2),NCONNECTIONS_VTK(NM2)/))
         DO IPROC=0,N_MPI_PROCESSES-1
            IF (MY_RANK/=IPROC) THEN
               CALL MPI_SEND(NINFO,6,MPI_INTEGER,IPROC,MY_RANK,MPI_COMM_WORLD,IERR)
            ENDIF
         ENDDO
      ENDIF
   ENDDO MESH_LOOP2


   NPOINTS_TOTAL=0
   NCELLS_TOTAL=0
   NCELLS_ACCUM=0
   NPOINTS_ACCUM=0
   NOFFSETS_TOTAL=0
   NOFFSETS_ACCUM = 0
   NCONN_TOTAL = 0
   NCONN_ACCUM = 0
   NCELLS_MAX = 1
   NPOINTS_MAX = 1
   NCONN_MAX = 1
   NOFFSETS_MAX = 1
   NPOINTS_START = 0
   
   MESH_LOOP_HDF_COUNT: DO NMNM=1,2*NMESHES
      NPOINTS_TOTAL = NPOINTS_TOTAL + NPOINTS_VTK(NMNM) !+ NVERTS*3
      NCELLS_TOTAL = NCELLS_TOTAL + NCELLS_VTK(NMNM) !+ NFACES*3
      NOFFSETS_TOTAL = NOFFSETS_TOTAL + NCELLS_VTK(NMNM) + 1
      NCONN_TOTAL = NCONN_TOTAL + NCONNECTIONS_VTK(NMNM)
      NCELLS_MAX = MAX(NCELLS_MAX,NCELLS_VTK(NMNM))
      NPOINTS_MAX = MAX(NPOINTS_MAX,NPOINTS_VTK(NMNM))
      NOFFSETS_MAX = MAX(NOFFSETS_MAX,NCELLS_VTK(NMNM) + 1)
      NCONN_MAX = MAX(NCONN_MAX,NCONNECTIONS_VTK(NMNM))
   ENDDO MESH_LOOP_HDF_COUNT
   
   CALL POINT_TO_MESH(1)
   MDIM = (/H5S_UNLIMITED_F/)
   DO N=1,N_BNDF
      BF => BOUNDARY_FILE(N)
      IF (BF%CELL_CENTERED) THEN
         START1 = NCELLS_MAX
         EXTENT1 = (/0_HSIZE_T/)
         CALL PARALLEL_INIT_F32(HDF_BNDF_G2, BF%SMOKEVIEW_LABEL(1:30), HDF_BNDF_CRP_LIST, 1,&
            START1,EXTENT1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM) ! Data
      ELSE
         START1 = NPOINTS_MAX
         EXTENT1 = (/0_HSIZE_T/)
         CALL PARALLEL_INIT_F32(HDF_BNDF_G4, BF%SMOKEVIEW_LABEL(1:30), HDF_BNDF_CRP_LIST, 1,&
            START1,EXTENT1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM) ! Data
      ENDIF
      CALL VTK_DSET_RELEASE(DSET_ID, ERROR)
   ENDDO
   
   START1 = NCONN_MAX
   EXTENT1 = NCONN_TOTAL
   CALL PARALLEL_INIT_I32(HDF_BNDF_G1, TRIM("Connectivity"), HDF_BNDF_CRP_LIST, 1,&
      START1,EXTENT1, DSET_ID_CON, HDF_BNDF_PLIST_ID) ! Connectivity
   START1 = NOFFSETS_MAX
   EXTENT1 = NOFFSETS_TOTAL
   CALL PARALLEL_INIT_I32(HDF_BNDF_G1, TRIM("Offsets"), HDF_BNDF_CRP_LIST, 1,&
      START1,EXTENT1, DSET_ID_OFF, HDF_BNDF_PLIST_ID) ! Offsets
   START2(1:2) = INT((/3, NPOINTS_MAX/),HSIZE_T)
   EXTENT2(1:2) = INT((/3, NPOINTS_TOTAL/),HSIZE_T)
   CALL PARALLEL_INIT_F32(HDF_BNDF_G1, TRIM("Points"), HDF_BNDF_CRP_LIST, 2,&
      START2,EXTENT2, DSET_ID_PTS, HDF_BNDF_PLIST_ID) ! Points
   START1 = NCELLS_MAX
   EXTENT1 = NCELLS_TOTAL
   CALL PARALLEL_INIT_U8(HDF_BNDF_G1, TRIM("Types"), HDF_BNDF_CRP_LIST, 1,&
      START1,EXTENT1, DSET_ID_TYP, HDF_BNDF_PLIST_ID) ! Types
   
   ! Fill metadata
   N_WRITTEN=0
   MESH_LOOP_HDF: DO NM=1,NMESHES
      NM1 = 2*NM-1
      NM2 = 2*NM
      NPOINTS_START = NPOINTS_ACCUM
      NCELLS_START = NCELLS_ACCUM
      NCONN_START = NCONN_ACCUM
      NOFFSETS_START = NOFFSETS_ACCUM
      LAST_OFFSET_VALUE = 0
      IF (PROCESS(NM)/=MY_RANK) THEN
         NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS_VTK(NM1) + NPOINTS_VTK(NM2)
         NCELLS_ACCUM = NCELLS_ACCUM + NCELLS_VTK(NM1) + NCELLS_VTK(NM2)
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS_VTK(NM1)+1
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS_VTK(NM2)+1
         NCONN_ACCUM = NCONN_ACCUM + NCONNECTIONS_VTK(NM1) + NCONNECTIONS_VTK(NM2)
         CYCLE MESH_LOOP_HDF
      ENDIF
      CALL POINT_TO_MESH(NM)
      
      ! Build OBST boundary geometry
      IF (MESHES(NM)%N_PATCH>0) THEN
         ALLOCATE(ALL_VERTICES(3,NPOINTS_VTK(NM1)))
         ALLOCATE(ALL_CONNECT(NCONNECTIONS_VTK(NM1)))
         ALLOCATE(ALL_OFFSETS(NCELLS_VTK(NM1)+1))
         ALLOCATE(ALL_VTKC_TYPE(NCELLS_VTK(NM1)))
         ALL_OFFSETS(NOFFSETS_ACCUM-NOFFSETS_START+1) = 0
         PATCH_LOOP2: DO IP=1,N_PATCH
            PA => PATCH(IP)
            IF (PA%OBST_INDEX<=0) CYCLE PATCH_LOOP2
            ! Initialize piece
            CALL BUILD_VTK_SOLID_PHASE_GEOMETRY(NM, PA, PA_NCELLS, PA_NPOINTS,&
               X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)
            ALLOCATE(VERTICES(3,PA_NPOINTS))
            DO II=1,PA_NPOINTS
               VERTICES(1:3,II) = (/X_PTS(II),Y_PTS(II),Z_PTS(II)/)
            ENDDO
            ALL_VERTICES(:,NPOINTS_ACCUM-NPOINTS_START+1:NPOINTS_ACCUM-NPOINTS_START+PA_NPOINTS) = VERTICES
            ALL_CONNECT(NCONN_ACCUM-NCONN_START+1:NCONN_ACCUM-NCONN_START+PA_NCELLS*4) = CONNECT + NPOINTS_ACCUM-NPOINTS_START
            OFFSETS = OFFSETS + LAST_OFFSET_VALUE
            LAST_OFFSET_VALUE = OFFSETS(SIZE(OFFSETS))
            ALL_OFFSETS(NOFFSETS_ACCUM-NOFFSETS_START+2:NOFFSETS_ACCUM-NOFFSETS_START+PA_NCELLS+1) = OFFSETS
            ALL_VTKC_TYPE(NCELLS_ACCUM-NCELLS_START+1:NCELLS_ACCUM-NCELLS_START+PA_NCELLS) = VTKC_TYPE
            NCONN_ACCUM = NCONN_ACCUM + PA_NCELLS*4
            NCELLS_ACCUM = NCELLS_ACCUM + PA_NCELLS
            NPOINTS_ACCUM = NPOINTS_ACCUM + PA_NPOINTS
            NOFFSETS_ACCUM = NOFFSETS_ACCUM + PA_NCELLS
            CALL DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)
            DEALLOCATE(VERTICES)
         ENDDO PATCH_LOOP2
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + 1
         
         ! Write connectivity data to file
         START1 = NCONN_START
         EXTENT1 = NCONNECTIONS_VTK(NM1)
         CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, HDF_BNDF_PLIST_ID, START1, EXTENT1, ALL_CONNECT)
            
         ! Write offsets data to file
         START1 = NOFFSETS_START
         EXTENT1 = NCELLS_VTK(NM1)+1
         CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, HDF_BNDF_PLIST_ID, START1, EXTENT1, ALL_OFFSETS)
         
         ! Write point data to file
         START2(1:2) = INT((/0, NPOINTS_START/),HSIZE_T)
         EXTENT2(1:2) = INT((/3, NPOINTS_VTK(NM1)/),HSIZE_T)
         CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, HDF_BNDF_PLIST_ID, START2, EXTENT2, ALL_VERTICES)
         
         ! Write types data to file
         START1 = NCELLS_START
         EXTENT1 = NCELLS_VTK(NM1)
         CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, HDF_BNDF_PLIST_ID, START1, EXTENT1, ALL_VTKC_TYPE)
         
         DEALLOCATE(ALL_VERTICES)
         DEALLOCATE(ALL_CONNECT)
         DEALLOCATE(ALL_OFFSETS)
         DEALLOCATE(ALL_VTKC_TYPE)
         N_WRITTEN = N_WRITTEN + 1
      ELSE
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + 1
      ENDIF
         
      ! Build GEOM patch geometry
      IF ((MESHES(NM)%N_INTERNAL_CFACE_CELLS>0).AND.(.TRUE.)) THEN
         NPOINTS_START = NPOINTS_ACCUM
         NCELLS_START = NCELLS_ACCUM
         NCONN_START = NCONN_ACCUM
         NOFFSETS_START = NOFFSETS_ACCUM
         CALL GET_GEOMSIZES('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS)
         IF (NVERTS>0 .AND. NFACES>0) THEN
            ALLOCATE(ALL_OFFSETS(NCELLS_VTK(NM2)+1))
            ALLOCATE(VERTS(3*NVERTS))
            ALLOCATE(FACES(3*NFACES))
            ALLOCATE(LOCATIONS(NFACES))
            ALLOCATE(SURFIND(NFACES))
            ALLOCATE(GEOMIND(NFACES))
            CALL GET_GEOMINFO('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS,&
                              VERTS,FACES,LOCATIONS,SURFIND=SURFIND,GEOMIND=GEOMIND)
            CALL BUILD_VTK_GEOM_GEOMETRY(VERTS, FACES, PA_NCELLS, PA_NPOINTS, X_PTS, Y_PTS, Z_PTS, CONNECT, OFFSETS, VTKC_TYPE)
            ALLOCATE(VERTICES(3,PA_NPOINTS))
            DO II=1,PA_NPOINTS
               VERTICES(1:3,II) = (/X_PTS(II),Y_PTS(II),Z_PTS(II)/)
            ENDDO
            ALL_OFFSETS(1) = 0
            ALL_OFFSETS(2:NCELLS_VTK(NM2)+1) = OFFSETS
            NCONN_ACCUM = NCONN_ACCUM + PA_NCELLS*3
            NCELLS_ACCUM = NCELLS_ACCUM + PA_NCELLS
            NPOINTS_ACCUM = NPOINTS_ACCUM + PA_NPOINTS
            NOFFSETS_ACCUM = NOFFSETS_ACCUM + PA_NCELLS + 1
            
            ! Write connectivity data to file
            START1 = NCONN_START
            EXTENT1 = NCONNECTIONS_VTK(NM2)
            CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, HDF_BNDF_PLIST_ID, START1, EXTENT1, CONNECT)
               
            ! Write offsets data to file
            START1 = NOFFSETS_START+1
            EXTENT1 = NCELLS_VTK(NM2)+1
            CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, HDF_BNDF_PLIST_ID, START1, EXTENT1, ALL_OFFSETS)
            
            ! Write point data to file
            START2(1:2) = INT((/0, NPOINTS_START/),HSIZE_T)
            EXTENT2(1:2) = INT((/3, NPOINTS_VTK(NM2)/),HSIZE_T)
            CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, HDF_BNDF_PLIST_ID, START2, EXTENT2, VERTICES)
            
            ! Write types data to file
            START1 = NCELLS_START
            EXTENT1 = NCELLS_VTK(NM2)
            CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, HDF_BNDF_PLIST_ID, START1, EXTENT1, VTKC_TYPE)
            
            DEALLOCATE(VERTS)
            DEALLOCATE(FACES)
            DEALLOCATE(LOCATIONS)
            DEALLOCATE(SURFIND)
            DEALLOCATE(GEOMIND)
            DEALLOCATE(VERTICES)
            DEALLOCATE(ALL_OFFSETS)
            CALL DEALLOCATE_VTK_GAS_PHASE_GEOMETRY(X_PTS,Y_PTS,Z_PTS,OFFSETS,VTKC_TYPE,CONNECT)
            N_WRITTEN = N_WRITTEN + 1
         ENDIF
      ELSE
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + 1
      ENDIF
   ENDDO MESH_LOOP_HDF
   
   ! Write fake data in processes that have less meshes than max process
   DO NM=1,MAXVAL(MESHES_PER_PROCESS)*2
      IF (N_WRITTEN < MAXVAL(MESHES_PER_PROCESS)*2) THEN
   
         ! Write connectivity data to file
         START1 = 0
         EXTENT1 = 0
         IDATA_OUT1 = 0
         CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
         
         ! Write offsets data to file
         CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, HDF_BNDF_PLIST_ID, START1, EXTENT1, IDATA_OUT1)
      
         ! Write point data to file
         START2(1:2) = INT((/0, 0/),HSIZE_T)
         EXTENT2(1:2) = INT((/3, 0/),HSIZE_T)
         FDATA_OUT3 = 0
         CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, HDF_BNDF_PLIST_ID, START2, EXTENT2, FDATA_OUT3)
         
         ! Write types data to file
         DATA_OUT1_U8 = 0
         CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, HDF_BNDF_PLIST_ID, START1, EXTENT1, DATA_OUT1_U8)
         N_WRITTEN = N_WRITTEN + 1
      ENDIF
   ENDDO
      
   ! Close VTKHDF interface
   CALL VTK_DSET_RELEASE(DSET_ID_CON, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_OFF, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_PTS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_TYP, ERROR)

END SUBROUTINE INITIALIZE_VTKHDF_BNDF


























!> \brief Write every locally owned mesh's slice data in a single collective operation
!>
!> \param DATASET  Name of the point-data dataset
!> \param GROUP_ID4 Group holding the point data
!> \param NPOINTS  Number of points contributed by each mesh, global mesh order
!> \param BASE_OFFSET Offset of this time step within the (appended) dataset
!> \param N_LOCAL  Number of locally owned meshes that contributed data
!> \param WRITE_NM Global mesh indices of those meshes, in increasing order
!> \param DATA     Their data, concatenated in the same order
!>
!> Every rank calls this exactly once per dataset per time step, so the number of
!> collective HDF5 operations is independent of the number of meshes.

SUBROUTINE WRITE_VTKHDF_SLICE_DATA_MULTI(DATASET,GROUP_ID4,CRP_LIST,PLIST_ID,NPOINTS,BASE_OFFSET,&
                                         N_LOCAL,WRITE_NM,DATA)
   CHARACTER(*), INTENT(IN) :: DATASET
   INTEGER(HID_T), INTENT(IN) :: GROUP_ID4
   INTEGER(HID_T), INTENT(INOUT) :: CRP_LIST, PLIST_ID
   INTEGER(IB32), DIMENSION(1:NMESHES), INTENT(IN) :: NPOINTS
   INTEGER(HSIZE_T), INTENT(IN) :: BASE_OFFSET
   INTEGER, INTENT(IN) :: N_LOCAL
   INTEGER, DIMENSION(:), INTENT(IN) :: WRITE_NM
   REAL(FB), DIMENSION(*), INTENT(IN) :: DATA
   INTEGER(HID_T) :: DSET_ID, MEMSPACE, DATASPACE
   INTEGER(HSIZE_T) :: OFFSET_ACCUM(NMESHES), NPOINTS_MAX, NPOINTS_TOTAL, NLOC
   INTEGER(HSIZE_T), DIMENSION(1) :: START1, EXTENT1
   INTEGER :: NM, I, ERROR

   ! Offset of each mesh's block within one time step, and the sizing used at creation

   NPOINTS_TOTAL = 0_HSIZE_T
   NPOINTS_MAX = 0_HSIZE_T
   DO NM=1,NMESHES
      OFFSET_ACCUM(NM) = NPOINTS_TOTAL
      NPOINTS_TOTAL = NPOINTS_TOTAL + INT(NPOINTS(NM),HSIZE_T)
      NPOINTS_MAX = MAX(NPOINTS_MAX,INT(NPOINTS(NM),HSIZE_T))
   ENDDO

   START1 = NPOINTS_MAX
   EXTENT1 = NPOINTS_TOTAL
   CALL PARALLEL_INIT_F32(GROUP_ID4, DATASET, CRP_LIST, 1, START1, EXTENT1, DSET_ID, PLIST_ID)

   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SSELECT_NONE_F(DATASPACE, ERROR)

   NLOC = 0_HSIZE_T
   DO I=1,N_LOCAL
      NM = WRITE_NM(I)
      START1  = OFFSET_ACCUM(NM) + BASE_OFFSET
      EXTENT1 = INT(NPOINTS(NM),HSIZE_T)
      IF (EXTENT1(1)==0_HSIZE_T) CYCLE
      CALL H5SSELECT_HYPERSLAB_F(DATASPACE, H5S_SELECT_OR_F, START1, EXTENT1, ERROR)
      NLOC = NLOC + EXTENT1(1)
   ENDDO

   EXTENT1 = NLOC
   CALL H5SCREATE_SIMPLE_F(1, EXTENT1, MEMSPACE, ERROR)
   IF (NLOC==0_HSIZE_T) CALL H5SSELECT_NONE_F(MEMSPACE, ERROR)
   CALL H5DWRITE_F(DSET_ID, H5T_IEEE_F32LE, DATA, EXTENT1, ERROR, &
                   MEM_SPACE_ID = MEMSPACE, FILE_SPACE_ID = DATASPACE, XFER_PRP = PLIST_ID)

   CALL H5SCLOSE_F(MEMSPACE, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID, ERROR)

END SUBROUTINE WRITE_VTKHDF_SLICE_DATA_MULTI
















!> \brief Write one mesh's contribution to a slice plane's grid
















!>
















!> \param SLCFNAME Name of the slice plane
















!> \param SL3D True if the slice spans three dimensions
















!> \param NM Mesh number
















!> \param NCELLS Number of cells each mesh contributes
















!> \param NPOINTS Number of points each mesh contributes
















!> \param NCONNECTIONS Number of connectivity entries each mesh contributes
















!> \param NTSL Terrain slice counter, used to index K_AGL_SLICE
















!> \param PLIST_ID Data transfer property list for this file
















!> \param CRP_LIST Dataset creation property list for this file
















!> \param GROUP_ID1 VTKHDF group of this file

















SUBROUTINE WRITE_VTKHDF_SLICE_CELL_FILE_NOOPEN(SLCFNAME,SL3D,NM,NCELLS,NPOINTS,NCONNECTIONS,NTSL,&
   PLIST_ID,CRP_LIST,GROUP_ID1)
   CHARACTER(*), INTENT(IN) :: SLCFNAME
   INTEGER, INTENT(IN) :: NM, NTSL
   LOGICAL, INTENT(IN) :: SL3D
   INTEGER(HID_T), INTENT(IN) :: GROUP_ID1  ! Identifiers
   INTEGER(HID_T), INTENT(INOUT) :: PLIST_ID  ! Identifiers
   INTEGER(HID_T), INTENT(INOUT) :: CRP_LIST       ! Identifiers
   INTEGER(HID_T) :: DSET_ID_CON, DSET_ID_NCELLS    ! Dataset identifiers
   INTEGER(HID_T) :: DSET_ID_NPTS, DSET_ID_OFF, DSET_ID_PTS, DSET_ID_TYP       ! Dataset identifiers
   INTEGER :: NPOINTS_TOTAL, NCELLS_TOTAL, NOFFSETS_TOTAL
   INTEGER :: NCONN_TOTAL, NPIECES
   INTEGER :: NCONN_ACCUM, NPIECES_ACCUM
   INTEGER :: NCONN_MAX, NOFFSETS_MAX, NCELLS_MAX, NPOINTS_MAX
   INTEGER :: NMNM, ERROR, NC, NP, IQ, NQT
   TYPE (MESH_TYPE), POINTER :: M
   REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: VERTICES
   INTEGER(IB32), ALLOCATABLE, DIMENSION(:) :: CONNECT, OFFSETS
   INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: VTKC_TYPE
   INTEGER(IB32), DIMENSION(1:NMESHES), INTENT(OUT) :: NCELLS, NPOINTS
   INTEGER(IB32) :: NCELLS_ACCUM, NPOINTS_ACCUM, NOFFSETS_ACCUM
   TYPE(SLICE_TYPE), POINTER :: SL
   INTEGER, INTENT(OUT) :: NCONNECTIONS
   INTEGER(HSIZE_T), DIMENSION(1) :: START1(1), EXTENT1(1), START2(2), EXTENT2(2)
   INTEGER(HSIZE_T), DIMENSION(1) :: START3(1), EXTENT3(1), START4(1), EXTENT4(4)
   
   ! Read number of cells
   CALL H5DOPEN_F(GROUP_ID1, "NumberOfCells", DSET_ID_NCELLS, ERROR)
   START1 = NMESHES
   CALL H5DREAD_F(DSET_ID_NCELLS, H5T_STD_I32LE, NCELLS, START1, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NCELLS, ERROR)

   ! Read number of points
   CALL H5DOPEN_F(GROUP_ID1, "NumberOfPoints", DSET_ID_NPTS, ERROR)
   CALL H5DREAD_F(DSET_ID_NPTS, H5T_STD_I32LE, NPOINTS, START1, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_NPTS, ERROR)
   
   IF (SL3D) THEN
      NCONNECTIONS=8 ! 3-D slice
   ELSE
      NCONNECTIONS=4 ! 2-D slice
   ENDIF
   
   NPOINTS_TOTAL=0
   NCELLS_TOTAL=0
   NCELLS_ACCUM=0
   NPOINTS_ACCUM=0
   NOFFSETS_TOTAL=0
   NOFFSETS_ACCUM = 0
   NCONN_TOTAL = 0
   NPIECES = 0
   NCONN_ACCUM = 0
   NPIECES_ACCUM = 0
   NCELLS_MAX = 0
   NPOINTS_MAX = 0
   NCONN_MAX = 0
   NOFFSETS_MAX = 0
   
   MESH_LOOP_HDF_COUNT: DO NMNM=1,NMESHES
      NPOINTS_TOTAL = NPOINTS_TOTAL + NPOINTS(NMNM) !+ NVERTS*3
      NCELLS_TOTAL = NCELLS_TOTAL + NCELLS(NMNM) !+ NFACES*3
      NOFFSETS_TOTAL = NOFFSETS_TOTAL + NCELLS(NMNM) + 1
      NCONN_TOTAL = NCONN_TOTAL + NCELLS(NMNM)*NCONNECTIONS
      NPIECES = NPIECES + 1
      NCELLS_MAX = MAX(NCELLS_MAX,NCELLS(NMNM))
      NPOINTS_MAX = MAX(NPOINTS_MAX,NPOINTS(NMNM))
      NOFFSETS_MAX = MAX(NOFFSETS_MAX,NCELLS(NMNM) + 1)
      NCONN_MAX = MAX(NCONN_MAX,NCELLS(NMNM)*NCONNECTIONS)
   ENDDO MESH_LOOP_HDF_COUNT
   
   START1 = NCONN_MAX
   EXTENT1 = NCONN_TOTAL
   CALL PARALLEL_INIT_I32(GROUP_ID1, TRIM("Connectivity"), CRP_LIST, 1,&
      START1, EXTENT1, DSET_ID_CON, PLIST_ID) ! Connectivity
   START1 = NOFFSETS_MAX
   EXTENT1 = NOFFSETS_TOTAL
   CALL PARALLEL_INIT_I32(GROUP_ID1, TRIM("Offsets"), CRP_LIST, 1,&
      START1, EXTENT1, DSET_ID_OFF, PLIST_ID) ! Offsets
   START2(1:2) = INT((/3, NPOINTS_MAX/),HSIZE_T)
   EXTENT2(1:2) = INT((/3, NPOINTS_TOTAL/),HSIZE_T)
   CALL PARALLEL_INIT_F32(GROUP_ID1, TRIM("Points"), CRP_LIST, 2,&
      START2, EXTENT2, DSET_ID_PTS, PLIST_ID) ! Points
   START1 = NCELLS_MAX
   EXTENT1 = NCELLS_TOTAL
   CALL PARALLEL_INIT_U8(GROUP_ID1, TRIM("Types"), CRP_LIST, 1,&
      START1, EXTENT1, DSET_ID_TYP, PLIST_ID) ! Types
   
   CALL POINT_TO_MESH(NM)
   NQT = MESHES(1)%N_SLCF_VTK
   MESH_LOOP_HDF: DO NMNM=1,NMESHES
      IF (NMNM/=NM) THEN
         NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS(NMNM)
         NCELLS_ACCUM = NCELLS_ACCUM + NCELLS(NMNM)
         !IF (NCELLS(NMNM) > 0) NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS(NMNM) + 1
         NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS(NMNM) + 1
         NCONN_ACCUM = NCONN_ACCUM + NCELLS(NMNM)*NCONNECTIONS
         NPIECES_ACCUM = NPIECES_ACCUM+1
         CYCLE MESH_LOOP_HDF
      ENDIF
      M => MESHES(NMNM)
      
      IF (PROCESS(NMNM)/=MY_RANK) THEN
         ALLOCATE(VERTICES(3,0))
         ALLOCATE(OFFSETS(0))
         ALLOCATE(VTKC_TYPE(0))
         ALLOCATE(CONNECT(0))
         START1 = NCONN_ACCUM
         EXTENT1 = 0
         START2(1:2) = INT((/0,NPOINTS_ACCUM/),HSIZE_T)
         EXTENT2(1:2) = INT((/3, 0/),HSIZE_T)
         START3 = NOFFSETS_ACCUM
         EXTENT3 = 0
         START4 = NCELLS_ACCUM
         EXTENT4 = 0
      ELSEIF (NCELLS(NMNM)==0) THEN
         ALLOCATE(VERTICES(3,NCELLS(NMNM)))
         ALLOCATE(OFFSETS(NCELLS(NMNM)))
         ALLOCATE(VTKC_TYPE(NCELLS(NMNM)))
         ALLOCATE(CONNECT(NCELLS(NMNM)*NCONNECTIONS))
         START1 = NCONN_ACCUM
         EXTENT1 = NCELLS(NMNM)*NCONNECTIONS
         START2(1:2) = INT((/0,NPOINTS_ACCUM/),HSIZE_T)
         EXTENT2(1:2) = INT((/3, NPOINTS(NMNM)/),HSIZE_T)
         START3 = NOFFSETS_ACCUM
         EXTENT3 = 0 !NCELLS(NMNM)+1
         START4 = NCELLS_ACCUM
         EXTENT4 = 0 !NCELLS(NMNM)
         VERTICES = 0
         OFFSETS = 0
         VTKC_TYPE = 0
         CONNECT = 0
      ELSE
         QUANTITY_LOOPB: DO IQ=1,NQT
            SL => SLICE(IQ)
            IF (TRIM(SL%SLCF_NAME)/=TRIM(SLCFNAME)) CYCLE QUANTITY_LOOPB
            CALL BUILD_VTK_SLICE_GEOMETRY2(NMNM, SL, NTSL, NC, NP, VERTICES, CONNECT, OFFSETS, VTKC_TYPE)
            EXIT
         ENDDO QUANTITY_LOOPB
         START1 = NCONN_ACCUM
         EXTENT1 = NCELLS(NMNM)*NCONNECTIONS
         START2(1:2) = INT((/0,NPOINTS_ACCUM/),HSIZE_T)
         EXTENT2(1:2) = INT((/3, NPOINTS(NMNM)/),HSIZE_T)
         START3 = NOFFSETS_ACCUM
         EXTENT3 = NCELLS(NMNM)+1
         START4 = NCELLS_ACCUM
         EXTENT4 = NCELLS(NMNM)
      ENDIF
      
      ! Write connectivity data to file
      CALL PARALLEL_WRITE_I32(1, DSET_ID_CON, PLIST_ID, START1, EXTENT1, CONNECT)
      
      ! Write point data to file
      CALL PARALLEL_WRITE_F32(2, DSET_ID_PTS, PLIST_ID, START2, EXTENT2, VERTICES)
      
      ! Write offsets data to file
      CALL PARALLEL_WRITE_I32(1, DSET_ID_OFF, PLIST_ID, START3, EXTENT3, OFFSETS)
      
      ! Write types data to file
      CALL PARALLEL_WRITE_U8(1, DSET_ID_TYP, PLIST_ID, START4, EXTENT4, VTKC_TYPE)
      
      DEALLOCATE(VERTICES)
      DEALLOCATE(OFFSETS)
      DEALLOCATE(CONNECT)
      DEALLOCATE(VTKC_TYPE)
      NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS(NMNM)
      NCELLS_ACCUM = NCELLS_ACCUM + NCELLS(NMNM)
      !IF (NCELLS(NM)>0) NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS(NM) + 1
      NOFFSETS_ACCUM = NOFFSETS_ACCUM + NCELLS(NMNM) + 1
      NCONN_ACCUM = NCONN_ACCUM + NCELLS(NMNM)*NCONNECTIONS
      NPIECES_ACCUM = NPIECES_ACCUM+1
   ENDDO MESH_LOOP_HDF
   
   ! Close HDF5 datasets
   CALL VTK_DSET_RELEASE(DSET_ID_CON, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_OFF, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_PTS, ERROR)
   CALL VTK_DSET_RELEASE(DSET_ID_TYP, ERROR)
   !CALL H5PCLOSE_F(PLIST_ID, ERROR)
   
   
END SUBROUTINE WRITE_VTKHDF_SLICE_CELL_FILE_NOOPEN





!> \brief Return the shared collective data-transfer property list, creating it once

!> \brief Look up an already open dataset handle
!>
!> \param GROUP_ID Group the dataset lives in
!> \param SNAME Dataset name
!> \param DSET_ID Open handle, set only when the lookup succeeds

LOGICAL FUNCTION VTK_DSET_CACHE_FIND(GROUP_ID,SNAME,DSET_ID)
INTEGER(HID_T), INTENT(IN) :: GROUP_ID
CHARACTER(*), INTENT(IN) :: SNAME
INTEGER(HID_T), INTENT(OUT) :: DSET_ID
INTEGER :: I
DSET_ID = -1_HID_T
VTK_DSET_CACHE_FIND = .FALSE.
IF (LEN_TRIM(SNAME)>VTK_DSET_NAME_LEN) RETURN
DO I=1,VTK_DSET_CACHE_N
   IF (VTK_DSET_CACHE_GROUP(I)==GROUP_ID .AND. VTK_DSET_CACHE_NAME(I)==SNAME) THEN
      DSET_ID = VTK_DSET_CACHE_ID(I)
      VTK_DSET_CACHE_FIND = .TRUE.
      RETURN
   ENDIF
ENDDO
END FUNCTION VTK_DSET_CACHE_FIND


!> \brief Hand a newly opened dataset to the cache, which then owns it
!>
!> \param GROUP_ID Group the dataset lives in
!> \param SNAME Dataset name
!> \param DSET_ID Open handle

SUBROUTINE VTK_DSET_CACHE_ADD(GROUP_ID,SNAME,DSET_ID)
INTEGER(HID_T), INTENT(IN) :: GROUP_ID, DSET_ID
CHARACTER(*), INTENT(IN) :: SNAME
IF (LEN_TRIM(SNAME)>VTK_DSET_NAME_LEN) RETURN
IF (VTK_DSET_CACHE_N>=VTK_DSET_CACHE_MAX) RETURN
VTK_DSET_CACHE_N = VTK_DSET_CACHE_N + 1
VTK_DSET_CACHE_GROUP(VTK_DSET_CACHE_N) = GROUP_ID
VTK_DSET_CACHE_NAME(VTK_DSET_CACHE_N)  = SNAME
VTK_DSET_CACHE_ID(VTK_DSET_CACHE_N)    = DSET_ID
END SUBROUTINE VTK_DSET_CACHE_ADD


!> \brief Release a dataset handle obtained from PARALLEL_INIT_*
!>
!> \param DSET_ID Handle to release
!> \param ERROR HDF5 error flag
!>
!> Cached handles stay open until the file closes; anything else is closed here, so
!> call sites that open a dataset themselves still behave as they did.

SUBROUTINE VTK_DSET_RELEASE(DSET_ID,ERROR)
INTEGER(HID_T), INTENT(IN) :: DSET_ID
INTEGER, INTENT(OUT) :: ERROR
INTEGER :: I
ERROR = 0
DO I=1,VTK_DSET_CACHE_N
   IF (VTK_DSET_CACHE_ID(I)==DSET_ID) RETURN
ENDDO
CALL H5DCLOSE_F(DSET_ID, ERROR)
END SUBROUTINE VTK_DSET_RELEASE


!> \brief Close every cached dataset handle
!>
!> Called before a VTKHDF file closes.  The cache spans every VTKHDF file and they are
!> closed together, so there is nothing to gain from purging one file at a time.

SUBROUTINE VTK_DSET_CACHE_PURGE()
INTEGER :: I, ERROR
DO I=1,VTK_DSET_CACHE_N
   CALL H5DCLOSE_F(VTK_DSET_CACHE_ID(I), ERROR)
   VTK_DSET_CACHE_GROUP(I) = -1_HID_T
   VTK_DSET_CACHE_ID(I)    = -1_HID_T
   VTK_DSET_CACHE_NAME(I)  = ''
ENDDO
VTK_DSET_CACHE_N = 0
END SUBROUTINE VTK_DSET_CACHE_PURGE


!> \brief Is this one of the Steps bookkeeping groups?
!>
!> \param GROUP_ID Group to test
!>
!> The datasets under VTKHDF/Steps and VTKHDF/Steps/PointDataOffsets take one integer
!> per output time.  Running the filter pipeline over a whole chunk to add four bytes
!> costs more than the bytes are worth, and leaves the chunk dirty so every H5Fflush
!> compresses and writes it again.  Only called when a dataset is first created.

LOGICAL FUNCTION VTK_IS_STEPS_GROUP(GROUP_ID)
INTEGER(HID_T), INTENT(IN) :: GROUP_ID
CHARACTER(256) :: GNAME
INTEGER(SIZE_T) :: NAMELEN
INTEGER :: ERROR
GNAME = ''
CALL H5IGET_NAME_F(GROUP_ID, GNAME, INT(LEN(GNAME),SIZE_T), NAMELEN, ERROR)
VTK_IS_STEPS_GROUP = ERROR>=0 .AND. INDEX(GNAME,'/Steps')>0
END FUNCTION VTK_IS_STEPS_GROUP


!> \brief Return the data transfer property list shared by every VTKHDF write


!>


!> \param PLIST_ID Data transfer property list


!>


!> Created on first use and kept for the run.  Creating one per dataset access,


!> as this once did, leaked a property list on every write.



SUBROUTINE GET_VTK_DXPL(PLIST_ID)
   INTEGER(HID_T), INTENT(OUT) :: PLIST_ID
   INTEGER :: ERROR
   IF (VTK_DXPL_ID<0_HID_T) THEN
      CALL H5PCREATE_F(H5P_DATASET_XFER_F, VTK_DXPL_ID, ERROR)
      ! Collective transfer is meaningless, and rejected, without the MPI-IO driver
      IF (N_MPI_PROCESSES>1) CALL H5PSET_DXPL_MPIO_F(VTK_DXPL_ID, H5FD_MPIO_COLLECTIVE_F, ERROR)
   ENDIF
   PLIST_ID = VTK_DXPL_ID
END SUBROUTINE GET_VTK_DXPL


!> \brief Build the dataset creation property list shared by every VTKHDF dataset
!>
!> \param CRP_LIST Dataset creation property list, created here
!> \param RANK Dataset rank
!> \param CHUNK Chunk dimensions
!>
!> The shuffle filter transposes the bytes of each element so that deflate sees
!> runs of equal high-order bytes.  It costs almost nothing and is worth far more
!> than a higher deflate level on the geometry datasets, which are long runs of
!> slowly increasing integers: on an eight mesh 50x50x20 case it takes the
!> Connectivity dataset of the Smoke3D file from 3.8 MB to 0.2 MB.

SUBROUTINE GET_VTK_DCPL(CRP_LIST,RANK,CHUNK,TINY_CHUNK)
INTEGER(HID_T), INTENT(OUT) :: CRP_LIST
INTEGER, INTENT(IN) :: RANK
INTEGER(HSIZE_T), DIMENSION(RANK), INTENT(IN) :: CHUNK
LOGICAL, INTENT(IN) :: TINY_CHUNK
INTEGER :: ERROR
CALL H5PCREATE_F(H5P_DATASET_CREATE_F, CRP_LIST, ERROR)
IF (VTK_COMPRESSION_LEVEL>0 .AND. .NOT.TINY_CHUNK) THEN
   CALL H5PSET_SHUFFLE_F(CRP_LIST, ERROR)   ! must precede deflate in the filter pipeline
   CALL H5PSET_DEFLATE_F(CRP_LIST, VTK_COMPRESSION_LEVEL, ERROR)
ENDIF
CALL H5PSET_CHUNK_F(CRP_LIST, RANK, CHUNK, ERROR)
END SUBROUTINE GET_VTK_DCPL


!> \brief Clamp a requested chunk size into a range HDF5 handles efficiently


!>


!> \param RANK Dataset rank


!> \param DDIM Requested chunk dimensions


!> \param CHUNK Chunk dimensions to use (out)


!>


!> A chunk of one element costs a B-tree entry and a filter pipeline pass each;


!> a chunk of a whole mesh blows past the chunk cache.



PURE SUBROUTINE VTK_CHUNK_DIMS(RANK,DDIM,CHUNK)
   INTEGER, INTENT(IN) :: RANK
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: DDIM
   INTEGER(HSIZE_T), DIMENSION(RANK), INTENT(OUT) :: CHUNK
   INTEGER :: I
   DO I=1,RANK
      CHUNK(I) = DDIM(I)
   ENDDO
   CHUNK(RANK) = MIN(MAX(CHUNK(RANK),VTK_CHUNK_MIN),VTK_CHUNK_MAX)
END SUBROUTINE VTK_CHUNK_DIMS


!> \brief Open, or create on first use, a 32 bit real dataset


!>


!> \param GROUP_ID Group the dataset lives in


!> \param SNAME Dataset name


!> \param CRP_LIST Dataset creation property list, used when the dataset is created


!> \param RANK Dataset rank


!> \param DDIM Requested chunk dimensions


!> \param CDIM Initial dataset dimensions


!> \param DSET_ID Open dataset handle (out)


!> \param PLIST_ID Data transfer property list for the write that follows (out)


!> \param MDIM Maximum dimensions, present when the dataset grows over time


!> \param NOFILTER Create the dataset without the filter pipeline



SUBROUTINE PARALLEL_INIT_F32(GROUP_ID, SNAME, CRP_LIST, RANK, DDIM, CDIM, DSET_ID, PLIST_ID, MDIM, NOFILTER)
   INTEGER(HID_T), INTENT(IN) :: GROUP_ID               ! Memory identifiers
   INTEGER(HID_T), INTENT(INOUT) :: CRP_LIST, PLIST_ID  ! Memory identifiers
   INTEGER, INTENT(IN)     ::   RANK    ! Dataset rank
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: CDIM, DDIM
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN), OPTIONAL :: MDIM
   LOGICAL, INTENT(IN), OPTIONAL :: NOFILTER
   CHARACTER(LEN=*), INTENT(IN) :: SNAME
   INTEGER(HID_T), INTENT(OUT) :: DSET_ID  ! Memory identifiers
   INTEGER(HID_T) :: DATASPACE
   INTEGER(HSIZE_T), DIMENSION(RANK) :: CHUNK
   INTEGER     ::   ERROR ! Error flag
   LOGICAL :: LINK_EXISTS, TINY_CHUNK
   IF (VTK_DSET_CACHE_FIND(GROUP_ID, SNAME, DSET_ID)) THEN
      CALL GET_VTK_DXPL(PLIST_ID)
      RETURN
   ENDIF
   CALL H5LEXISTS_F(GROUP_ID, SNAME, LINK_EXISTS, ERROR)
   IF (.NOT.LINK_EXISTS) THEN
      TINY_CHUNK = VTK_IS_STEPS_GROUP(GROUP_ID)
      IF (PRESENT(NOFILTER)) TINY_CHUNK = TINY_CHUNK .OR. NOFILTER
      IF (PRESENT(MDIM)) THEN
         CALL H5SCREATE_SIMPLE_F(RANK, CDIM, DATASPACE, ERROR, MDIM)
         CALL VTK_CHUNK_DIMS(RANK, DDIM, CHUNK)
      ELSE
         CALL H5SCREATE_SIMPLE_F(RANK, CDIM, DATASPACE, ERROR)
         CHUNK(1:RANK) = DDIM(1:RANK)
      ENDIF
      CALL GET_VTK_DCPL(CRP_LIST, RANK, CHUNK, TINY_CHUNK)
      CALL H5DCREATE_F(GROUP_ID, SNAME, H5T_IEEE_F32LE, DATASPACE, DSET_ID, ERROR, CRP_LIST)
      CALL H5PCLOSE_F(CRP_LIST, ERROR)
      CALL GET_VTK_DXPL(PLIST_ID)
      CALL H5SCLOSE_F(DATASPACE, ERROR)
      CALL VTK_DSET_CACHE_ADD(GROUP_ID, SNAME, DSET_ID)
   ELSE
      CALL H5DOPEN_F(GROUP_ID, SNAME, DSET_ID, ERROR)
      CALL GET_VTK_DXPL(PLIST_ID)
      CALL VTK_DSET_CACHE_ADD(GROUP_ID, SNAME, DSET_ID)
   ENDIF
END SUBROUTINE PARALLEL_INIT_F32

!> \brief Open, or create on first use, a 32 bit integer dataset

!>

!> \param GROUP_ID Group the dataset lives in

!> \param SNAME Dataset name

!> \param CRP_LIST Dataset creation property list, used when the dataset is created

!> \param RANK Dataset rank

!> \param DDIM Requested chunk dimensions

!> \param CDIM Initial dataset dimensions

!> \param DSET_ID Open dataset handle (out)

!> \param PLIST_ID Data transfer property list for the write that follows (out)

!> \param MDIM Maximum dimensions, present when the dataset grows over time


SUBROUTINE PARALLEL_INIT_I32(GROUP_ID, SNAME, CRP_LIST, RANK, DDIM, CDIM, DSET_ID, PLIST_ID, MDIM)
   INTEGER(HID_T), INTENT(IN) :: GROUP_ID               ! Memory identifiers
   INTEGER(HID_T), INTENT(INOUT) :: CRP_LIST, PLIST_ID  ! Memory identifiers
   INTEGER, INTENT(IN)     ::   RANK    ! Dataset rank
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: CDIM, DDIM
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN), OPTIONAL :: MDIM
   CHARACTER(LEN=*), INTENT(IN) :: SNAME
   INTEGER(HID_T), INTENT(OUT) :: DSET_ID  ! Memory identifiers
   INTEGER(HID_T) :: DATASPACE
   INTEGER(HSIZE_T), DIMENSION(RANK) :: CHUNK
   INTEGER     ::   ERROR ! Error flag
   LOGICAL :: LINK_EXISTS, TINY_CHUNK
   IF (VTK_DSET_CACHE_FIND(GROUP_ID, SNAME, DSET_ID)) THEN
      CALL GET_VTK_DXPL(PLIST_ID)
      RETURN
   ENDIF
   CALL H5LEXISTS_F(GROUP_ID, SNAME, LINK_EXISTS, ERROR)
   IF (.NOT.LINK_EXISTS) THEN
      TINY_CHUNK = VTK_IS_STEPS_GROUP(GROUP_ID)
      IF (PRESENT(MDIM)) THEN
         CALL H5SCREATE_SIMPLE_F(RANK, CDIM, DATASPACE, ERROR, MDIM)
         CALL VTK_CHUNK_DIMS(RANK, DDIM, CHUNK)
      ELSE
         CALL H5SCREATE_SIMPLE_F(RANK, CDIM, DATASPACE, ERROR)
         CHUNK(1:RANK) = DDIM(1:RANK)
      ENDIF
      CALL GET_VTK_DCPL(CRP_LIST, RANK, CHUNK, TINY_CHUNK)
      CALL H5DCREATE_F(GROUP_ID, SNAME, H5T_STD_I32LE, DATASPACE, DSET_ID, ERROR, CRP_LIST)
      CALL H5PCLOSE_F(CRP_LIST, ERROR)
      CALL GET_VTK_DXPL(PLIST_ID)
      CALL H5SCLOSE_F(DATASPACE, ERROR)
      CALL VTK_DSET_CACHE_ADD(GROUP_ID, SNAME, DSET_ID)
   ELSE
      CALL H5DOPEN_F(GROUP_ID, SNAME, DSET_ID, ERROR)
      CALL GET_VTK_DXPL(PLIST_ID)
      CALL VTK_DSET_CACHE_ADD(GROUP_ID, SNAME, DSET_ID)
   ENDIF
END SUBROUTINE PARALLEL_INIT_I32

!> \brief Open, or create on first use, an unsigned 8 bit integer dataset

!>

!> \param GROUP_ID Group the dataset lives in

!> \param SNAME Dataset name

!> \param CRP_LIST Dataset creation property list, used when the dataset is created

!> \param RANK Dataset rank

!> \param DDIM Requested chunk dimensions

!> \param CDIM Initial dataset dimensions

!> \param DSET_ID Open dataset handle (out)

!> \param PLIST_ID Data transfer property list for the write that follows (out)

!> \param MDIM Maximum dimensions, present when the dataset grows over time


SUBROUTINE PARALLEL_INIT_U8(GROUP_ID, SNAME, CRP_LIST, RANK, DDIM, CDIM, DSET_ID, PLIST_ID, MDIM)
   INTEGER(HID_T), INTENT(IN) :: GROUP_ID               ! Memory identifiers
   INTEGER(HID_T), INTENT(INOUT) :: CRP_LIST, PLIST_ID  ! Memory identifiers
   INTEGER, INTENT(IN)     ::   RANK    ! Dataset rank
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: CDIM, DDIM
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN), OPTIONAL :: MDIM
   CHARACTER(LEN=*), INTENT(IN) :: SNAME
   INTEGER(HID_T), INTENT(OUT) :: DSET_ID  ! Memory identifiers
   INTEGER(HID_T) :: DATASPACE
   INTEGER(HSIZE_T), DIMENSION(RANK) :: CHUNK
   INTEGER     ::   ERROR ! Error flag
   LOGICAL :: LINK_EXISTS, TINY_CHUNK
   IF (VTK_DSET_CACHE_FIND(GROUP_ID, SNAME, DSET_ID)) THEN
      CALL GET_VTK_DXPL(PLIST_ID)
      RETURN
   ENDIF
   CALL H5LEXISTS_F(GROUP_ID, SNAME, LINK_EXISTS, ERROR)
   IF (.NOT.LINK_EXISTS) THEN
      TINY_CHUNK = VTK_IS_STEPS_GROUP(GROUP_ID)
      IF (PRESENT(MDIM)) THEN
         CALL H5SCREATE_SIMPLE_F(RANK, CDIM, DATASPACE, ERROR, MDIM)
         CALL VTK_CHUNK_DIMS(RANK, DDIM, CHUNK)
      ELSE
         CALL H5SCREATE_SIMPLE_F(RANK, CDIM, DATASPACE, ERROR)
         CHUNK(1:RANK) = DDIM(1:RANK)
      ENDIF
      CALL GET_VTK_DCPL(CRP_LIST, RANK, CHUNK, TINY_CHUNK)
      CALL H5DCREATE_F(GROUP_ID, SNAME, H5T_STD_U8LE, DATASPACE, DSET_ID, ERROR, CRP_LIST)
      CALL H5PCLOSE_F(CRP_LIST, ERROR)
      CALL GET_VTK_DXPL(PLIST_ID)
      CALL H5SCLOSE_F(DATASPACE, ERROR)
      CALL VTK_DSET_CACHE_ADD(GROUP_ID, SNAME, DSET_ID)
   ELSE 
      CALL H5DOPEN_F(GROUP_ID, SNAME, DSET_ID, ERROR)
      CALL GET_VTK_DXPL(PLIST_ID)
      CALL VTK_DSET_CACHE_ADD(GROUP_ID, SNAME, DSET_ID)
   ENDIF
END SUBROUTINE PARALLEL_INIT_U8

!> \brief Append integers to the end of a rank one dataset

!>

!> \param DSET_ID Open dataset handle

!> \param PLIST_ID Data transfer property list

!> \param INDATA Values to append

!> \param N Number of values


SUBROUTINE APPEND_RANK1_DATASET_I32(DSET_ID,PLIST_ID,INDATA,N)
   INTEGER, INTENT(IN)     ::   N    ! Dataset rank
   INTEGER(IB32), DIMENSION(N), INTENT(IN) :: INDATA
   INTEGER(HID_T), INTENT(INOUT) :: PLIST_ID  ! Memory identifiers
   INTEGER(HSIZE_T), DIMENSION(1) :: SIZE1, MDIM, CDIM, OFFSET, COUNT
   INTEGER(HID_T), INTENT(IN) :: DSET_ID  ! Memory identifiers
   INTEGER(HID_T) :: DATASPACE
   INTEGER     ::   ERROR ! Error flag
   
   COUNT = SIZE(INDATA)
   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE, CDIM, MDIM, ERROR)
   SIZE1 = CDIM(1)+COUNT
   OFFSET = CDIM
   CALL H5SCLOSE_F(DATASPACE, ERROR)
   CALL H5DSET_EXTENT_F(DSET_ID, SIZE1, ERROR)
   CALL PARALLEL_WRITE_I32(1, DSET_ID, PLIST_ID, OFFSET, COUNT, INDATA)
END SUBROUTINE APPEND_RANK1_DATASET_I32

!> \brief Append reals to the end of a rank one dataset

!>

!> \param DSET_ID Open dataset handle

!> \param PLIST_ID Data transfer property list

!> \param INDATA Values to append

!> \param N Number of values


SUBROUTINE APPEND_RANK1_DATASET_F32(DSET_ID,PLIST_ID,INDATA,N)
   INTEGER, INTENT(IN)     ::   N    ! Dataset rank
   REAL(FB), DIMENSION(N), INTENT(IN) :: INDATA
   INTEGER(HID_T), INTENT(INOUT) :: PLIST_ID  ! Memory identifiers
   INTEGER(HSIZE_T), DIMENSION(1) :: SIZE1, MDIM, CDIM, OFFSET, COUNT
   INTEGER(HID_T), INTENT(IN) :: DSET_ID  ! Memory identifiers
   INTEGER(HID_T) :: DATASPACE
   INTEGER     ::   ERROR ! Error flag
   
   COUNT = SIZE(INDATA)
   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE, CDIM, MDIM, ERROR)
   SIZE1 = CDIM(1)+COUNT
   OFFSET = CDIM
   CALL H5SCLOSE_F(DATASPACE, ERROR)
   CALL H5DSET_EXTENT_F(DSET_ID, SIZE1, ERROR)
   CALL PARALLEL_WRITE_F32(1, DSET_ID, PLIST_ID, OFFSET, COUNT, INDATA)
END SUBROUTINE APPEND_RANK1_DATASET_F32

!> \brief Write reals into a hyperslab of a dataset

!>

!> \param RANK Dataset rank

!> \param DSET_ID Open dataset handle

!> \param PLIST_ID Data transfer property list

!> \param OFFSET Where the hyperslab starts

!> \param COUNT Extent of the hyperslab

!> \param DATA Values to write


SUBROUTINE PARALLEL_WRITE_F32(RANK, DSET_ID,&
   PLIST_ID, OFFSET, COUNT, DATA)
   INTEGER(HID_T), INTENT(IN) :: DSET_ID, PLIST_ID  ! Memory identifiers
   INTEGER, INTENT(IN)     ::   RANK    ! Dataset rank
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: OFFSET, COUNT ! Attribute dimension
   REAL(FB), DIMENSION(*), INTENT(IN) :: DATA
   INTEGER     ::   ERROR ! Error flag
   INTEGER(HID_T) :: MEMSPACE, DATASPACE
   
   CALL H5SCREATE_SIMPLE_F (RANK, COUNT, MEMSPACE, ERROR)
   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SSELECT_HYPERSLAB_F(DATASPACE, H5S_SELECT_SET_F, OFFSET, COUNT, ERROR)
   CALL H5DWRITE_F(DSET_ID, H5T_IEEE_F32LE, DATA, COUNT, ERROR,&
      MEM_SPACE_ID = MEMSPACE, FILE_SPACE_ID = DATASPACE, XFER_PRP = PLIST_ID)
   CALL H5SCLOSE_F(MEMSPACE, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
END SUBROUTINE PARALLEL_WRITE_F32

!> \brief Write 32 bit integers into a hyperslab of a dataset

!>

!> \param RANK Dataset rank

!> \param DSET_ID Open dataset handle

!> \param PLIST_ID Data transfer property list

!> \param OFFSET Where the hyperslab starts

!> \param COUNT Extent of the hyperslab

!> \param DATA Values to write


SUBROUTINE PARALLEL_WRITE_I32(RANK, DSET_ID,&
   PLIST_ID, OFFSET, COUNT, DATA)
   INTEGER(HID_T), INTENT(IN) :: DSET_ID, PLIST_ID  ! Memory identifiers
   INTEGER, INTENT(IN)     ::   RANK    ! Dataset rank
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: OFFSET, COUNT ! Attribute dimension
   INTEGER(IB32), DIMENSION(*), INTENT(IN) :: DATA
   INTEGER     ::   ERROR ! Error flag
   INTEGER(HID_T) :: MEMSPACE, DATASPACE
   
   CALL H5SCREATE_SIMPLE_F (RANK, COUNT, MEMSPACE, ERROR)
   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SSELECT_HYPERSLAB_F(DATASPACE, H5S_SELECT_SET_F, OFFSET, COUNT, ERROR)
   CALL H5DWRITE_F(DSET_ID, H5T_STD_I32LE, DATA, COUNT, ERROR, &
        MEM_SPACE_ID = MEMSPACE, FILE_SPACE_ID = DATASPACE, XFER_PRP = PLIST_ID)
   CALL H5SCLOSE_F(MEMSPACE, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
END SUBROUTINE PARALLEL_WRITE_I32

!> \brief Write unsigned 8 bit integers into a hyperslab of a dataset

!>

!> \param RANK Dataset rank

!> \param DSET_ID Open dataset handle

!> \param PLIST_ID Data transfer property list

!> \param OFFSET Where the hyperslab starts

!> \param COUNT Extent of the hyperslab

!> \param DATA Values to write


SUBROUTINE PARALLEL_WRITE_U8(RANK, DSET_ID,&
   PLIST_ID, OFFSET, COUNT, DATA)
   INTEGER(HID_T), INTENT(IN) :: DSET_ID, PLIST_ID  ! Memory identifiers
   INTEGER, INTENT(IN)     ::   RANK    ! Dataset rank
   INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: OFFSET, COUNT ! Attribute dimension
   INTEGER(IB8), DIMENSION(*), INTENT(IN) :: DATA
   INTEGER     ::   ERROR ! Error flag
   INTEGER(HID_T) :: MEMSPACE, DATASPACE
   
   CALL H5SCREATE_SIMPLE_F (RANK, COUNT, MEMSPACE, ERROR)
   CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
   CALL H5SSELECT_HYPERSLAB_F(DATASPACE, H5S_SELECT_SET_F, OFFSET, COUNT, ERROR)
   CALL H5DWRITE_F(DSET_ID, H5T_STD_U8LE, DATA, COUNT, ERROR, &
        MEM_SPACE_ID = MEMSPACE, FILE_SPACE_ID = DATASPACE, XFER_PRP = PLIST_ID)
   CALL H5SCLOSE_F(MEMSPACE, ERROR)
   CALL H5SCLOSE_F(DATASPACE, ERROR)
END SUBROUTINE PARALLEL_WRITE_U8

!> \brief Write the VTKHDF format version as an attribute

!>

!> \param GROUP_ID Group to attach the attribute to

!> \param ADIMS Attribute dimensions

!> \param ARANK Attribute rank

!> \param ANAME Attribute name

!> \param ATTR_DATA Version numbers


SUBROUTINE ADD_VERSION(GROUP_ID,ADIMS,ARANK,ANAME,ATTR_DATA)
   INTEGER(HID_T) :: ATTR_ID, ASPACE_ID, GROUP_ID             ! Identifiers
   INTEGER(HSIZE_T), DIMENSION(1), INTENT(IN) :: ADIMS        ! Attribute dimension
   INTEGER, INTENT(IN)     ::   ARANK                         ! Attribute rank
   INTEGER :: ERROR                                           ! Error flag
   CHARACTER(LEN=*), INTENT(IN) :: ANAME                      ! Attribute name
   INTEGER(HSIZE_T), DIMENSION(:), INTENT(IN) :: ATTR_DATA(:) ! Attribute data
   
   CALL H5SCREATE_SIMPLE_F(ARANK, ADIMS, ASPACE_ID, ERROR)
   CALL H5ACREATE_F(GROUP_ID, TRIM(ANAME), H5T_STD_I32LE, ASPACE_ID, ATTR_ID, ERROR)
   CALL H5AWRITE_F(ATTR_ID, H5T_STD_I32LE, ATTR_DATA, ADIMS, ERROR)
   CALL H5ACLOSE_F(ATTR_ID, ERROR)
   CALL H5SCLOSE_F(ASPACE_ID, ERROR)
END SUBROUTINE ADD_VERSION

!> \brief Write a 32 bit real attribute

!>

!> \param GROUP_ID Group to attach the attribute to

!> \param ADIMS Attribute dimensions

!> \param ARANK Attribute rank

!> \param ANAME Attribute name

!> \param ATTR_DATA Attribute values


SUBROUTINE ADD_ATTRIBUTE_F32(GROUP_ID,ADIMS,ARANK,ANAME,ATTR_DATA)
   INTEGER(HID_T) :: ATTR_ID, ASPACE_ID, GROUP_ID             ! Identifiers
   INTEGER(HSIZE_T), DIMENSION(1), INTENT(IN) :: ADIMS        ! Attribute dimension
   INTEGER, INTENT(IN)     ::   ARANK                         ! Attribute rank
   INTEGER :: ERROR                                           ! Error flag
   CHARACTER(LEN=*), INTENT(IN) :: ANAME                      ! Attribute name
   REAL(FB), DIMENSION(:), INTENT(IN) :: ATTR_DATA            ! Attribute data

   CALL H5SCREATE_SIMPLE_F(ARANK, ADIMS, ASPACE_ID, ERROR)
   CALL H5ACREATE_F(GROUP_ID, TRIM(ANAME), H5T_IEEE_F32LE, ASPACE_ID, ATTR_ID, ERROR)
   CALL H5AWRITE_F(ATTR_ID, H5T_IEEE_F32LE, ATTR_DATA, ADIMS, ERROR)
   CALL H5ACLOSE_F(ATTR_ID, ERROR)
   CALL H5SCLOSE_F(ASPACE_ID, ERROR)
END SUBROUTINE ADD_ATTRIBUTE_F32

!> \brief Write a character attribute

!>

!> \param GROUP_ID Group to attach the attribute to

!> \param DATA_DIMS Attribute dimensions

!> \param ANAME Attribute name

!> \param ATTR_DATA Attribute value

!> \param ALEN Length of the attribute string


SUBROUTINE ADD_ATTRIBUTE_CHAR(GROUP_ID,DATA_DIMS,ANAME,ATTR_DATA,ALEN)
   INTEGER(HID_T) :: ATTR_ID, ASPACE_ID, GROUP_ID, ATYPE_ID   ! Identifiers
   INTEGER(HSIZE_T), DIMENSION(1), INTENT(IN) :: DATA_DIMS    ! Attribute dimension
   INTEGER(SIZE_T), INTENT(IN) :: ALEN     ! Length of the attribute string
   INTEGER     ::   ERROR                                     ! Error flag
   CHARACTER(LEN=*), INTENT(IN) :: ANAME
   CHARACTER(LEN=*), INTENT(IN) :: ATTR_DATA
   
   CALL H5SCREATE_F(H5S_SCALAR_F, ASPACE_ID, ERROR)
   CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, ATYPE_ID, ERROR)
   CALL H5TSET_SIZE_F(ATYPE_ID, ALEN, ERROR)
   CALL H5TSET_STRPAD_F(ATYPE_ID, H5T_STR_NULLPAD_F, ERROR)
   CALL H5ACREATE_F(GROUP_ID, TRIM(ANAME), ATYPE_ID, ASPACE_ID, ATTR_ID, ERROR)
   CALL H5AWRITE_F(ATTR_ID, ATYPE_ID, ATTR_DATA, DATA_DIMS, ERROR)
   CALL H5ACLOSE_F(ATTR_ID, ERROR)
   CALL H5TCLOSE_F(ATYPE_ID, ERROR)
   CALL H5SCLOSE_F(ASPACE_ID, ERROR)
END SUBROUTINE ADD_ATTRIBUTE_CHAR

!> \brief Write an integer attribute, replacing it if it is already there

!>

!> \param GROUP_ID Group to attach the attribute to

!> \param DATA_DIMS Attribute dimensions

!> \param ANAME Attribute name

!> \param ATTR_DATA Attribute value

!>

!> The number of output times written so far is kept this way, so it is rewritten

!> at every output time.


SUBROUTINE ADD_ATTRIBUTE_INT(GROUP_ID,DATA_DIMS,ANAME,ATTR_DATA)
   INTEGER(HID_T) :: ATTR_ID, ASPACE_ID, GROUP_ID   ! Identifiers
   INTEGER(HSIZE_T), DIMENSION(1), INTENT(IN) :: DATA_DIMS    ! Attribute dimension
   INTEGER     ::   ERROR                                     ! Error flag
   CHARACTER(LEN=*), INTENT(IN) :: ANAME
   INTEGER(HID_T), INTENT(IN) :: ATTR_DATA
   LOGICAL :: LINK_EXISTS
   
   CALL H5SCREATE_F(H5S_SCALAR_F, ASPACE_ID, ERROR)
   CALL H5AEXISTS_F(GROUP_ID, TRIM(ANAME), LINK_EXISTS, ERROR)
   IF (.NOT.LINK_EXISTS) THEN
      CALL H5ACREATE_F(GROUP_ID, TRIM(ANAME), H5T_STD_I64LE, ASPACE_ID, ATTR_ID, ERROR)
      CALL H5AWRITE_F(ATTR_ID, H5T_STD_I64LE, ATTR_DATA, DATA_DIMS, ERROR)
   ELSE
      CALL H5AOPEN_NAME_F(GROUP_ID, TRIM(ANAME), ATTR_ID, ERROR)
      CALL H5AWRITE_F(ATTR_ID, H5T_STD_I64LE, ATTR_DATA, DATA_DIMS, ERROR)
   ENDIF
   CALL H5ACLOSE_F(ATTR_ID, ERROR)
   CALL H5SCLOSE_F(ASPACE_ID, ERROR)
END SUBROUTINE ADD_ATTRIBUTE_INT

#endif


!> \brief The colour every particle of a class is drawn in
!>
!> \param N Particle class index
!>
!> Taken from the class, or from its surface when the class does not set one.  It does
!> not vary between particles of the class or over time, which is why the VTKHDF files
!> record it once rather than per point per output time.

FUNCTION GET_PART_CLASS_COLOR(N) RESULT(RGB)
INTEGER, INTENT(IN) :: N
REAL(FB), DIMENSION(3) :: RGB
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC_LOCAL
LPC_LOCAL => LAGRANGIAN_PARTICLE_CLASS(N)
IF (LPC_LOCAL%RGB(1)==-1) THEN
   RGB = REAL(SURFACE(LPC_LOCAL%SURF_INDEX)%RGB,FB)/255._FB
ELSE
   RGB = REAL(LPC_LOCAL%RGB,FB)/255._FB
ENDIF
END FUNCTION GET_PART_CLASS_COLOR


!> \brief Write a Python script that loads the VTKHDF output into ParaView


!>


!> \param NMESHES Number of meshes


!>


!> The script finds the files next to itself, so it works whether ParaView is


!> running locally or against a remote server.



SUBROUTINE WRITE_PARAVIEW_STATE_FILE(NMESHES)
USE OUTPUT_CLOCKS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

INTEGER, INTENT(IN) :: NMESHES
TYPE (MESH_TYPE), POINTER :: M
REAL(EB) :: CX,CY,CZ,XMN,XMX,YMN,YMX,ZMN,ZMX,CAM_R,CAM_D
INTEGER :: NM,N
REAL(FB) :: RGB_PART(3)
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

XMX = -HUGE(EB)
YMX = -HUGE(EB)
ZMX = -HUGE(EB)
XMN = HUGE(EB)
YMN = HUGE(EB)
ZMN = HUGE(EB)
DO NM=1,NMESHES
   CALL POINT_TO_MESH(NM)
   M => MESHES(NM)
   XMX = MAX(XMX, M%XF)
   XMN = MIN(XMN, M%XS)
   YMX = MAX(YMX, M%YF)
   YMN = MIN(YMN, M%YS)
   ZMX = MAX(ZMX, M%ZF)
   ZMN = MIN(ZMN, M%ZS)
ENDDO

CX = (XMX+XMN)/2
CY = (YMX+YMN)/2
CZ = (ZMX+ZMN)/2

! Stand the camera off along (-1,-1,1) far enough for the whole domain to be in
! frame.  _DisableFirstRenderCameraReset() below stops ParaView from framing the
! scene itself, so without an explicit position the camera is left at its default
! inside the domain and the view opens empty.
CAM_R = 0.5_EB*SQRT((XMX-XMN)**2+(YMX-YMN)**2+(ZMX-ZMN)**2)
IF (CAM_R<=0._EB) CAM_R = 1._EB
CAM_D = 3.9_EB*CAM_R/SQRT(3._EB)

OPEN(LU_PARAVIEW,FILE=FN_PARAVIEW,FORM='FORMATTED', STATUS='REPLACE',ACTION='WRITE')
WRITE(LU_PARAVIEW,'(A)') '#Script to import FDS generated data for visualization in Paraview'
WRITE(LU_PARAVIEW,'(A)') 'import os'
WRITE(LU_PARAVIEW,'(A)') 'import glob'
WRITE(LU_PARAVIEW,'(A)') "def writeSeries(files, times, outfile):"
WRITE(LU_PARAVIEW,'(A)') "    with open(outfile, 'w') as f:"
WRITE(LU_PARAVIEW,'(A,A,A,A,A)') '        f.write(',"'",&
                             '{\n  "file-series-version" : "1.0",\n  "files" : [',"'",')'
WRITE(LU_PARAVIEW,'(A)') "        for time, file in zip(times, files):"
WRITE(LU_PARAVIEW,'(A,A,A,A,A)') '            f.write(',"'",'    { "name" : "%s", "time" : %0.2f },\n',&
                              "'",'%(file, time))'
WRITE(LU_PARAVIEW,'(A)') "        f.write('    ]\n  }\n')"
WRITE(LU_PARAVIEW,'(A)') "def parseTimes(files, ext):"
WRITE(LU_PARAVIEW,'(A)') "    times = [float(file.split('_')[-1].split(ext)[0])/100 for file in files]"
WRITE(LU_PARAVIEW,'(A)') "    return times"

WRITE(LU_PARAVIEW,'(A,A,A)') "chid = '",TRIM(CHID),"'"
WRITE(LU_PARAVIEW,'(A)') 'T_Begin = 0.0'
IF (DT_VTK_SPECIFIED > 0) THEN
   WRITE(LU_PARAVIEW,'(A,F15.3)') 'T_End = ',(T_END-T_BEGIN)/DT_VTK_SPECIFIED
ELSE
   WRITE(LU_PARAVIEW,'(A,F15.3)') 'T_End = ',REAL(NFRAMES,FB)
ENDIF

WRITE(LU_PARAVIEW,'(A,F15.3,A,F15.3,A,F15.3,A)') 'CenterOfRotation = [',CX,',',CY,',',CZ,']'
WRITE(LU_PARAVIEW,'(A,F15.3,A,F15.3,A,F15.3,A)') 'CameraFocalPoint = [',CX,',',CY,',',CZ,']'
WRITE(LU_PARAVIEW,'(A,F15.3,A,F15.3,A,F15.3,A)') 'CameraPosition = [',CX-CAM_D,',',CY-CAM_D,',',CZ+CAM_D,']'
WRITE(LU_PARAVIEW,'(A)') 'diff = [abs(x-y) for x,y in zip(CenterOfRotation,CameraFocalPoint)]'
WRITE(LU_PARAVIEW,'(A)') 'if max(diff) < 0.1:'
WRITE(LU_PARAVIEW,'(A)') '    CameraFocalPoint[0] = CameraFocalPoint[0] + 0.1'
WRITE(LU_PARAVIEW,'(A)') '    CameraFocalPoint[1] = CameraFocalPoint[1] + 0.1'
WRITE(LU_PARAVIEW,'(A)') '    CameraFocalPoint[2] = CameraFocalPoint[2] + 0.1'
WRITE(LU_PARAVIEW,'(A)') 'CameraFocalPoint = [x+0.01 if (abs(x) < 0.01) else x for x in CameraFocalPoint]'
WRITE(LU_PARAVIEW,'(A)') 'import paraview'
WRITE(LU_PARAVIEW,'(A)') 'from paraview.simple import *'
WRITE(LU_PARAVIEW,'(A)') 'paraview.simple._DisableFirstRenderCameraReset()'
WRITE(LU_PARAVIEW,'(A)') 'version = paraview.simple.GetParaViewVersion()'
WRITE(LU_PARAVIEW,'(A)') 'version_num = version.major + version.minor/100'
WRITE(LU_PARAVIEW,'(A)') 'if version_num <= 5.11:'
WRITE(LU_PARAVIEW,'(A)') "    piecewisefunction = 'PiecewiseFunction'"
WRITE(LU_PARAVIEW,'(A)') "    axesactor = 'GridAxes3DActor'"
WRITE(LU_PARAVIEW,'(A)') "    gridaxesrep = 'GridAxesRepresentation'"
WRITE(LU_PARAVIEW,'(A)') "    polaraxesrep = 'PolarAxesRepresentation'"
WRITE(LU_PARAVIEW,'(A)') 'else:'
WRITE(LU_PARAVIEW,'(A)') "    piecewisefunction = 'Piecewise Function'"
WRITE(LU_PARAVIEW,'(A)') "    axesactor = 'Grid Axes 3D Actor'"
WRITE(LU_PARAVIEW,'(A)') "    gridaxesrep = 'Grid Axes Representation'"
WRITE(LU_PARAVIEW,'(A)') "    polaraxesrep = 'Polar Axes Representation'"
WRITE(LU_PARAVIEW,'(A)') 'materialLibrary1 = GetMaterialLibrary()'

WRITE(LU_PARAVIEW,'(A)') "renderView1 = CreateView('RenderView')"
WRITE(LU_PARAVIEW,'(A)') "renderView1.AxesGrid = axesactor"
WRITE(LU_PARAVIEW,'(A)') "renderView1.CenterOfRotation = CenterOfRotation"
WRITE(LU_PARAVIEW,'(A)') "renderView1.StereoType = 'Crystal Eyes'"
WRITE(LU_PARAVIEW,'(A)') "renderView1.CameraFocalPoint = CameraFocalPoint"
WRITE(LU_PARAVIEW,'(A)') "renderView1.CameraPosition = CameraPosition"
WRITE(LU_PARAVIEW,'(A)') "renderView1.CameraViewUp = [0.0, 0.0, 1.0]"
WRITE(LU_PARAVIEW,'(A)') "paraview.simple.LoadPalette('WhiteBackground')"
WRITE(LU_PARAVIEW,'(A)') "renderView1.BackEnd = 'OSPRay raycaster'"
WRITE(LU_PARAVIEW,'(A)') "renderView1.OSPRayMaterialLibrary = materialLibrary1"
WRITE(LU_PARAVIEW,'(A)') "SetActiveView(None)"
WRITE(LU_PARAVIEW,'(A)') "layout1 = CreateLayout(name='Layout #1')"
WRITE(LU_PARAVIEW,'(A)') "layout1.AssignView(0, renderView1)"
WRITE(LU_PARAVIEW,'(A)') "SetActiveView(renderView1)"
WRITE(LU_PARAVIEW,'(A)') "# ----------------------------------------------------------------"
WRITE(LU_PARAVIEW,'(A)') "# setup the data processing pipelines"
WRITE(LU_PARAVIEW,'(A)') "# ----------------------------------------------------------------"
WRITE(LU_PARAVIEW,'(A)') "remoteConnection = servermanager.ActiveConnection.IsRemote()"
WRITE(LU_PARAVIEW,'(A)') "if remoteConnection:"
WRITE(LU_PARAVIEW,'(A)') "    indir = r'"//TRIM(WORKING_DIR)//"'"
WRITE(LU_PARAVIEW,'(A)') "    sep = '/'"
WRITE(LU_PARAVIEW,'(A)') "    uri = servermanager.ActiveConnection.GetURI()"
WRITE(LU_PARAVIEW,'(A)') "else:"
WRITE(LU_PARAVIEW,'(A)') "    indir = os.path.dirname(os.path.realpath(__file__))"
WRITE(LU_PARAVIEW,'(A)') "    sep = os.sep"
WRITE(LU_PARAVIEW,'(A)') "    uri = None"
WRITE(LU_PARAVIEW,'(A)') "    # this script is written into VTK_DIR, so step back up to the case"
WRITE(LU_PARAVIEW,'(A)') "    # root and give indir the meaning it has for a remote connection"
WRITE(LU_PARAVIEW,'(A)') "    _rdir = r'"//TRIM(VTK_DIR)//"'.replace('/',sep).replace(chr(92),sep).strip(sep)"
WRITE(LU_PARAVIEW,'(A)') "    if _rdir and not os.path.isabs(_rdir) and indir.endswith(sep+_rdir):"
WRITE(LU_PARAVIEW,'(A)') "        indir = indir[:-(len(_rdir)+1)]"
WRITE(LU_PARAVIEW,'(A)') "rdir = r'"//TRIM(VTK_DIR)//"'"
WRITE(LU_PARAVIEW,'(A)') "if rdir == '':"
WRITE(LU_PARAVIEW,'(A)') "    namespace=indir+sep+chid"
WRITE(LU_PARAVIEW,'(A)') "else:"
WRITE(LU_PARAVIEW,'(A)') "    namespace=indir+sep+rdir+sep+chid"
WRITE(LU_PARAVIEW,'(A)') "pxm = servermanager.ProxyManager()"
WRITE(LU_PARAVIEW,'(A)') "directory_proxy = pxm.NewProxy('misc', 'ListDirectory')"
WRITE(LU_PARAVIEW,'(A)') "directory_proxy.List(indir+sep+rdir)"
WRITE(LU_PARAVIEW,'(A)') 'if version_num < 5.12:'
WRITE(LU_PARAVIEW,'(A)') "    directory_proxy.UpdatePropertyINFOrmation()"
WRITE(LU_PARAVIEW,'(A)') 'else:'
WRITE(LU_PARAVIEW,'(A)') "    directory_proxy.UpdatePropertyInformation()"
WRITE(LU_PARAVIEW,'(A,A)') "fileList = sorted(servermanager.VectorProperty(",&
                               "directory_proxy,directory_proxy.GetProperty('FileList')))"
WRITE(LU_PARAVIEW,'(A,A)') "directoryList = servermanager.VectorProperty(",&
                               "directory_proxy,directory_proxy.GetProperty('DirectoryList'))"
WRITE(LU_PARAVIEW,'(A)') "directory_proxy_root = pxm.NewProxy('misc', 'ListDirectory')"
WRITE(LU_PARAVIEW,'(A)') "directory_proxy_root.List(indir+sep)"
WRITE(LU_PARAVIEW,'(A)') 'if version_num < 5.12:'
WRITE(LU_PARAVIEW,'(A)') "    directory_proxy_root.UpdatePropertyINFOrmation()"
WRITE(LU_PARAVIEW,'(A)') 'else:'
WRITE(LU_PARAVIEW,'(A)') "    directory_proxy_root.UpdatePropertyInformation()"
WRITE(LU_PARAVIEW,'(A,A)') "fileList_root = sorted(servermanager.VectorProperty(",&
                               "directory_proxy_root,directory_proxy_root.GetProperty('FileList')))"
WRITE(LU_PARAVIEW,'(A,A)') "directoryList_root = servermanager.VectorProperty(",&
                               "directory_proxy_root,directory_proxy_root.GetProperty('DirectoryList'))"

WRITE(LU_PARAVIEW,'(A)') "# add geometry data"
WRITE(LU_PARAVIEW,'(A)') "if chid + '_GEOM.vtkhdf' in fileList_root:"
WRITE(LU_PARAVIEW,'(A,A)') "    geom = VTKHDFReader(registrationName='Geometry',",&
                      "FileName=[indir + sep + chid + '_GEOM.vtkhdf'])"
WRITE(LU_PARAVIEW,'(A)') "    geomDisplay = Show(geom, renderView1, 'UnstructuredGridRepresentation')"
WRITE(LU_PARAVIEW,'(A)') "    geomDisplay.MapScalars = 0"
WRITE(LU_PARAVIEW,'(A)') "    geomDisplay.Representation = 'Surface'"
WRITE(LU_PARAVIEW,'(A)') "    gcOLORTF2D = GetTransferFunction2D('Color')"
WRITE(LU_PARAVIEW,'(A)') "    geomColor = GetColorTransferFunction('Color')"
WRITE(LU_PARAVIEW,'(A)') "    geomColor.TransferFunction2D = gcOLORTF2D"
WRITE(LU_PARAVIEW,'(A)') "    geomColor.RGBPoints = [1.13, 0.23, 0.30, 0.75, 1.13, 0.87, 0.87, 0.87, 1.13, 0.71, 0.02, 0.15]"
WRITE(LU_PARAVIEW,'(A)') "    geomColor.ScalarRangeInitialized = 1.0"
WRITE(LU_PARAVIEW,'(A)') "    geomDisplay.ColorArrayName = ['CELLS', 'Color']"
WRITE(LU_PARAVIEW,'(A)') "    geomDisplay.LookupTable = geomColor"

WRITE(LU_PARAVIEW,'(A)') "# create a new 'STL Reader'"
WRITE(LU_PARAVIEW,'(A)') "if chid + '.stl' in fileList_root:"
WRITE(LU_PARAVIEW,'(A)') "    casestl = STLReader(registrationName='GeometrySTL', FileNames=[indir+sep+chid+'.stl'])"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay = Show(casestl, renderView1, 'GeometryRepresentation')"
WRITE(LU_PARAVIEW,'(A)') "    # trace defaults for the display properties."
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.Representation = 'Surface'"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.ColorArrayName = [None, '']"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.SelectTCoordArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.SelectNormalArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.SelectTangentArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.OSPRayScaleFunction = piecewisefunction"
!WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.Assembly = ''"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.SelectOrientationVectors = 'None'"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.ScaleFactor = 1.5493113040924074"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.SelectScaleArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.GlyphType = 'Arrow'"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.GlyphTableIndexArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.GaussianRadius = 0.07746556520462036"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.SetScaleArray = [None, '']"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.ScaleTransferFunction = piecewisefunction"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.OpacityArray = [None, '']"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.OpacityTransferFunction = piecewisefunction"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.DataAxesGrid = gridaxesrep"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.PolarAxes = polaraxesrep"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.SelectInputVectors = [None, '']"
WRITE(LU_PARAVIEW,'(A)') "    stlDisplay.WriteLog = ''"
WRITE(LU_PARAVIEW,'(A)') "# Load data files"
!WRITE(LU_PARAVIEW,'(A,A)') "sm3dFiles = [indir+sep+rdir+sep+x for x in fileList ",&
!                            "if ('_SM3D_' in x) and ('.vtkhdf' in x)]"
WRITE(LU_PARAVIEW,'(A)') "sm3dFiles = [indir+sep+rdir+sep+chid+'_SM3D.vtkhdf']"
WRITE(LU_PARAVIEW,'(A)') "if os.path.exists(sm3dFiles[0]) is False: sm3dFiles = []"
WRITE(LU_PARAVIEW,'(A,A)') "sl2dxFiles = [indir+sep+rdir+sep+x for x in fileList ",&
                            "if ('_X_' in x) and ('.vtkhdf' in x)]"
WRITE(LU_PARAVIEW,'(A,A)') "sl2dyFiles = [indir+sep+rdir+sep+x for x in fileList ",&
                            "if ('_Y_' in x) and ('.vtkhdf' in x)]"
WRITE(LU_PARAVIEW,'(A,A)') "sl2dzFiles = [indir+sep+rdir+sep+x for x in fileList ",&
                            "if ('_Z_' in x) and ('.vtkhdf' in x)]"
WRITE(LU_PARAVIEW,'(A,A)') "sl2daFiles = [indir+sep+rdir+sep+x for x in fileList ",&
                            "if ('_AGL_' in x) and ('.vtkhdf' in x)]"
WRITE(LU_PARAVIEW,'(A,A)') "sl3dFiles = [indir+sep+rdir+sep+x for x in fileList ",&
                            "if ('_SL3D_' in x) and ('.vtkhdf' in x)]"
WRITE(LU_PARAVIEW,'(A)') "bndfFiles = [indir+sep+rdir+sep+chid+'_BNDF.vtkhdf']"
WRITE(LU_PARAVIEW,'(A,A)') "partFiles = [indir+sep+rdir+sep+x for x in fileList ",&
                               "if ('_PART_' in x) and ('.vtkhdf' in x)]"
!WRITE(LU_PARAVIEW,'(A)') "bndfFiles = sorted(glob.glob(namespace+'_BNDF_*.pvtu'))"
!WRITE(LU_PARAVIEW,'(A)') "sm3dFiles = sorted(glob.glob(namespace+'_SM3D_*.pvtu'))"
!WRITE(LU_PARAVIEW,'(A)') "sl2dxFiles = sorted(glob.glob(namespace+'_X_*.pvtu'))"
!WRITE(LU_PARAVIEW,'(A)') "sl2dyFiles = sorted(glob.glob(namespace+'_Y_*.pvtu'))"
!WRITE(LU_PARAVIEW,'(A)') "sl2dzFiles = sorted(glob.glob(namespace+'_Z_*.pvtu'))"
!WRITE(LU_PARAVIEW,'(A)') "sl3dFiles = sorted(glob.glob(namespace+'_SL3D_*.pvtu'))"
!WRITE(LU_PARAVIEW,'(A)') "partFiles = sorted(glob.glob(namespace+'_PART_*.pvtp'))"

WRITE(LU_PARAVIEW,'(A)') "# Add boundary data"
WRITE(LU_PARAVIEW,'(A)') "if len(bndfFiles) > 0:"
WRITE(LU_PARAVIEW,'(A)') "    if remoteConnection:"
WRITE(LU_PARAVIEW,'(A,A)') "        BoundaryData = VTKHDFReader(",&
                                        "registrationName='Boundary', FileName=bndfFiles)"
WRITE(LU_PARAVIEW,'(A)') "    else:"
!WRITE(LU_PARAVIEW,'(A)') "        bndfFiles = [rdir + x.split(sep)[-1] for x in bndfFiles]"
!WRITE(LU_PARAVIEW,'(A)') "        times = parseTimes(bndfFiles, '.vtkhdf')"
!WRITE(LU_PARAVIEW,'(A)') "        outname = indir+sep+'bndf.vtkhdf.series'"
!WRITE(LU_PARAVIEW,'(A)') "        writeSeries(bndfFiles, times, outname)"
WRITE(LU_PARAVIEW,'(A,A)') "        BoundaryData = VTKHDFReader(",&
                                        "registrationName='Boundary', FileName=bndfFiles)"
WRITE(LU_PARAVIEW,'(A)') "# Add smoke 3d data"
WRITE(LU_PARAVIEW,'(A)') "if len(sm3dFiles) > 0:"
WRITE(LU_PARAVIEW,'(A)') "    if remoteConnection:"
WRITE(LU_PARAVIEW,'(A,A)') "        sm3dData = VTKHDFReader(",&
                                     "registrationName='Raw Smoke 3D', FileName=sm3dFiles)"
WRITE(LU_PARAVIEW,'(A)') "    else:"
!WRITE(LU_PARAVIEW,'(A)') "        sm3dFiles = [rdir + x.split(sep)[-1] for x in sm3dFiles]"
!WRITE(LU_PARAVIEW,'(A)') "        times = parseTimes(sm3dFiles, '.vtkhdf')"
!WRITE(LU_PARAVIEW,'(A)') "        outname = indir+sep+'sm3d.vtkhdf.series'"
!WRITE(LU_PARAVIEW,'(A)') "        writeSeries(sm3dFiles, times, outname)"
!WRITE(LU_PARAVIEW,'(A,A)') "        sm3dData = VTKHDFReader(",&
!                                     "registrationName='Raw Smoke 3D', FileName=[outname])"
WRITE(LU_PARAVIEW,'(A,A)') "        sm3dData = VTKHDFReader(",&
                                     "registrationName='Raw Smoke 3D', FileName=sm3dFiles)"
WRITE(LU_PARAVIEW,'(A)') "    smokeName = None"
WRITE(LU_PARAVIEW,'(A)') "    fireName = None"
WRITE(LU_PARAVIEW,'(A)') "    for s in sm3dData.PointArrayStatus:"
WRITE(LU_PARAVIEW,'(A)') "        if ('smoke' in s.lower() or 'soot' in s.lower()):"
WRITE(LU_PARAVIEW,'(A)') "            smokeName = s"
WRITE(LU_PARAVIEW,'(A)') "        if ('hrrpuv' in s.lower()):"
WRITE(LU_PARAVIEW,'(A)') "            fireName = s"
WRITE(LU_PARAVIEW,'(A)') "# Add 3d slice data"
WRITE(LU_PARAVIEW,'(A)') "if len(sl3dFiles) > 0:"
WRITE(LU_PARAVIEW,'(A)') "    if remoteConnection:"
WRITE(LU_PARAVIEW,'(A,A)') "        sl3dData = VTKHDFReader(",&
                                        "registrationName='Raw 3D Slice', FileName=sl3dFiles)"
WRITE(LU_PARAVIEW,'(A)') "    else:"
!WRITE(LU_PARAVIEW,'(A)') "        sl3dFiles = [rdir + x.split(sep)[-1] for x in sl3dFiles]"
!WRITE(LU_PARAVIEW,'(A)') "        times = parseTimes(sl3dFiles, '.vtkhdf')"
!WRITE(LU_PARAVIEW,'(A)') "        outname = indir+sep+'sl3d.vtkhdf.series'"
!WRITE(LU_PARAVIEW,'(A)') "        writeSeries(sl3dFiles, times, outname)"
WRITE(LU_PARAVIEW,'(A,A)') "        sl3dData = VTKHDFReader(",&
                                     "registrationName='Raw 3D Slice', FileName=sl3dFiles)"
WRITE(LU_PARAVIEW,'(A,A)') "    sl3dImage = ResampleToImage(",&
                                    "registrationName='Sampled 3D Slice', Input=sl3dData)"
WRITE(LU_PARAVIEW,'(A,A)') "    sl3dSlice = Slice(",&
                                    "registrationName='3D Slice Extraction', Input=sl3dImage)"
WRITE(LU_PARAVIEW,'(A)') "    sl3dSlice.SliceType = 'Plane'"
WRITE(LU_PARAVIEW,'(A)') "    sl3dSlice.HyperTreeGridSlicer = 'Plane'"
WRITE(LU_PARAVIEW,'(A)') "    sl3dSlice.SliceOffsetValues = [0.0]"
!WRITE(LU_PARAVIEW,'(A)') "    sl3dSlice.PointMergeMethod = 'Uniform Binning'"
WRITE(LU_PARAVIEW,'(A)') "    sl3dSlice.SliceType.Origin = CenterOfRotation"
WRITE(LU_PARAVIEW,'(A)') "    sl3dSlice.HyperTreeGridSlicer.Origin = CenterOfRotation"

WRITE(LU_PARAVIEW,'(A)') "# Add 2d slice data"
WRITE(LU_PARAVIEW,'(A,A)') "for sl2dFiles, axis_name in zip([sl2dxFiles,sl2dyFiles,sl2dzFiles,sl2daFiles],",&
                           "    ['X','Y','Z','AGL']):"
WRITE(LU_PARAVIEW,'(A)') "    if len(sl2dFiles) > 0:"
WRITE(LU_PARAVIEW,'(A,A)') "        slcfTypes = slcfTypes = [x.split(chid+'_'+axis_name+'_')[1]",&
                                  ".replace('.vtkhdf','') for x in sl2dFiles]"
WRITE(LU_PARAVIEW,'(A)') "        slcfTypes = [x.replace('neg_','-').replace('pos_','') for x in slcfTypes]"
WRITE(LU_PARAVIEW,'(A)') "        uniqueSlcfTypes = sorted(list(set(slcfTypes)))"
WRITE(LU_PARAVIEW,'(A)') "        for slcfType in uniqueSlcfTypes:"
WRITE(LU_PARAVIEW,'(A)') "            axis=float(slcfType)/100"
WRITE(LU_PARAVIEW,'(A)') "            if remoteConnection:"
WRITE(LU_PARAVIEW,'(A,A)') "                slcf_files = sorted([x for x,y in zip(sl2dFiles, slcfTypes)",&
                                                "if y == slcfType])"
WRITE(LU_PARAVIEW,'(A,A)') "                sl2dData = VTKHDFReader(",&
                                                "registrationName='%s=%0.4f'%(axis_name,axis), FileName=slcf_files)"
WRITE(LU_PARAVIEW,'(A)') "            else:"
WRITE(LU_PARAVIEW,'(A,A)') "                slcf_files = sorted([indir+sep+rdir + x.split(sep)[-1] for x,y in ",&
                                                "zip(sl2dFiles, slcfTypes) if y == slcfType])"
!WRITE(LU_PARAVIEW,'(A)') "                times = parseTimes(slcf_files, '.vtkhdf')"
!WRITE(LU_PARAVIEW,'(A)') "                outname = indir+sep+'sl2d-'+slcfType.replace(' ','-')+'.vtkhdf.series'"
!WRITE(LU_PARAVIEW,'(A)') "                writeSeries(slcf_files, times, outname)"
WRITE(LU_PARAVIEW,'(A,A)') "                sl2dData = VTKHDFReader(",&
                                                "registrationName='%s=%0.4f'%(axis_name,axis), FileName=slcf_files)"
WRITE(LU_PARAVIEW,'(A)') "# Add particle data"
! Every particle of a class is drawn in the class's colour, which the writer no longer
! repeats per point per output time, so name the colours here.  The file also carries
! its own in VTKHDF/FieldData/COLOR.

WRITE(LU_PARAVIEW,'(A)') "partColors = {}"
DO N=1,N_LAGRANGIAN_CLASSES
   RGB_PART = GET_PART_CLASS_COLOR(N)
   WRITE(LU_PARAVIEW,'(A,A,A,F6.3,A,F6.3,A,F6.3,A)') "partColors['",TRIM(LAGRANGIAN_PARTICLE_CLASS(N)%ID),&
      "'] = [",RGB_PART(1),", ",RGB_PART(2),", ",RGB_PART(3),"]"
ENDDO

WRITE(LU_PARAVIEW,'(A)') "if len(partFiles) > 0:"
WRITE(LU_PARAVIEW,'(A)') "    for partFile in partFiles:"
WRITE(LU_PARAVIEW,'(A)') "        partType = partFile.split('_PART_')[1].rsplit('.vtkhdf', 1)[0]"
WRITE(LU_PARAVIEW,'(A,A)') "        partData = VTKHDFReader(",&
                                        "registrationName='Particle: '+partType, FileName=[partFile])"
WRITE(LU_PARAVIEW,'(A)') "        partDisplay = Show(partData, renderView1, 'GeometryRepresentation')"
WRITE(LU_PARAVIEW,'(A)') "        partDisplay.Representation = 'Point Gaussian'"
WRITE(LU_PARAVIEW,'(A)') "        partDisplay.ColorArrayName = [None, '']"
WRITE(LU_PARAVIEW,'(A)') "        partRGB = partColors.get(partType, [0.7, 0.7, 0.7])"
WRITE(LU_PARAVIEW,'(A)') "        partDisplay.AmbientColor = partRGB"
WRITE(LU_PARAVIEW,'(A)') "        partDisplay.DiffuseColor = partRGB"
WRITE(LU_PARAVIEW,'(A)') "        partDisplay.GaussianRadius = 0.05"
WRITE(LU_PARAVIEW,'(A)') "        partDisplay.ShaderPreset = 'Plain circle'"

WRITE(LU_PARAVIEW,'(A)') "if len(sm3dFiles) > 0:"
WRITE(LU_PARAVIEW,'(A)') "    # ----------------------------------------------------------------"
WRITE(LU_PARAVIEW,'(A)') "    # setup filter for fire 3d data"
WRITE(LU_PARAVIEW,'(A)') "    # ----------------------------------------------------------------"
WRITE(LU_PARAVIEW,'(A)') "    if fireName is not None:"
WRITE(LU_PARAVIEW,'(A)') "        fireImage = ResampleToImage(registrationName='Fire', Input=sm3dData)"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay = Show(fireImage, renderView1, 'UniformGridRepresentation')"
WRITE(LU_PARAVIEW,'(A)') "        # get 2D transfer function for 'HRRPUV'"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVTF2D = GetTransferFunction2D(fireName, separate=True)"
WRITE(LU_PARAVIEW,'(A)') "        # get color transfer function/color map for 'HRRPUV'"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVLUT = GetColorTransferFunction(fireName, separate=True)"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVLUT.TransferFunction2D = hRRPUVTF2D"
WRITE(LU_PARAVIEW,'(A,A)') "        hRRPUVLUT.RGBPoints = [0.0, 0.0, 0.0, 0.0, 25.4, 0.9, 0.0, 0.0, 76.2, ",&
                                  "0.9, 0.9, 0.0, 254.0, 1.0, 1.0, 1.0]"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVLUT.ColorSpace = 'RGB'"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVLUT.NanColor = [0.0, 0.5, 1.0]"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVLUT.ScalarRangeInitialized = 1.0"
WRITE(LU_PARAVIEW,'(A)') "        # get opacity transfer function/opacity map for 'HRRPUV'"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVPWF = GetOpacityTransferFunction(fireName, separate=True)"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVPWF.Points = [0.0, 0.0, 0.5, 0.0, 254.0, 1.0, 0.5, 0.0]"
WRITE(LU_PARAVIEW,'(A)') "        hRRPUVPWF.ScalarRangeInitialized = 1"
WRITE(LU_PARAVIEW,'(A)') "        # trace defaults for the display properties."
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.Representation = 'Volume'"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.ColorArrayName = ['POINTS', fireName]"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.LookupTable = hRRPUVLUT"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SelectTCoordArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SelectNormalArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SelectTangentArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.OSPRayScaleArray = fireName"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.OSPRayScaleFunction = piecewisefunction"
!WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.Assembly = ''"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SelectOrientationVectors = 'None'"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.ScaleFactor = 3.0"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SelectScaleArray = fireName"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.GlyphType = 'Arrow'"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.GlyphTableIndexArray = fireName"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.GaussianRadius = 0.15"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SetScaleArray = ['POINTS', fireName]"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.ScaleTransferFunction = piecewisefunction"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.OpacityArray = ['POINTS', fireName]"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.OpacityTransferFunction = piecewisefunction"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.DataAxesGrid = gridaxesrep"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.PolarAxes = polaraxesrep"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.ScalarOpacityUnitDistance = 0.44"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.ScalarOpacityFunction = hRRPUVPWF"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.TransferFunction2D = hRRPUVTF2D"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.OpacityArrayName = ['POINTS', fireName]"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.ColorArray2Name = ['POINTS', fireName]"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SliceFunction = 'Plane'"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.Slice = 49"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SelectInputVectors = ['POINTS', '']"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.WriteLog = ''"
WRITE(LU_PARAVIEW,'(A)') "        # init the piecewisefunction selected for 'ScaleTransferFunction'"
WRITE(LU_PARAVIEW,'(A,A)') "        fireImageDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.0, ",&
                                  "1.0, 0.5, 0.0]"
WRITE(LU_PARAVIEW,'(A)') "        # init the piecewisefunction selected for 'OpacityTransferFunction'"
WRITE(LU_PARAVIEW,'(A,A)') "        fireImageDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.0, ",&
                                  "1.0, 0.5, 0.0]"
WRITE(LU_PARAVIEW,'(A)') "        # init the 'Plane' selected for 'SliceFunction'"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.SliceFunction.Origin = CenterOfRotation"
WRITE(LU_PARAVIEW,'(A)') "        fireImageDisplay.UseSeparateColorMap = True"
WRITE(LU_PARAVIEW,'(A)') "    # ----------------------------------------------------------------"
WRITE(LU_PARAVIEW,'(A)') "    # setup filter for smoke 3d data"
WRITE(LU_PARAVIEW,'(A)') "    # ----------------------------------------------------------------"
WRITE(LU_PARAVIEW,'(A)') "    if smokeName is not None:"
WRITE(LU_PARAVIEW,'(A)') "        smokeImage = ResampleToImage(registrationName='Smoke', Input=sm3dData)"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay = Show(smokeImage, renderView1, 'UniformGridRepresentation')"
WRITE(LU_PARAVIEW,'(A)') "        # get 2D transfer function for 'SOOTDENSITY'"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYTF2D = GetTransferFunction2D(smokeName)"
WRITE(LU_PARAVIEW,'(A)') "        # get color transfer function/color map for 'SOOTDENSITY'"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYLUT = GetColorTransferFunction(smokeName)"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYLUT.TransferFunction2D = sOOTDENSITYTF2D"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYLUT.RGBPoints = [0.0, 0.0, 0.0, 0.0, 254.0, 0.0, 0.0, 0.0]"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYLUT.ColorSpace = 'RGB'"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYLUT.NanColor = [1.0, 0.0, 0.0]"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYLUT.ScalarRangeInitialized = 1.0"
WRITE(LU_PARAVIEW,'(A)') "        # get opacity transfer function/opacity map for 'SOOTDENSITY'"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYPWF = GetOpacityTransferFunction(smokeName)"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 0.25, 0.5, 0.0, 254.0, 1.0, 0.5, 0.0]"
WRITE(LU_PARAVIEW,'(A)') "        sOOTDENSITYPWF.ScalarRangeInitialized = 1"
WRITE(LU_PARAVIEW,'(A)') "        # trace defaults for the display properties."
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.Representation = 'Volume'"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.ColorArrayName = ['POINTS', smokeName]"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.LookupTable = sOOTDENSITYLUT"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.SelectTCoordArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.SelectNormalArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.SelectTangentArray = 'None'"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.OSPRayScaleArray = smokeName"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.OSPRayScaleFunction = piecewisefunction"
!WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.Assembly = ''"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.SelectOrientationVectors = 'None'"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.ScaleFactor = 3.0"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.SelectScaleArray = smokeName"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.GlyphType = 'Arrow'"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.GlyphTableIndexArray = smokeName"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.GaussianRadius = 0.15"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.SetScaleArray = ['POINTS', smokeName]"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.ScaleTransferFunction = piecewisefunction"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.OpacityArray = ['POINTS', smokeName]"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.OpacityTransferFunction = piecewisefunction"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.DataAxesGrid = gridaxesrep"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.PolarAxes = polaraxesrep"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.ScalarOpacityFunction = sOOTDENSITYPWF"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.ScalarOpacityUnitDistance = 1.2957383373220388"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.OpacityArrayName = ['POINTS', smokeName]"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.SelectInputVectors = [None, '']"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.WriteLog = ''"
WRITE(LU_PARAVIEW,'(A)') "        # init the piecewisefunction selected for 'ScaleTransferFunction'"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.ScaleTransferFunction.Points =[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]"
WRITE(LU_PARAVIEW,'(A)') "        # init the piecewisefunction selected for 'OpacityTransferFunction'"
WRITE(LU_PARAVIEW,'(A)') "        smokeImageDisplay.OpacityTransferFunction.Points =[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]"

WRITE(LU_PARAVIEW,'(A)') "# ----------------------------------------------------------------"
WRITE(LU_PARAVIEW,'(A)') "# setup animation scene, tracks and keyframes"
WRITE(LU_PARAVIEW,'(A)') "# note: the Get..() functions create a new object, if needed"
WRITE(LU_PARAVIEW,'(A)') "# ----------------------------------------------------------------"

WRITE(LU_PARAVIEW,'(A)') "# get time animation track"
WRITE(LU_PARAVIEW,'(A)') "timeAnimationCue1 = GetTimeTrack()"

WRITE(LU_PARAVIEW,'(A)') "# initialize the animation scene"

WRITE(LU_PARAVIEW,'(A)') "# get the time-keeper"
WRITE(LU_PARAVIEW,'(A)') "timeKeeper1 = GetTimeKeeper()"
WRITE(LU_PARAVIEW,'(A)') "# initialize the timekeeper"
WRITE(LU_PARAVIEW,'(A)') "# initialize the animation track"
WRITE(LU_PARAVIEW,'(A)') "# get animation scene"
WRITE(LU_PARAVIEW,'(A)') "animationScene1 = GetAnimationScene()"

WRITE(LU_PARAVIEW,'(A)') "# initialize the animation scene"
WRITE(LU_PARAVIEW,'(A)') "animationScene1.ViewModules = renderView1"
WRITE(LU_PARAVIEW,'(A)') "animationScene1.Cues = timeAnimationCue1"
WRITE(LU_PARAVIEW,'(A)') "animationScene1.AnimationTime = T_Begin"
!WRITE(LU_PARAVIEW,'(A)') "animationScene1.EndTime = T_End"
WRITE(LU_PARAVIEW,'(A)') "animationScene1.PlayMode = 'Snap To TimeSteps'"

WRITE(LU_PARAVIEW,'(A)') "# ----------------------------------------------------------------"
WRITE(LU_PARAVIEW,'(A)') "# restore active source"
WRITE(LU_PARAVIEW,'(A)') "if len(sm3dFiles) > 0:"
WRITE(LU_PARAVIEW,'(A)') "    SetActiveSource(sm3dData)"
WRITE(LU_PARAVIEW,'(A)') "# ----------------------------------------------------------------"

CLOSE(LU_PARAVIEW)

T_USED(7) = T_USED(7) + CURRENT_TIME() - TNOW

END SUBROUTINE WRITE_PARAVIEW_STATE_FILE





#ifdef WITH_HDF5

!------------------------------------------------------------------------------
! VTKHDF output drivers.  These were previously carried in dump.f90; they live
! here so that all HDF5-specific output code is in one module.
!------------------------------------------------------------------------------

!> \brief Append one output time to the VTKHDF files

!>

!> \param T Current simulation time (s)

!> \param DT Current time step size (s)

!> \param INITIAL_DUMP True on the call that creates the files

!>

!> Called at every time step.  Each output has its own clock, so most calls write

!> nothing; the ones that do write are followed by a flush of just those files.


SUBROUTINE DUMP_VTK_MESH_OUTPUTS_SERIES(T,DT,INITIAL_DUMP)

USE COMP_FUNCTIONS, ONLY : CURRENT_TIME
REAL(EB) :: TNOW
REAL(EB), INTENT(IN) :: T,DT
LOGICAL, INTENT(IN) :: INITIAL_DUMP
LOGICAL :: WROTE_SLCF,WROTE_SM3D,WROTE_BNDF,WROTE_PART
INTEGER :: NM

TNOW = CURRENT_TIME()

! Which files this call writes to, so that only those are flushed at the end of it

WROTE_SLCF = .FALSE. ; WROTE_SM3D = .FALSE. ; WROTE_BNDF = .FALSE. ; WROTE_PART = .FALSE.

IF (INITIAL_DUMP) THEN
   CALL INITIALIZE_SMOKE3D_VTKHDF_SERIES()
   CALL INITIALIZE_SLCF_VTKHDF_SERIES()
   CALL INITIALIZE_BNDF_VTKHDF_SERIES()
   CALL INITIALIZE_VTKHDF_PART()
   WROTE_SLCF = .TRUE. ; WROTE_SM3D = .TRUE. ; WROTE_BNDF = .TRUE. ; WROTE_PART = .TRUE.
ENDIF

! VTK 3-D slices
IF (T>=SL3D_VTK_CLOCK(SL3D_VTK_COUNTER(LOWER_MESH_INDEX)) .OR. STOP_STATUS==INSTABILITY_STOP) THEN
   IF (.NOT.VTK_KEEPOPEN) CALL OPEN_VTKHDF_SLICE()
   CALL DUMP_SLCF_VTK(T,DT,0)
   WROTE_SLCF = .TRUE.
   IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_SLICE()
   SL3D_VTK_COUNTER(LOWER_MESH_INDEX) = SL3D_VTK_COUNTER(LOWER_MESH_INDEX) + 1
ENDIF

! VTK 2-D slices
IF (T>=SLCF_VTK_CLOCK(SLCF_VTK_COUNTER(LOWER_MESH_INDEX)) .OR. STOP_STATUS==INSTABILITY_STOP) THEN
   IF (.NOT.VTK_KEEPOPEN) CALL OPEN_VTKHDF_SLICE()
   CALL DUMP_SLCF_VTK(T,DT,1)
   WROTE_SLCF = .TRUE.
   IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_SLICE()
   SLCF_VTK_COUNTER(LOWER_MESH_INDEX) = SLCF_VTK_COUNTER(LOWER_MESH_INDEX) + 1
ENDIF

! VTK Smoke 3D slices
IF (T>=SM3D_VTK_CLOCK(SM3D_VTK_COUNTER(LOWER_MESH_INDEX)) .OR. STOP_STATUS==INSTABILITY_STOP) THEN
   IF (N_SMOKE3D > 0) THEN
      IF (.NOT.VTK_KEEPOPEN) CALL OPEN_VTKHDF_SMOKE3D()
      CALL DUMP_SMOKE3D_VTKHDF(T,DT)
      WROTE_SM3D = .TRUE.
      IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_SMOKE3D()
   ENDIF
   !CALL H5FFLUSH_F(HDF_SM3D_FILE_ID,H5F_SCOPE_GLOBAL_F,ERROR)
   SM3D_VTK_COUNTER(LOWER_MESH_INDEX) = SM3D_VTK_COUNTER(LOWER_MESH_INDEX) + 1
ENDIF

! VTK Boundary data
IF (T>=BNDF_VTK_CLOCK(BNDF_VTK_COUNTER(LOWER_MESH_INDEX))) THEN
   IF (.NOT.VTK_KEEPOPEN) CALL OPEN_VTKHDF_BNDF()
   CALL DUMP_BNDF_VTKHDF(T,DT)
   WROTE_BNDF = .TRUE.
   IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_BNDF()
   BNDF_VTK_COUNTER(LOWER_MESH_INDEX) = BNDF_VTK_COUNTER(LOWER_MESH_INDEX) + 1
ENDIF

! VTK Particle data
IF (T>=PART_VTK_CLOCK(PART_VTK_COUNTER(LOWER_MESH_INDEX))) THEN
   IF (.NOT.VTK_KEEPOPEN) CALL OPEN_VTKHDF_PART()
   CALL DUMP_PART_VTKHDF(T)
   WROTE_PART = .TRUE.
   IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_PART()
   PART_VTK_COUNTER(LOWER_MESH_INDEX) = PART_VTK_COUNTER(LOWER_MESH_INDEX) + 1
ENDIF

! Make what has been written so far readable on disk.  Without this a run that is killed
! leaves files that HDF5 cannot open at all, because with MPI-IO the superblock and the
! group and object headers are not written until the file is closed.

IF (VTK_KEEPOPEN .AND. FLUSH_FILE_BUFFERS) CALL FLUSH_VTKHDF(WROTE_SLCF,WROTE_SM3D,WROTE_BNDF,WROTE_PART)

! Spreadsheet data is written here only when the Smokeview path is switched off.
! When WRITE_SMV is true, DUMP_MESH_OUTPUTS writes it instead.

IF (.NOT.WRITE_SMV) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL DUMP_MESH_SPREADSHEET_OUTPUTS(T,NM)
   ENDDO
ENDIF
T_USED(7) = T_USED(7) + CURRENT_TIME() - TNOW
END SUBROUTINE DUMP_VTK_MESH_OUTPUTS_SERIES













!> \brief Dump Lagrangian particle data to CHID.prt5
!>
!> \param T Current simulation time (s)







SUBROUTINE INITIALIZE_SMOKE3D_VTKHDF_SERIES()

! Only generate file if at least one SMOKE3D quantity requested
IF (N_SMOKE3D > 0) THEN
   CALL INITIALIZE_VTKHDF_SMOKE3D()
ELSE
   RETURN
ENDIF

IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_SMOKE3D()

END SUBROUTINE INITIALIZE_SMOKE3D_VTKHDF_SERIES

!> \brief Set up the boundary VTKHDF output for the run


SUBROUTINE INITIALIZE_BNDF_VTKHDF_SERIES()

CHARACTER(FN_LENGTH) :: FILENAME

WRITE(FILENAME,'(A,A,A)') "",TRIM(VTK_DIR)//TRIM(CHID),'_BNDF.vtkhdf'
CALL INITIALIZE_VTKHDF_BNDF(FILENAME)
IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_BNDF()

END SUBROUTINE INITIALIZE_BNDF_VTKHDF_SERIES


!> \brief Set up the slice VTKHDF output for the run, one file per unique slice plane



SUBROUTINE INITIALIZE_SLCF_VTKHDF_SERIES()

CHARACTER(FN_LENGTH) :: FILENAME,SLCFNAME
INTEGER :: IQ,NTSL

ALLOCATE(HDF_SLCF_FILE_ID(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_PLIST_ID(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_CRP_LIST(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G1(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G2(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G3(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G4(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G5(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G6(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G7(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_COUNTER(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(HDF_SLCF_G1_NCELLS(MESHES(1)%N_UNIQUE_SLCF,NMESHES))
ALLOCATE(HDF_SLCF_G1_NPOINTS(MESHES(1)%N_UNIQUE_SLCF,NMESHES))
HDF_SLCF_G1_NCELLS = 0_IB32
HDF_SLCF_G1_NPOINTS = 0_IB32
NTSL = 0
UNIQUE_LOOPF: DO IQ=1,MESHES(1)%N_UNIQUE_SLCF
   SLCFNAME = MESHES(1)%UNIQUE_SLICE_NAMES(IQ)
   WRITE(FILENAME,'(A,A,A,A)') TRIM(VTK_DIR)//TRIM(CHID),'_',&
      TRIM(SLCFNAME),'.vtkhdf'
   IF (MESHES(1)%UNIQUE_SLCF_AGL(IQ) > 0) NTSL = NTSL + 1
   CALL INITIALIZE_VTKHDF_SLICE(FILENAME,SLCFNAME,IQ,NTSL)
   HDF_SLCF_COUNTER(IQ)=0
ENDDO UNIQUE_LOOPF

IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_SLICE()

END SUBROUTINE INITIALIZE_SLCF_VTKHDF_SERIES





!> \brief Append one output time of Smoke3D data to the VTKHDF file
!>
!> \param T Current simulation time (s)
!> \param DT Current time step size (s)

SUBROUTINE DUMP_SMOKE3D_VTKHDF(T,DT)

USE HDF5
USE ISOSMOKE, ONLY: SMOKE3D_TO_FILE
REAL(EB), INTENT(IN) :: T,DT
INTEGER :: NM
INTEGER  :: I,J,K,N,IFACT,NC,NP,I1=0,J1=0,K1=0,I2,J2,K2,NX,NY,NZ
REAL(FB) :: DXX
REAL(EB), POINTER, DIMENSION(:,:,:) :: FF
REAL(FB), ALLOCATABLE, DIMENSION(:) :: QQ_PACK
INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: QQ_PACK_INT
INTEGER, ALLOCATABLE, DIMENSION(:) :: WRITE_NM
INTEGER :: N_LOCAL, N_LOCAL_TOT, IOFF
TYPE(SMOKE3D_TYPE), POINTER :: S3
REAL(FB) :: FACTOR,VAL_FDS,VAL_SMV,TEMP_MIN
INTEGER(HID_T) :: DSET_ID       ! Identifiers
INTEGER     ::   ERROR ! Error flag
REAL(FB), DIMENSION(1) :: VTK_T
INTEGER(HSIZE_T) :: BASE_OFFSET
INTEGER(HSIZE_T), DIMENSION(1) :: CDIM1,DDIM1,MDIM
INTEGER(IB32), DIMENSION(1) :: INDATA(1)

! Miscellaneous settings

DRY   = .FALSE.

IF (N_SMOKE3D > 0) THEN
   !WRITE(FILENAME,'(A,A,A)') "",TRIM(VTK_DIR)//TRIM(CHID),'_SM3D.vtkhdf'
ELSE
   RETURN
ENDIF

ALLOCATE(WRITE_NM(MAX(1,UPPER_MESH_INDEX-LOWER_MESH_INDEX+1)))
ALLOCATE(QQ_PACK_INT(1024))

CDIM1=(/1_HSIZE_T/)
DDIM1=(/0_HSIZE_T/)
MDIM = (/H5S_UNLIMITED_F/)
! Write time step Information
CALL PARALLEL_INIT_F32(HDF_SM3D_G5, TRIM("Values"), HDF_SM3D_CRP_LIST, 1, CDIM1, DDIM1,&
   DSET_ID, HDF_SM3D_PLIST_ID, MDIM) ! Values = Timesteps
VTK_T = REAL(T,FB)
CALL APPEND_RANK1_DATASET_F32(DSET_ID,HDF_SM3D_PLIST_ID,VTK_T,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
INDATA=(/0_IB32/)
CALL PARALLEL_INIT_I32(HDF_SM3D_G5, TRIM("ConnectivityIdOffsets"), HDF_SM3D_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_SM3D_PLIST_ID, MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SM3D_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_SM3D_G5, TRIM("PointOffsets"), HDF_SM3D_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_SM3D_PLIST_ID, MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SM3D_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_SM3D_G5, TRIM("CellOffsets"), HDF_SM3D_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_SM3D_PLIST_ID, MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SM3D_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_SM3D_G5, TRIM("PartOffsets"), HDF_SM3D_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_SM3D_PLIST_ID, MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SM3D_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_SM3D_G5, TRIM("NumberOfParts"), HDF_SM3D_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_SM3D_PLIST_ID, MDIM)
INDATA=(/NMESHES/)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SM3D_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

HDF_SM3D_COUNTER = HDF_SM3D_COUNTER + 1
CDIM1=(/1_HSIZE_T/)
CALL ADD_ATTRIBUTE_INT(HDF_SM3D_G5,CDIM1,"NSteps",HDF_SM3D_COUNTER)

!CALL VTK_DSET_RELEASE(DSET_ID,ERROR)
!CALL TEST_VTKHDF_EXTEND(TRIM("Values"),HDF_SM3D_G5,DSET_ID,HDF_SM3D_PLIST_ID)

! This set of commented code reads in the current time step information
! and writes it out for debugging purposes
!CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
!CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE, CDIM2, MDIM2, ERROR)
!ALLOCATE(DATA(CDIM2(1)))
!CALL H5DREAD_F(DSET_ID, H5T_IEEE_F32LE, DATA, CDIM2, ERROR)
!WRITE(*,*) "DATA READ", CDIM2, DATA
!DEALLOCATE(DATA)
!CALL H5SCLOSE_F(DATASPACE,ERROR)

!Extend the dataset.
!VTK_T = REAL(T,FB)
!CALL APPEND_RANK1_DATASET_F32(DSET_ID,HDF_SM3D_PLIST_ID,VTK_T,1)

! This set of commented code reads in the current time step information
! and writes it out for debugging purposes
!CALL H5DGET_SPACE_F(DSET_ID, DATASPACE, ERROR)
!CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE, CDIM2, MDIM2, ERROR)
!ALLOCATE(DATA(CDIM2(1)))
!CALL H5DREAD_F(DSET_ID, H5T_IEEE_F32LE, DATA, CDIM2, ERROR)
!WRITE(*,*) "DATA WRITE", CDIM2, DATA
!DEALLOCATE(DATA)
!CALL H5SCLOSE_F(DATASPACE,ERROR)
!CALL VTK_DSET_RELEASE(DSET_ID, ERROR)


! Write data

DATA_FILE_LOOP: DO N=1,N_SMOKE3D
   CALL POINT_TO_MESH(LOWER_MESH_INDEX)
   S3 => SMOKE3D_FILE(N)
   !WRITE(*,*) MY_RANK, N, S3%QUANTITY_INDEX
   IF (S3%QUANTITY_INDEX==0) CYCLE
   
   CALL EXTEND_SMOKE3D_VTKHDF(S3%SMOKEVIEW_LABEL(1:30), BASE_OFFSET)
   !WRITE(*,*) MY_RANK, N, TRIM(S3%SMOKEVIEW_LABEL(1:30)), BASE_OFFSET
   
   ! Write point data offsets.
   CALL PARALLEL_INIT_I32(HDF_SM3D_G7, TRIM(S3%SMOKEVIEW_LABEL(1:30)), HDF_SM3D_CRP_LIST, 1,&
      CDIM1, DDIM1, DSET_ID, HDF_SM3D_PLIST_ID, MDIM)
   INDATA=INT((/BASE_OFFSET/),IB32)
   CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SM3D_PLIST_ID,INDATA,1)
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)
   
   ! Loop over this rank's own meshes only, packing everything into one buffer.  The write
   ! itself is a single collective call made below, instead of one per mesh in the global
   ! mesh list on every rank.

   N_LOCAL = 0
   N_LOCAL_TOT = 0
   ALL_MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      DXX   = REAL(DX(1),FB)
      S3 => SMOKE3D_FILE(N)
      FF   => WORK3
      ! Current mesh is owned by this rank, and is contained in this mesh, dump real data
      ! Obtain Smoke3D output at cell centers
      CALL GET_SMOKE3D_QQ(S3,T,DT,NM,FF,QQ)

      ! Pack the data into a 1-D array and send to the routine that writes the file for Smokeview

      I2 = MESHES(NM)%IBAR
      J2 = MESHES(NM)%JBAR
      K2 = MESHES(NM)%KBAR
      NX = I2 + 1 - I1
      NY = J2 + 1 - J1
      NZ = K2 + 1 - K1
      ALLOCATE(QQ_PACK(MESHES(NM)%NP))

      IFACT = 1
      DO K=0,KBP1-1
         DO J=0,JBP1-1
            DO I=0,IBP1-1
               QQ_PACK(IFACT) = REAL(QQ(I,J,K,1))
               IFACT = IFACT + 1
            ENDDO
         ENDDO
      ENDDO
      NP = MESHES(NM)%NP
      NC = MESHES(NM)%NC

      IF (N_LOCAL_TOT+NP>SIZE(QQ_PACK_INT)) CALL GROW_QQ_PACK_INT(N_LOCAL_TOT+NP)
      IOFF = N_LOCAL_TOT
      IF (S3%DISPLAY_TYPE=='GAS') THEN

         FACTOR=-REAL(S3%MASS_EXTINCTION_COEFFICIENT,FB)*DXX
         DO I=1,NP
            VAL_FDS = MAX(0.0_FB,QQ_PACK(I))
            VAL_SMV = 254*(1.0_FB-EXP(FACTOR*VAL_FDS)) !-127
            QQ_PACK_INT(IOFF+I) = INT(NINT(VAL_SMV),IB8) !NINT(VAL_SMV)
         ENDDO


      ELSEIF (S3%DISPLAY_TYPE=='FIRE') THEN

         DO I=1,NP
            VAL_FDS = MIN(HRRPUV_MAX_SMV,MAX(0._FB,QQ_PACK(I)))
            VAL_SMV = 254*(VAL_FDS/HRRPUV_MAX_SMV) !-127
            QQ_PACK_INT(IOFF+I) = INT(NINT(VAL_SMV),IB8)
         ENDDO

      ELSEIF (S3%DISPLAY_TYPE=='TEMPERATURE') THEN

         TEMP_MIN = REAL(TMPA-TMPM,FB)
         DO I=1,NP
            VAL_FDS = MIN(TEMP_MAX_SMV,MAX(TEMP_MIN,QQ_PACK(I)))
            VAL_SMV = 254*((VAL_FDS-TEMP_MIN)/(TEMP_MAX_SMV-TEMP_MIN)) !-127
            QQ_PACK_INT(IOFF+I) = INT(NINT(VAL_SMV),IB8)
         ENDDO

      ELSE
         DO I=1,NP
            QQ_PACK_INT(IOFF+I) = INT(NINT(QQ_PACK(I)),IB8)
         ENDDO
      ENDIF
      DEALLOCATE(QQ_PACK)
      N_LOCAL = N_LOCAL + 1
      WRITE_NM(N_LOCAL) = NM
      N_LOCAL_TOT = N_LOCAL_TOT + NP
   ENDDO ALL_MESH_LOOP

   CALL ADD_DATA_TO_SMOKE3D_VTKHDF_MULTI(S3%SMOKEVIEW_LABEL(1:30),QQ_PACK_INT,&
      BASE_OFFSET,N_LOCAL,WRITE_NM)
ENDDO DATA_FILE_LOOP

DEALLOCATE(WRITE_NM)
DEALLOCATE(QQ_PACK_INT)

CONTAINS

!> \brief Grow the packing buffer, preserving what has already been written into it

SUBROUTINE GROW_QQ_PACK_INT(NREQ)
INTEGER, INTENT(IN) :: NREQ
INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: TMP
ALLOCATE(TMP(2*NREQ))
TMP(1:SIZE(QQ_PACK_INT)) = QQ_PACK_INT
CALL MOVE_ALLOC(TMP,QQ_PACK_INT)
END SUBROUTINE GROW_QQ_PACK_INT

END SUBROUTINE DUMP_SMOKE3D_VTKHDF


! \brief Write contour slices, Plot3D data, or 3d slices to a file
!>
!> \param T Current simulation time (s)
!> \param DT Current time step size (s)
!> \param NM Mesh number
!> \param IFRMT VTK 3D slice (IFRMT=0) VTK 2D slice (IFRMT=1)

SUBROUTINE DUMP_SLCF_VTK(T,DT,IFRMT)

USE MEMORY_FUNCTIONS, ONLY: RE_ALLOCATE_STRINGS
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
USE TRAN, ONLY : GET_IJK
USE ISOSMOKE, ONLY: SLICE_TO_RLEFILE
INTEGER :: NM,IFRMT
REAL(EB), INTENT(IN) :: T,DT
INTEGER :: I,J,K,I1,I2,J1,J2,K1,K2,IQ,IND,II, &
           IFACT,JFACT,KFACT,NX,NY,NZ,NTSL,NTSL_LOCAL,SLICEIND
REAL(EB), POINTER, DIMENSION(:,:,:) :: B,S,QUANTITY
LOGICAL :: VTK3D,SL3D
LOGICAL :: AGL_TERRAIN_SLICE,CC_CELL_CENTERED,CC_INTERP2FACES
CHARACTER(200) :: SLCFNAME,QTY
INTEGER(IB32), DIMENSION(1:NMESHES) :: NCELLS, NPOINTS
INTEGER(HID_T) :: DSET_ID       ! Identifiers
INTEGER     ::   ERROR ! Error flag
REAL(FB), DIMENSION(1) :: VTK_T
INTEGER(HSIZE_T) :: BASE_OFFSET_PTS, BASE_OFFSET_CELLS
INTEGER(HSIZE_T), DIMENSION(1) :: CDIM1(1), DDIM1(1), MDIM(1)
INTEGER(IB32), DIMENSION(1) :: INDATA(1)
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: SLICE_LOOKUP
INTEGER, ALLOCATABLE, DIMENSION(:) :: WRITE_NM
REAL(FB), ALLOCATABLE, DIMENSION(:) :: LOCAL_DATA
INTEGER :: N_LOCAL, N_LOCAL_TOT

SELECT CASE(IFRMT)
   CASE(0) ; VTK3D=.TRUE.
   CASE(1) ; VTK3D=.FALSE.
END SELECT

! The cell-corner averaging weights (B and S below) depend only on the mesh topology, so
! build them once per mesh here rather than once per (unique slice, quantity, mesh).  The
! same goes for the map from a (unique slice, quantity) pair to this mesh's slice index,
! which used to be a linear string search repeated inside the innermost loop.

CALL BUILD_SLICE_WEIGHTS_AND_LOOKUP

ALLOCATE(WRITE_NM(MAX(1,UPPER_MESH_INDEX-LOWER_MESH_INDEX+1)))
ALLOCATE(LOCAL_DATA(1024))

! Get time string for filename

NTSL = 0
BASE_OFFSET_PTS=0
BASE_OFFSET_CELLS=0
UNIQUE_LOOPF: DO IQ=1,MESHES(1)%N_UNIQUE_SLCF
   SL3D = MESHES(1)%UNIQUE_SLICE_IS_SL3D(IQ)
   IF (SL3D.AND..NOT.VTK3D) CYCLE UNIQUE_LOOPF
   IF (.NOT.SL3D.AND.VTK3D) CYCLE UNIQUE_LOOPF
   SLCFNAME = MESHES(1)%UNIQUE_SLICE_NAMES(IQ)
   IF (MESHES(1)%UNIQUE_SLCF_AGL(IQ) > 0) NTSL = NTSL + 1
   
   ! Write time step Information
   DDIM1=1
   CDIM1=0
   MDIM = (/H5S_UNLIMITED_F/)
   INDATA=0
   CALL PARALLEL_INIT_F32(HDF_SLCF_G5(IQ), TRIM("Values"), HDF_SLCF_CRP_LIST(IQ), 1, DDIM1,CDIM1,&
      DSET_ID, HDF_SLCF_PLIST_ID(IQ), MDIM) ! Values = Timesteps
      VTK_T = REAL(T,FB)
   CALL APPEND_RANK1_DATASET_F32(DSET_ID,HDF_SLCF_PLIST_ID(IQ),VTK_T,1)
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)
   
   ! Write step offsets. Assuming no offsets for meshing at this time.
   CALL PARALLEL_INIT_I32(HDF_SLCF_G5(IQ), TRIM("ConnectivityIdOffsets"), HDF_SLCF_CRP_LIST(IQ), 1,&
      DDIM1, CDIM1, DSET_ID, HDF_SLCF_PLIST_ID(IQ), MDIM)
   CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SLCF_PLIST_ID(IQ),INDATA,1)
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

   ! Write step offsets. Assuming no offsets for meshing at this time.
   CALL PARALLEL_INIT_I32(HDF_SLCF_G5(IQ), TRIM("PointOffsets"), HDF_SLCF_CRP_LIST(IQ), 1,&
      DDIM1, CDIM1, DSET_ID, HDF_SLCF_PLIST_ID(IQ), MDIM)
   CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SLCF_PLIST_ID(IQ),INDATA,1)
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

   ! Write step offsets. Assuming no offsets for meshing at this time.
   CALL PARALLEL_INIT_I32(HDF_SLCF_G5(IQ), TRIM("CellOffsets"), HDF_SLCF_CRP_LIST(IQ), 1,&
      DDIM1, CDIM1, DSET_ID, HDF_SLCF_PLIST_ID(IQ), MDIM)
   CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SLCF_PLIST_ID(IQ),INDATA,1)
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

   ! Write step offsets. Assuming no offsets for meshing at this time.
   CALL PARALLEL_INIT_I32(HDF_SLCF_G5(IQ), TRIM("PartOffsets"), HDF_SLCF_CRP_LIST(IQ), 1,&
      DDIM1, CDIM1, DSET_ID, HDF_SLCF_PLIST_ID(IQ), MDIM)
   CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SLCF_PLIST_ID(IQ),INDATA,1)
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

   ! Write step offsets. Assuming no offsets for meshing at this time.
   INDATA=(/NMESHES/)
   CALL PARALLEL_INIT_I32(HDF_SLCF_G5(IQ), TRIM("NumberOfParts"), HDF_SLCF_CRP_LIST(IQ), 1,&
      DDIM1, CDIM1, DSET_ID, HDF_SLCF_PLIST_ID(IQ), MDIM)
   CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SLCF_PLIST_ID(IQ),INDATA,1)
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

   HDF_SLCF_COUNTER(IQ) = HDF_SLCF_COUNTER(IQ) + 1
   CALL ADD_ATTRIBUTE_INT(HDF_SLCF_G5(IQ),DDIM1,"NSteps",HDF_SLCF_COUNTER(IQ))
   
   !CYCLE UNIQUE_LOOPF
   
   ALL_SLCF_LOOP: DO II=1,MESHES(1)%N_SLCF_VTK
      ! If this slice is not part of this unique location, cycle
      !SL => SLICE(II)
      !IF (TRIM(SL%SLCF_NAME)/=TRIM(SLCFNAME)) CYCLE ALL_SLCF_LOOP
      IF(TRIM(MESHES(1)%ALL_SLICE_NAMES(II))/=TRIM(MESHES(1)%UNIQUE_SLICE_NAMES(IQ))) CYCLE ALL_SLCF_LOOP
      QTY = MESHES(1)%ALL_SLICE_QUANTITIES(II)
      
      CALL EXTEND_SLICE_VTKHDF(QTY, IQ, BASE_OFFSET_PTS,NCELLS,NPOINTS)
      !WRITE(*,*) MY_RANK, TRIM(MESHES(1)%ALL_SLICE_NAMES(II)), TRIM(QTY), BASE_OFFSET_PTS, NCELLS, NPOINTS
      
      ! Write point data offsets.
      CALL PARALLEL_INIT_I32(HDF_SLCF_G7(IQ), TRIM(QTY), HDF_SLCF_CRP_LIST(IQ), 1,&
         DDIM1, CDIM1, DSET_ID, HDF_SLCF_PLIST_ID(IQ), MDIM)
      INDATA = INT((/BASE_OFFSET_PTS/),IB32)
      CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_SLCF_PLIST_ID(IQ),INDATA,1)
      CALL VTK_DSET_RELEASE(DSET_ID,ERROR)
      
      ! Gather this rank's contribution for all of the meshes it owns, then hand it to a
      ! single collective write.  The previous version opened the dataset and issued a
      ! collective write once per mesh in the global mesh list, on every rank, which made
      ! the cost of a slice dump grow with NMESHES on every rank.

      N_LOCAL = 0
      N_LOCAL_TOT = 0
      LOCAL_MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

         IF (MESHES(NM)%EMPTY_UNIQUE_SLICE(IQ)) CYCLE LOCAL_MESH_LOOP
         CALL POINT_TO_MESH(NM)
         
         ! Current mesh is owned by this rank, and is contained in this mesh, dump real data
         ! B and S (the cell-corner averaging weights) were built once per mesh above.

         B => WORK1
         S => WORK2

         ! If sprinkler diagnostic on, pre-compute various PARTICLE flux output

         IF (SLCF_PARTICLE_FLUX) CALL COMPUTE_PARTICLE_FLUXES ! TODO Not sure what we need for VTK here

         ! Determine slice or Plot3D indicies

         QUANTITY=>WORK7
         
         SLICEIND = SLICE_LOOKUP(IQ,II,NM)
         ! If this slice and quantity combo is not present in this mesh, contribute nothing
         IF (SLICEIND < 0) CYCLE LOCAL_MESH_LOOP
         
         ! Pack data for parallel write
         SL => SLICE(SLICEIND)
         DRY  = SL%DRY   ! GAS_PHASE_OUTPUT reads this, as it does in DUMP_SLCF
         IND  = SL%INDEX
         I1  = SL%I1
         I2  = SL%I2
         J1  = SL%J1
         J2  = SL%J2
         K1  = SL%K1
         K2  = SL%K2
         AGL_TERRAIN_SLICE = SL%TERRAIN_SLICE
         CC_CELL_CENTERED  = SL%CELL_CENTERED
         CC_INTERP2FACES   = .FALSE.
         IF (.NOT.CC_CELL_CENTERED .AND. TRIM(SL%SLICETYPE)/='STRUCTURED') CC_INTERP2FACES = .TRUE.

         ! Evaluate the quantity and average it onto the slice nodes, faces or edges.  This
         ! is the same routine DUMP_SLCF uses, so the .sf and .vtkhdf values agree.

         NTSL_LOCAL = NTSL - 1
         CALL GET_SLICE_QUANTITY(T,DT,NM,IND,I1,I2,J1,J2,K1,K2,SL%Y_INDEX,SL%Z_INDEX,SL%PART_INDEX,SL%VELO_INDEX,0,&
                                 SL%REAC_INDEX,AGL_TERRAIN_SLICE,CC_CELL_CENTERED,CC_INTERP2FACES,NTSL_LOCAL,&
                                 B,S,QUANTITY,QQ,1)

         NX = I2 + 1 - I1
         NY = J2 + 1 - J1
         NZ = K2 + 1 - K1

         IF (NPOINTS(NM)<=0) CYCLE LOCAL_MESH_LOOP

         IF (N_LOCAL_TOT+NPOINTS(NM)>SIZE(LOCAL_DATA)) CALL GROW_LOCAL_DATA(N_LOCAL_TOT+NPOINTS(NM))
         LOCAL_DATA(N_LOCAL_TOT+1:N_LOCAL_TOT+NPOINTS(NM)) = 0._FB
         DO K = K1, K2
            KFACT = (K-K1)*NY*NX
            DO J = J1, J2
               JFACT = (J-J1)*NX
               DO I = I1, I2
                  IFACT = (I-I1)
                  IF (1+IFACT+JFACT+KFACT>NPOINTS(NM)) CYCLE
                  LOCAL_DATA(N_LOCAL_TOT+1+IFACT+JFACT+KFACT) = QQ(I,J,K,1)
               ENDDO
            ENDDO
         ENDDO
         N_LOCAL = N_LOCAL + 1
         WRITE_NM(N_LOCAL) = NM
         N_LOCAL_TOT = N_LOCAL_TOT + NPOINTS(NM)

      ENDDO LOCAL_MESH_LOOP

      CALL WRITE_VTKHDF_SLICE_DATA_MULTI(TRIM(QTY),HDF_SLCF_G4(IQ),HDF_SLCF_CRP_LIST(IQ),&
         HDF_SLCF_PLIST_ID(IQ),NPOINTS,BASE_OFFSET_PTS,N_LOCAL,WRITE_NM,LOCAL_DATA)
      ! Write fake data in processes that have less meshes than max process
      !ALLOCATE(QQ_PACK(0))
      !N_WRITTEN = MESHES_PER_PROCESS(MY_RANK)
      !DO NM=1,MAXVAL(MESHES_PER_PROCESS)
      !   IF (N_WRITTEN < MAXVAL(MESHES_PER_PROCESS)) THEN
      !      !WRITE(*,*) "LINE1", QTY, NM, NCELLS
      !      !WRITE(*,*) "LINE2", NPOINTS
      !      !WRITE(*,*) "LINE3", BASE_OFFSET_PTS
      !      
      !      !CALL WRITE_VTKHDF_SLICE_DATA_FILE_NOOPEN(QTY,NM,NCELLS,NPOINTS,QQ_PACK,&
      !      !   HDF_SLCF_PLIST_ID(IQ),HDF_SLCF_CRP_LIST(IQ),HDF_SLCF_G4(IQ),BASE_OFFSET_PTS,.TRUE.)
      !      N_WRITTEN = N_WRITTEN + 1
      !   ENDIF
      !ENDDO
      !DEALLOCATE(QQ_PACK)
      ! Close VTKHDF interface
      !IF (.NOT.VTK_KEEPOPEN) THEN CALL CLOSE_VTKHDF(FILE_ID, GROUP_ID1,GROUP_ID2,GROUP_ID3,GROUP_ID4)
   ENDDO ALL_SLCF_LOOP
   
ENDDO UNIQUE_LOOPF

IF (ALLOCATED(SLICE_LOOKUP)) DEALLOCATE(SLICE_LOOKUP)
IF (ALLOCATED(WRITE_NM))     DEALLOCATE(WRITE_NM)
IF (ALLOCATED(LOCAL_DATA))   DEALLOCATE(LOCAL_DATA)

CONTAINS

!> \brief Grow the packing buffer, preserving what has already been written into it

SUBROUTINE GROW_LOCAL_DATA(N)
INTEGER, INTENT(IN) :: N
REAL(FB), ALLOCATABLE, DIMENSION(:) :: TMP
ALLOCATE(TMP(2*N))
TMP(1:SIZE(LOCAL_DATA)) = LOCAL_DATA
CALL MOVE_ALLOC(TMP,LOCAL_DATA)
END SUBROUTINE GROW_LOCAL_DATA


!> \brief Build the per-mesh cell-corner averaging weights and the slice index map
!>
!> B is 1 in every cell that takes part in the 8-cell corner average and 0 otherwise; S
!> is the reciprocal of the sum of the eight B values meeting at a corner.  Both depend
!> only on mesh topology, so they are built once per mesh per dump instead of once per
!> (unique slice, quantity, mesh) as before.  SLICE_LOOKUP replaces the linear string
!> search that used to run in the innermost loop.

SUBROUTINE BUILD_SLICE_WEIGHTS_AND_LOOKUP

INTEGER :: NML,ICL,IL,JL,KL,IQL,IIL,SIQL,NOML,IIOL,JJOL,KKOL,ICOL
REAL(EB) :: BSUML
REAL(EB), POINTER, DIMENSION(:,:,:) :: BL,SL_W
TYPE (MESH_TYPE), POINTER :: M2L

ALLOCATE(SLICE_LOOKUP(MAX(1,MESHES(1)%N_UNIQUE_SLCF),MAX(1,MESHES(1)%N_SLCF_VTK),&
                      LOWER_MESH_INDEX:UPPER_MESH_INDEX))
SLICE_LOOKUP = -1

MESH_SETUP_LOOP: DO NML=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NML)

   BL => WORK1
   BL = 1._EB

   DO ICL=1,CELL_COUNT(NML)
      IF (CELL(ICL)%SOLID) BL(CELL(ICL)%I,CELL(ICL)%J,CELL(ICL)%K) = 0._EB
      IF (CELL(ICL)%EXTERIOR) THEN
         IF (CELL(ICL)%EXTERIOR_EDGE) THEN
            BL(CELL(ICL)%I,CELL(ICL)%J,CELL(ICL)%K) = 0._EB
         ELSE
            CALL SEARCH_OTHER_MESHES(XC(CELL(ICL)%I),YC(CELL(ICL)%J),ZC(CELL(ICL)%K),NOML,IIOL,JJOL,KKOL)
            IF (NOML==0) THEN
               BL(CELL(ICL)%I,CELL(ICL)%J,CELL(ICL)%K) = 0._EB
            ELSE
               M2L => MESHES(NOML)
               ICOL = M2L%CELL_INDEX(IIOL,JJOL,KKOL)
               IF (M2L%CELL(ICOL)%SOLID) BL(CELL(ICL)%I,CELL(ICL)%J,CELL(ICL)%K) = 0._EB
            ENDIF
         ENDIF
      ENDIF
   ENDDO

   SL_W => WORK2
   SL_W = 0._EB

   DO KL=0,KBAR
      DO JL=0,JBAR
         DO IL=0,IBAR
            BSUML = BL(IL,JL,KL)+BL(IL+1,JL+1,KL+1)+BL(IL+1,JL,KL)+BL(IL,JL+1,KL)+ &
                    BL(IL,JL,KL+1)+BL(IL+1,JL+1,KL)+BL(IL+1,JL,KL+1)+BL(IL,JL+1,KL+1)
            IF (BSUML>0._EB) SL_W(IL,JL,KL) = 1._EB/BSUML
         ENDDO
      ENDDO
   ENDDO

   DO IQL=1,MESHES(1)%N_UNIQUE_SLCF
      DO IIL=1,MESHES(1)%N_SLCF_VTK
         IF (TRIM(MESHES(1)%ALL_SLICE_NAMES(IIL))/=TRIM(MESHES(1)%UNIQUE_SLICE_NAMES(IQL))) CYCLE
         DO SIQL=1,MESHES(1)%N_SLCF_VTK
            IF (TRIM(MESHES(NML)%SLICE(SIQL)%SLCF_NAME)/=TRIM(MESHES(1)%UNIQUE_SLICE_NAMES(IQL))) CYCLE
            IF (TRIM(MESHES(NML)%SLICE(SIQL)%SMOKEVIEW_LABEL(1:30))/= &
                TRIM(MESHES(1)%ALL_SLICE_QUANTITIES(IIL))) CYCLE
            SLICE_LOOKUP(IQL,IIL,NML) = SIQL
            EXIT
         ENDDO
      ENDDO
   ENDDO

ENDDO MESH_SETUP_LOOP

END SUBROUTINE BUILD_SLICE_WEIGHTS_AND_LOOKUP

END SUBROUTINE DUMP_SLCF_VTK


! \brief Dump boundary quantities into CHID_nn.bf file
!> \param T Current simulation time (s)
!> \param DT Current time step size (s)

SUBROUTINE DUMP_BNDF_VTKHDF(T,DT)
INTEGER(HID_T) :: DSET_ID    ! Dataset identifiers
REAL(EB), INTENT(IN) :: T,DT
INTEGER :: NF,IND,IC,IW,IP,NC,I1,I2,J1,J2,K1,K2
INTEGER :: IFACT,NFACES_CUTCELLS,NVERTS_CUTCELLS,NVERTS,NFACES
INTEGER :: NM
TYPE(PATCH_TYPE), POINTER :: PA
INTEGER :: NMNM, NM1, NM2, ERROR, PA_NCELLS, PA_NPOINTS
INTEGER(IB32), DIMENSION(1:2*NMESHES) :: NCELLS, NPOINTS
INTEGER, DIMENSION(1:NMESHES) :: NCELLS_OFFSET, NPOINTS_OFFSET
INTEGER :: NCELLS_ACCUM, NCELLS_START
INTEGER :: NPOINTS_ACCUM, NPOINTS_START
INTEGER :: N_BLK, N_BLK_TOT
INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCATIONS, BLK_START, BLK_COUNT
REAL(FB), ALLOCATABLE, DIMENSION(:) :: QQ_PACK, ALL_DATA, BLK_DATA
INTEGER(HSIZE_T), DIMENSION(1) :: CDIM1(1), DDIM1(1), MDIM(1) !, CDIM2(2), DDIM2(2)
REAL(FB), DIMENSION(1) :: VTK_T
INTEGER(HSIZE_T) :: BASE_OFFSET_CELLS, BASE_OFFSET_PTS
INTEGER(IB32), DIMENSION(1) :: INDATA(1)

ALLOCATE(BLK_START(64)) ; ALLOCATE(BLK_COUNT(64)) ; ALLOCATE(BLK_DATA(1024))

CDIM1=1
DDIM1=0
MDIM = (/H5S_UNLIMITED_F/)
INDATA=(/0_IB32/)
! Write time step Information
CALL PARALLEL_INIT_F32(HDF_BNDF_G5, TRIM("Values"), HDF_BNDF_CRP_LIST, 1, CDIM1, DDIM1,&
   DSET_ID, HDF_BNDF_PLIST_ID, MDIM) ! Values = Timesteps
VTK_T = REAL(T,FB)
CALL APPEND_RANK1_DATASET_F32(DSET_ID,HDF_BNDF_PLIST_ID,VTK_T,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_BNDF_G5, TRIM("ConnectivityIdOffsets"), HDF_BNDF_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_BNDF_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_BNDF_G5, TRIM("PointOffsets"), HDF_BNDF_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_BNDF_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_BNDF_G5, TRIM("CellOffsets"), HDF_BNDF_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_BNDF_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_BNDF_G5, TRIM("PartOffsets"), HDF_BNDF_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_BNDF_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

! Write step offsets. Assuming no offsets for meshing at this time.
CALL PARALLEL_INIT_I32(HDF_BNDF_G5, TRIM("NumberOfParts"), HDF_BNDF_CRP_LIST, 1,&
   CDIM1, DDIM1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM)
INDATA=(/2*NMESHES/)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_BNDF_PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

HDF_BNDF_COUNTER = HDF_BNDF_COUNTER + 1
CDIM1=1
CALL ADD_ATTRIBUTE_INT(HDF_BNDF_G5,CDIM1,"NSteps",HDF_BNDF_COUNTER)

! Extend bndf datasize
CALL EXTEND_BNDF_VTKHDF(BASE_OFFSET_CELLS, BASE_OFFSET_PTS)

! Write point and cell data offsets.
DO NF=1,N_BNDF
   BF => BOUNDARY_FILE(NF)
   IF (BF%CELL_CENTERED) THEN
      CALL PARALLEL_INIT_I32(HDF_BNDF_G6, TRIM(BF%SMOKEVIEW_LABEL(1:30)), HDF_BNDF_CRP_LIST, 1,&
         CDIM1, DDIM1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM)
      INDATA=INT((/BASE_OFFSET_CELLS/),IB32)
      CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_BNDF_PLIST_ID,INDATA,1)
   ELSE
      CALL PARALLEL_INIT_I32(HDF_BNDF_G7, TRIM(BF%SMOKEVIEW_LABEL(1:30)), HDF_BNDF_CRP_LIST, 1,&
         CDIM1, DDIM1, DSET_ID, HDF_BNDF_PLIST_ID, MDIM)
      INDATA=INT((/BASE_OFFSET_PTS/),IB32)
      CALL APPEND_RANK1_DATASET_I32(DSET_ID,HDF_BNDF_PLIST_ID,INDATA,1)
   ENDIF
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)
ENDDO

NCELLS = NCELLS_VTK
NPOINTS = NPOINTS_VTK

! Prefix sums of the per-part sizes.  Parts are laid out two per mesh: the odd part holds
! the OBST patches, the even part holds the GEOM (cut-cell) faces.

NPOINTS_ACCUM = 0
NCELLS_ACCUM = 0
DO NMNM=1,NMESHES
   NM1 = 2*NMNM-1
   NM2 = 2*NMNM
   NPOINTS_OFFSET(NMNM) = NPOINTS_ACCUM
   NCELLS_OFFSET(NMNM) = NCELLS_ACCUM
   NPOINTS_ACCUM = NPOINTS_ACCUM + NPOINTS(NM1) + NPOINTS(NM2)
   NCELLS_ACCUM = NCELLS_ACCUM + NCELLS(NM1) + NCELLS(NM2)
ENDDO

! One collective write per boundary quantity, covering every mesh this rank owns.  The
! previous version opened the dataset and wrote once per (mesh, quantity), then padded
! with dummy collective writes so that ranks owning fewer meshes stayed in step.

DATA_FILE_LOOP: DO NF=1,N_BNDF

   BF => BOUNDARY_FILE(NF)
   PY => PROPERTY(BF%PROP_INDEX)
   IND  = ABS(BF%INDEX)
   N_BLK = 0
   N_BLK_TOT = 0

   MESH_LOOP_HDF: DO NMNM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      NM1 = 2*NMNM-1
      NM2 = 2*NMNM
      NM = NMNM
      CALL POINT_TO_MESH(NM)
      NPOINTS_START = NPOINTS_OFFSET(NMNM)
      NCELLS_START = NCELLS_OFFSET(NMNM)
      NPOINTS_ACCUM = NPOINTS_START
      NCELLS_ACCUM = NCELLS_START
      ! Write OBST boundary data
      IF (MESHES(NM)%N_PATCH>0) THEN
         ! Fill ALL_DATA with values from all patches in mesh NM
         IF (BF%CELL_CENTERED) ALLOCATE(ALL_DATA(NCELLS(NM1)))
         IF (.NOT.BF%CELL_CENTERED) ALLOCATE(ALL_DATA(NPOINTS(NM1)))
         PATCH_LOOP1: DO IP=1,N_PATCH
            PA => PATCH(IP)
            IF (PA%OBST_INDEX<=0) CYCLE PATCH_LOOP1
            CALL GET_PA_POINTS_AND_CELLS(PA,PA_NPOINTS,PA_NCELLS)
            CALL PACK_VTK_BNDF(PA,BF,IND,PP,PPN,QQ_PACK)
            IF (BF%CELL_CENTERED) THEN
               ALL_DATA(NCELLS_ACCUM-NCELLS_START+1:NCELLS_ACCUM-NCELLS_START+PA_NCELLS) = QQ_PACK
            ELSE
               ALL_DATA(NPOINTS_ACCUM-NPOINTS_START+1:NPOINTS_ACCUM-NPOINTS_START+PA_NPOINTS) = QQ_PACK
            ENDIF
            DEALLOCATE(QQ_PACK)
            NCELLS_ACCUM = NCELLS_ACCUM + PA_NCELLS
            NPOINTS_ACCUM = NPOINTS_ACCUM + PA_NPOINTS
         ENDDO PATCH_LOOP1

         ! Stage ALL_DATA for the collective write below
         IF (BF%CELL_CENTERED) THEN
            CALL ADD_BLOCK(NCELLS_START + INT(BASE_OFFSET_CELLS), NCELLS(NM1), ALL_DATA)
         ELSE
            CALL ADD_BLOCK(NPOINTS_START + INT(BASE_OFFSET_PTS), NPOINTS(NM1), ALL_DATA)
         ENDIF
         DEALLOCATE(ALL_DATA)
      ENDIF
      
      ! Write GEOM boundary data
      IF (MESHES(NM)%N_INTERNAL_CFACE_CELLS>0) THEN
         CALL GET_GEOMSIZES('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS)
         IF (NVERTS>0 .AND. NFACES>0) THEN
            CALL PACK_VTK_GEOM(NM, BF%CELL_CENTERED, BF, ALL_DATA)
            ! Stage ALL_DATA for the collective write below
            IF (BF%CELL_CENTERED) THEN
               CALL ADD_BLOCK(NCELLS_ACCUM + INT(BASE_OFFSET_CELLS), NCELLS(NM2), ALL_DATA)
            ELSE
               CALL ADD_BLOCK(NPOINTS_ACCUM + INT(BASE_OFFSET_PTS), NPOINTS(NM2), ALL_DATA)
            ENDIF
            DEALLOCATE(ALL_DATA)
         ENDIF
      ENDIF
      
   ENDDO MESH_LOOP_HDF

   IF (BF%CELL_CENTERED) THEN
      CALL WRITE_VTKHDF_BNDF_DATA_MULTI(BF%SMOKEVIEW_LABEL(1:30),HDF_BNDF_G2,N_BLK,BLK_START,BLK_COUNT,BLK_DATA)
   ELSE
      CALL WRITE_VTKHDF_BNDF_DATA_MULTI(BF%SMOKEVIEW_LABEL(1:30),HDF_BNDF_G4,N_BLK,BLK_START,BLK_COUNT,BLK_DATA)
   ENDIF

ENDDO DATA_FILE_LOOP

DEALLOCATE(BLK_START) ; DEALLOCATE(BLK_COUNT) ; DEALLOCATE(BLK_DATA)

CONTAINS

!> \brief Stage one contiguous block of boundary data for the collective write

SUBROUTINE ADD_BLOCK(START,COUNT,VALS)
INTEGER, INTENT(IN) :: START,COUNT
REAL(FB), DIMENSION(:), INTENT(IN) :: VALS
INTEGER, ALLOCATABLE, DIMENSION(:) :: ITMP
REAL(FB), ALLOCATABLE, DIMENSION(:) :: RTMP
IF (COUNT<=0) RETURN
IF (N_BLK+1>SIZE(BLK_START)) THEN
   ALLOCATE(ITMP(2*(N_BLK+1))) ; ITMP(1:N_BLK)=BLK_START(1:N_BLK) ; CALL MOVE_ALLOC(ITMP,BLK_START)
   ALLOCATE(ITMP(2*(N_BLK+1))) ; ITMP(1:N_BLK)=BLK_COUNT(1:N_BLK) ; CALL MOVE_ALLOC(ITMP,BLK_COUNT)
ENDIF
IF (N_BLK_TOT+COUNT>SIZE(BLK_DATA)) THEN
   ALLOCATE(RTMP(2*(N_BLK_TOT+COUNT))) ; RTMP(1:N_BLK_TOT)=BLK_DATA(1:N_BLK_TOT)
   CALL MOVE_ALLOC(RTMP,BLK_DATA)
ENDIF
N_BLK = N_BLK + 1
BLK_START(N_BLK) = START
BLK_COUNT(N_BLK) = COUNT
BLK_DATA(N_BLK_TOT+1:N_BLK_TOT+COUNT) = VALS(1:MIN(COUNT,SIZE(VALS)))
IF (COUNT>SIZE(VALS)) BLK_DATA(N_BLK_TOT+SIZE(VALS)+1:N_BLK_TOT+COUNT) = 0._FB
N_BLK_TOT = N_BLK_TOT + COUNT
END SUBROUTINE ADD_BLOCK

   !> \brief Count the points and cells a boundary patch contributes

   !>

   !> \param PA Patch to count

   !> \param PA_NPOINTS Number of points (out)

   !> \param PA_NCELLS Number of cells (out)


   SUBROUTINE GET_PA_POINTS_AND_CELLS(PA,PA_NPOINTS,PA_NCELLS)
      TYPE(PATCH_TYPE), POINTER, INTENT(IN) :: PA
      INTEGER :: L1, L2, N1, N2, NX, NY, NZ
      INTEGER, INTENT(OUT) :: PA_NPOINTS, PA_NCELLS
      ! Fill point data
      SELECT CASE(ABS(PA%IOR))
         CASE(1) ; L1=PA%JG1 ; L2=PA%JG2 ; N1=PA%KG1 ; N2=PA%KG2 ; NX=1; NY=(L2-L1); NZ=(N2-N1)
         CASE(2) ; L1=PA%IG1 ; L2=PA%IG2 ; N1=PA%KG1 ; N2=PA%KG2 ; NX=(L2-L1); NY=1; NZ=(N2-N1)
         CASE(3) ; L1=PA%IG1 ; L2=PA%IG2 ; N1=PA%JG1 ; N2=PA%JG2 ; NX=(L2-L1); NY=(N2-N1); NZ=1
      END SELECT
      PA_NPOINTS = (L2-L1+2)*(N2-N1+2)
      PA_NCELLS = (L2-L1+1)*(N2-N1+1)
   ENDSUBROUTINE GET_PA_POINTS_AND_CELLS

   !> \brief Gather one patch's boundary values into a flat array for writing

   !>

   !> \param PA Patch to gather

   !> \param BF Boundary file the quantity belongs to

   !> \param IND Output quantity index

   !> \param PP Cell centered values on the patch

   !> \param PPN Node averaged values on the patch

   !> \param QQ_PACK Gathered values, allocated here (out)


   SUBROUTINE PACK_VTK_BNDF(PA,BF,IND,PP,PPN,QQ_PACK)
      IMPLICIT NONE
      TYPE(PATCH_TYPE), POINTER, INTENT(IN) :: PA
      TYPE(BOUNDARY_FILE_TYPE), POINTER, INTENT(IN) :: BF
      REAL(FB), POINTER, DIMENSION(:,:), INTENT(IN) :: PP,PPN
      INTEGER, INTENT(IN) :: IND
      INTEGER :: ISUM,I,J,K,L,L1,L2,N,N1,N2
      REAL(FB), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: QQ_PACK
      INTEGER :: NPOINTS, NCELLS

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

      SELECT CASE(ABS(PA%IOR))
         CASE(1) ; L1=PA%JG1 ; L2=PA%JG2 ; N1=PA%KG1 ; N2=PA%KG2;
         CASE(2) ; L1=PA%IG1 ; L2=PA%IG2 ; N1=PA%KG1 ; N2=PA%KG2;
         CASE(3) ; L1=PA%IG1 ; L2=PA%IG2 ; N1=PA%JG1 ; N2=PA%JG2;
      END SELECT
      NPOINTS = (L2-L1+2)*(N2-N1+2)
      NCELLS = (L2-L1+1)*(N2-N1+1)
      IFACT = 1
      IF (.NOT.BF%CELL_CENTERED) THEN
         ALLOCATE(QQ_PACK(NPOINTS))
         DO N = N1-1,N2
            DO L = L1-1,L2
               QQ_PACK(IFACT)=PPN(L,N)
               IFACT = IFACT+1
            ENDDO
         ENDDO
      ELSE
         ALLOCATE(QQ_PACK(NCELLS))
         DO L = L1,L2
            DO N = N1,N2
               QQ_PACK(IFACT)=PP(L,N)
               IFACT = IFACT+1
            ENDDO
         ENDDO
      ENDIF

   END SUBROUTINE PACK_VTK_BNDF

   !> \brief Gather the boundary values on a mesh's triangulated geometry

   !>

   !> \param NM Mesh number

   !> \param CELL_CENTERED True if the quantity is written on cells rather than points

   !> \param BF Boundary file the quantity belongs to

   !> \param VALS Gathered values, allocated here (out)


   SUBROUTINE PACK_VTK_GEOM(NM,CELL_CENTERED,BF,VALS)
      INTEGER, INTENT(IN) :: NM
      LOGICAL, INTENT(IN) :: CELL_CENTERED
      TYPE(BOUNDARY_FILE_TYPE), POINTER, INTENT(IN) :: BF
      REAL(FB), ALLOCATABLE, DIMENSION(:) :: VALS, VERTS, VERT_VALS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: FACES, VERT_UNIQUE
      INTEGER :: IND
      INTEGER :: NFACES, NVERTS, NVALS

      IND  = ABS(BF%INDEX)
      CALL POINT_TO_MESH(NM)
      CALL GET_GEOMSIZES('INBOUND_FACES',0,0,0,0,0,0,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS)

      IF (CELL_CENTERED) THEN
         NVALS = NFACES
         ALLOCATE(VALS(NFACES))
         ! get values at geometry faces
         CALL GET_GEOMVALS(CELL_CENTERED,CELL_CENTERED,'INBOUND_FACES',&
                          I1,I2,J1,J2,K1,K2,NFACES,NFACES_CUTCELLS,VALS,&
                          IND,BF%Y_INDEX,BF%Z_INDEX,BF%PART_INDEX,0,0,BF%PROP_INDEX,0,T,DT,NM,NF)
      ELSE
         NVALS = NVERTS
         ALLOCATE(VALS(MAX(NVERTS,NFACES)))

         ! get values at geometry nodes
         ALLOCATE(VERTS(3*NVERTS))
         ALLOCATE(FACES(3*NFACES))
         ALLOCATE(LOCATIONS(NFACES))
         ALLOCATE(VERT_UNIQUE(NVERTS))
         ALLOCATE(VERT_VALS(NVERTS))

         CALL GET_GEOMVALS(CELL_CENTERED,CELL_CENTERED,'INBOUND_FACES',&
                           I1,I2,J1,J2,K1,K2,NFACES,NFACES_CUTCELLS,VALS,&
                           IND,BF%Y_INDEX,BF%Z_INDEX,BF%PART_INDEX,0,0,BF%PROP_INDEX,0,T,DT,NM,NF)

         ! these two routines need to be moved and called only once
         CALL GET_GEOMINFO('INBOUND_FACES',I1,I2,J1,J2,K1,K2,NVERTS,NVERTS_CUTCELLS,NFACES,NFACES_CUTCELLS,VERTS,FACES,LOCATIONS)
         CALL MAKE_UNIQUE_VERT_ARRAY(VERTS, VERT_UNIQUE, NVERTS)

         CALL AVERAGE_FACE_VALUES(VERT_UNIQUE, VERT_VALS, NVERTS, FACES, VALS, NFACES)
         VALS(1:NVERTS) = VERT_VALS(1:NVERTS)

         DEALLOCATE(VERTS)
         DEALLOCATE(FACES)
         DEALLOCATE(LOCATIONS)
         DEALLOCATE(VERT_UNIQUE)
         DEALLOCATE(VERT_VALS)
      ENDIF
   END SUBROUTINE PACK_VTK_GEOM
END SUBROUTINE DUMP_BNDF_VTKHDF

#endif


!------------------------------------------------------------------------------
! MPI exchanges needed to replicate the per-mesh sizing information that the VTKHDF
! writers need on every rank.  Previously carried in main.f90.
!------------------------------------------------------------------------------

#ifdef WITH_HDF5
!> \brief Exchange information needed for output VTK geometry data

SUBROUTINE EXCHANGE_NOBST_INFO

USE GLOBAL_CONSTANTS, ONLY: N_MPI_PROCESSES,MY_RANK,PROCESS,NMESHES
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MPI_F08
INTEGER :: NM,IERR
INTEGER, ALLOCATABLE, DIMENSION(:) :: INT_EXC
REAL(EB) :: TNOW

IF (N_MPI_PROCESSES==1) THEN
   DO NM=1,NMESHES
      MESHES(NM)%NP = NP(NM)
      MESHES(NM)%NC = NC(NM)
   ENDDO
   RETURN
ENDIF

TNOW = CURRENT_TIME()

! Each rank knows N_OBST, NP and NC for the meshes it owns and contributes zero for the
! rest, so a single MAX reduction leaves every rank with the full picture.  This used to
! be a rank-0 gather followed by a rank-0 scatter, both written as O(NMESHES*N_MPI_PROCESSES)
! blocking point-to-point exchanges.

ALLOCATE(INT_EXC(3*NMESHES)) ; INT_EXC = 0

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   INT_EXC(3*NM-2) = MESHES(NM)%N_OBST
   INT_EXC(3*NM-1) = NP(NM)
   INT_EXC(3*NM  ) = NC(NM)
ENDDO

CALL MPI_ALLREDUCE(MPI_IN_PLACE,INT_EXC,3*NMESHES,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)

DO NM=1,NMESHES
   MESHES(NM)%N_OBST = MAX(MESHES(NM)%N_OBST,INT_EXC(3*NM-2))
   MESHES(NM)%NP     = INT_EXC(3*NM-1)
   MESHES(NM)%NC     = INT_EXC(3*NM  )
ENDDO

DEALLOCATE(INT_EXC)

T_USED(11) = T_USED(11) + CURRENT_TIME() - TNOW

END SUBROUTINE EXCHANGE_NOBST_INFO
#endif

!> \brief Number of Smoke3D points in a mesh

!>

!> \param NM Mesh number


FUNCTION NP(NM)
   INTEGER, INTENT(IN) :: NM
   INTEGER :: NX, NY, NZ, NP
   NX = SIZE(MESHES(NM)%X)
   NY = SIZE(MESHES(NM)%Y)
   NZ = SIZE(MESHES(NM)%Z)
   NP = NX*NY*NZ
END FUNCTION NP

!> \brief Number of Smoke3D cells in a mesh

!>

!> \param NM Mesh number


FUNCTION NC(NM)
   INTEGER, INTENT(IN) :: NM
   INTEGER :: NX, NY, NZ, NC
   NX = SIZE(MESHES(NM)%X)
   NY = SIZE(MESHES(NM)%Y)
   NZ = SIZE(MESHES(NM)%Z)
   NC = (NX-1)*(NY-1)*(NZ-1)
END FUNCTION NC


!> \brief Exchange information needed for output VTK slice data

SUBROUTINE EXCHANGE_NSLICE_INFO

USE GLOBAL_CONSTANTS, ONLY: N_MPI_PROCESSES,MY_RANK,PROCESS,NMESHES
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MPI_F08
TYPE (MPI_STATUS) :: MPISTATUS
INTEGER :: NM,IERR,IPROC,I
INTEGER, PARAMETER :: PADDING = 256
LOGICAL, ALLOCATABLE, DIMENSION(:) :: EMPTY_UNIQUE_SLICE, UNIQUE_SLICE_IS_SL3D
CHARACTER(200), ALLOCATABLE, DIMENSION(:) :: ALL_SLICE_QUANTITIES
CHARACTER(:), ALLOCATABLE :: STRING_BUF
REAL(EB) :: TNOW

IF (N_MPI_PROCESSES==1) RETURN

TNOW = CURRENT_TIME()

! Copy SLICE components to mesh 1 for VTK output
IF (MY_RANK==0) ALLOCATE(EMPTY_UNIQUE_SLICE(MESHES(1)%N_UNIQUE_SLCF))
ALLOCATE(ALL_SLICE_QUANTITIES(MESHES(1)%N_SLCF_O*3))
! Whether or not a slice is contained within the mesh
IF (MESHES(1)%N_UNIQUE_SLCF>0) THEN
   ! Get empty slice information from other meshes to rank 0
   DO NM=1,NMESHES
      DO IPROC=1,N_MPI_PROCESSES-1
         IF (MY_RANK==IPROC.AND.IPROC==PROCESS(NM)) THEN
            IF (MY_RANK/=0) THEN
               CALL MPI_SEND(MESHES(NM)%EMPTY_UNIQUE_SLICE,MESHES(1)%N_UNIQUE_SLCF,MPI_LOGICAL,0,PROCESS(NM),MPI_COMM_WORLD,IERR)
            ENDIF
         ELSEIF (MY_RANK==0) THEN
            IF (PROCESS(NM)==0) CYCLE
            IF (PROCESS(NM)/=IPROC) CYCLE
            CALL MPI_RECV(EMPTY_UNIQUE_SLICE(1),MESHES(1)%N_UNIQUE_SLCF,MPI_LOGICAL,IPROC,IPROC,MPI_COMM_WORLD,MPISTATUS,IERR)
            MESHES(NM)%EMPTY_UNIQUE_SLICE = EMPTY_UNIQUE_SLICE
         ENDIF
      ENDDO
   ENDDO
   ! Send slice topology info from rank 0 to other ranks
   IF (.NOT.ALLOCATED(MESHES(1)%UNIQUE_SLICE_IS_SL3D)) THEN
      ALLOCATE(MESHES(1)%UNIQUE_SLICE_IS_SL3D(MESHES(1)%N_UNIQUE_SLCF))
   ENDIF
   ALLOCATE(UNIQUE_SLICE_IS_SL3D(MESHES(1)%N_UNIQUE_SLCF))
   IF (MY_RANK==0) UNIQUE_SLICE_IS_SL3D = MESHES(1)%UNIQUE_SLICE_IS_SL3D
   CALL MPI_BCAST(UNIQUE_SLICE_IS_SL3D,MESHES(1)%N_UNIQUE_SLCF,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
   MESHES(1)%UNIQUE_SLICE_IS_SL3D = UNIQUE_SLICE_IS_SL3D
   DEALLOCATE(UNIQUE_SLICE_IS_SL3D)
   ! Get slice quantities from other ranks
   DO NM=1,NMESHES
      DO IPROC=1,N_MPI_PROCESSES-1
         IF (MY_RANK==IPROC.AND.IPROC==PROCESS(NM)) THEN
            IF (MY_RANK/=0) THEN
               ALLOCATE(CHARACTER(LEN=MESHES(1)%N_SLCF_O*(200+PADDING)*3) :: STRING_BUF)
               DO I=1,MESHES(1)%N_SLCF_O*3
                  STRING_BUF((I-1)*(200+PADDING)+1:(I-1)*(200+PADDING)+200) = MESHES(1)%ALL_SLICE_QUANTITIES(I)
               ENDDO
               CALL MPI_SEND(STRING_BUF,MESHES(1)%N_SLCF_O*(200+PADDING)*3,MPI_CHARACTER,0,PROCESS(NM),MPI_COMM_WORLD,IERR)
               DEALLOCATE(STRING_BUF)
            ENDIF
         ELSEIF (MY_RANK==0) THEN
            IF (PROCESS(NM)==0) CYCLE
            IF (PROCESS(NM)/=IPROC) CYCLE
            ALLOCATE(CHARACTER(MESHES(1)%N_SLCF_O*(200+PADDING)*3) :: STRING_BUF)
            CALL MPI_RECV(STRING_BUF,MESHES(1)%N_SLCF_O*(200+PADDING)*3,MPI_CHARACTER,IPROC,IPROC,MPI_COMM_WORLD,MPISTATUS,IERR)
            QUANTITY_LOOP: DO I=1,MESHES(1)%N_SLCF_O*3
               ALL_SLICE_QUANTITIES(I) = STRING_BUF((I-1)*(200+PADDING)+1:(I-1)*(200+PADDING)+200)
               IF (ALL_SLICE_QUANTITIES(I)=="") CYCLE QUANTITY_LOOP
               MESHES(1)%ALL_SLICE_QUANTITIES(I) = ALL_SLICE_QUANTITIES(I)
            ENDDO QUANTITY_LOOP
            DEALLOCATE(STRING_BUF)
         ENDIF
      ENDDO
   ENDDO
   ! Send slice quantity info from rank 0 to other ranks
   ALLOCATE(CHARACTER(LEN=MESHES(1)%N_SLCF_O*(200+PADDING)*3) :: STRING_BUF)
   IF (MY_RANK==0) THEN
      DO I=1,MESHES(1)%N_SLCF_O*3
         STRING_BUF((I-1)*(200+PADDING)+1:(I-1)*(200+PADDING)+200) = MESHES(1)%ALL_SLICE_QUANTITIES(I)
      ENDDO
   ENDIF
   CALL MPI_BCAST(STRING_BUF,MESHES(1)%N_SLCF_O*(200+PADDING)*3,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
   IF (MY_RANK/=0) THEN
      QUANTITY_LOOP2: DO I=1,MESHES(1)%N_SLCF_O*3
         ALL_SLICE_QUANTITIES(I) = STRING_BUF((I-1)*(200+PADDING)+1:(I-1)*(200+PADDING)+200)
         MESHES(1)%ALL_SLICE_QUANTITIES(I) = ALL_SLICE_QUANTITIES(I)
      ENDDO QUANTITY_LOOP2
   ENDIF
   DEALLOCATE(STRING_BUF)
ENDIF

IF (MY_RANK==0) DEALLOCATE(EMPTY_UNIQUE_SLICE)
DEALLOCATE(ALL_SLICE_QUANTITIES)

T_USED(11) = T_USED(11) + CURRENT_TIME() - TNOW

END SUBROUTINE EXCHANGE_NSLICE_INFO

!> \brief Exchange information needed for output VTK boundary data

SUBROUTINE EXCHANGE_NPATCH_INFO

USE GLOBAL_CONSTANTS, ONLY: N_MPI_PROCESSES,MY_RANK,PROCESS,NMESHES
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MPI_F08
INTEGER :: NM,IERR
INTEGER, ALLOCATABLE, DIMENSION(:) :: N_PATCH_ALL
REAL(EB) :: TNOW

IF (N_MPI_PROCESSES==1) RETURN

TNOW = CURRENT_TIME()

! One MAX reduction in place of the previous O(NMESHES*N_MPI_PROCESSES) rank-0 gather.
! Every rank now ends up with N_PATCH for every mesh, not just rank 0.

ALLOCATE(N_PATCH_ALL(NMESHES)) ; N_PATCH_ALL = 0

DO NM=1,NMESHES
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   N_PATCH_ALL(NM) = MESHES(NM)%N_PATCH
ENDDO

CALL MPI_ALLREDUCE(MPI_IN_PLACE,N_PATCH_ALL,NMESHES,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)

DO NM=1,NMESHES
   MESHES(NM)%N_PATCH = N_PATCH_ALL(NM)
ENDDO

DEALLOCATE(N_PATCH_ALL)

T_USED(11) = T_USED(11) + CURRENT_TIME() - TNOW

END SUBROUTINE EXCHANGE_NPATCH_INFO


#ifdef WITH_HDF5

!------------------------------------------------------------------------------
! Lagrangian particle output
!
! Unlike the slice, boundary and Smoke3D files, the particle geometry changes at every
! output time: particles are inserted and removed, and the ones that survive move.  The
! file is therefore a temporal VTKHDF UnstructuredGrid whose geometry arrays (Points,
! Connectivity, Offsets, Types and the per-part counts) are appended at every step
! alongside the point data.  Steps/PointOffsets, Steps/CellOffsets,
! Steps/ConnectivityIdOffsets and Steps/PartOffsets record where each step begins.  Each
! particle is one VTK_VERTEX cell, and connectivity is part-local, as it is for the
! other VTKHDF outputs.
!------------------------------------------------------------------------------

!> \brief Return the current length of the slowest-varying dimension of a dataset

SUBROUTINE GET_DATASET_LENGTH(GROUP_ID,SNAME,LENGTH)
INTEGER(HID_T), INTENT(IN) :: GROUP_ID
CHARACTER(*), INTENT(IN) :: SNAME
INTEGER(HSIZE_T), INTENT(OUT) :: LENGTH
INTEGER(HID_T) :: DSET_ID,DATASPACE
INTEGER(HSIZE_T), DIMENSION(2) :: CDIM,MDIM
INTEGER :: ERROR,RANK
LOGICAL :: LINK_EXISTS

LENGTH = 0_HSIZE_T
IF (.NOT.VTK_DSET_CACHE_FIND(GROUP_ID,SNAME,DSET_ID)) THEN
   CALL H5LEXISTS_F(GROUP_ID,SNAME,LINK_EXISTS,ERROR)
   IF (.NOT.LINK_EXISTS) RETURN
   CALL H5DOPEN_F(GROUP_ID,SNAME,DSET_ID,ERROR)
   CALL VTK_DSET_CACHE_ADD(GROUP_ID,SNAME,DSET_ID)
ENDIF
CALL H5DGET_SPACE_F(DSET_ID,DATASPACE,ERROR)
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(DATASPACE,RANK,ERROR)
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE,CDIM,MDIM,ERROR)
LENGTH = CDIM(RANK)
CALL H5SCLOSE_F(DATASPACE,ERROR)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

END SUBROUTINE GET_DATASET_LENGTH


!> \brief Append one integer to a rank-1 dataset under the Steps group

SUBROUTINE APPEND_STEP_INDEX(GROUP_ID,SNAME,IVAL,CRP_LIST,PLIST_ID)
INTEGER(HID_T), INTENT(IN) :: GROUP_ID
INTEGER(HID_T), INTENT(INOUT) :: CRP_LIST,PLIST_ID
CHARACTER(*), INTENT(IN) :: SNAME
INTEGER(IB32), INTENT(IN) :: IVAL
INTEGER(HID_T) :: DSET_ID
INTEGER(HSIZE_T), DIMENSION(1) :: CDIM1,DDIM1,MDIM
INTEGER(IB32), DIMENSION(1) :: INDATA
INTEGER :: ERROR

CDIM1 = (/1_HSIZE_T/) ; DDIM1 = (/0_HSIZE_T/) ; MDIM = (/H5S_UNLIMITED_F/)
INDATA = IVAL
CALL PARALLEL_INIT_I32(GROUP_ID,SNAME,CRP_LIST,1,CDIM1,DDIM1,DSET_ID,PLIST_ID,MDIM)
CALL APPEND_RANK1_DATASET_I32(DSET_ID,PLIST_ID,INDATA,1)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

END SUBROUTINE APPEND_STEP_INDEX


!> \brief Extend a rank-1 dataset by N_TOTAL and write this rank's blocks in one call
!>
!> \param GROUP_ID Group holding the dataset
!> \param SNAME Dataset name
!> \param CHUNK Chunk size used if the dataset has to be created
!> \param N_TOTAL Number of elements added by all ranks together
!> \param N_BLK Number of blocks this rank contributes
!> \param BLK_START Start of each block, relative to the previous end of the dataset
!> \param BLK_COUNT Length of each block
!> \param IDATA Integer data, blocks concatenated (pass a zero-size array if N_BLK is 0)
!> \param FDATA Real data, as an alternative to IDATA
!> \param BDATA One-byte integer data, as an alternative to IDATA

SUBROUTINE PART_APPEND_1D(GROUP_ID,SNAME,CHUNK,N_TOTAL,N_BLK,BLK_START,BLK_COUNT,CRP_LIST,PLIST_ID,&
                          IDATA,FDATA,BDATA)
INTEGER(HID_T), INTENT(IN) :: GROUP_ID
INTEGER(HID_T), INTENT(INOUT) :: CRP_LIST,PLIST_ID
CHARACTER(*), INTENT(IN) :: SNAME
INTEGER, INTENT(IN) :: CHUNK,N_TOTAL,N_BLK
INTEGER, DIMENSION(:), INTENT(IN) :: BLK_START,BLK_COUNT
INTEGER, DIMENSION(*), INTENT(IN), OPTIONAL :: IDATA
REAL(FB), DIMENSION(*), INTENT(IN), OPTIONAL :: FDATA
INTEGER(IB8), DIMENSION(*), INTENT(IN), OPTIONAL :: BDATA
INTEGER(HID_T) :: DSET_ID,MEMSPACE,DATASPACE
INTEGER(HSIZE_T), DIMENSION(1) :: CDIM,MDIM,SIZE1,START1,EXTENT1
INTEGER(HSIZE_T) :: BASE,NLOC
INTEGER :: I,ERROR

CDIM = (/INT(MAX(1,CHUNK),HSIZE_T)/)
MDIM = (/H5S_UNLIMITED_F/)
START1 = (/0_HSIZE_T/)
IF (PRESENT(BDATA)) THEN
   CALL PARALLEL_INIT_U8(GROUP_ID,SNAME,CRP_LIST,1,CDIM,START1,DSET_ID,PLIST_ID,MDIM)
ELSEIF (PRESENT(FDATA)) THEN
   CALL PARALLEL_INIT_F32(GROUP_ID,SNAME,CRP_LIST,1,CDIM,START1,DSET_ID,PLIST_ID,MDIM)
ELSE
   CALL PARALLEL_INIT_I32(GROUP_ID,SNAME,CRP_LIST,1,CDIM,START1,DSET_ID,PLIST_ID,MDIM)
ENDIF

CALL H5DGET_SPACE_F(DSET_ID,DATASPACE,ERROR)
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE,CDIM,MDIM,ERROR)
CALL H5SCLOSE_F(DATASPACE,ERROR)
BASE = CDIM(1)
SIZE1 = BASE + INT(N_TOTAL,HSIZE_T)
CALL H5DSET_EXTENT_F(DSET_ID,SIZE1,ERROR)

CALL H5DGET_SPACE_F(DSET_ID,DATASPACE,ERROR)
CALL H5SSELECT_NONE_F(DATASPACE,ERROR)
NLOC = 0_HSIZE_T
DO I=1,N_BLK
   IF (BLK_COUNT(I)<=0) CYCLE
   START1  = BASE + INT(BLK_START(I),HSIZE_T)
   EXTENT1 = INT(BLK_COUNT(I),HSIZE_T)
   CALL H5SSELECT_HYPERSLAB_F(DATASPACE,H5S_SELECT_OR_F,START1,EXTENT1,ERROR)
   NLOC = NLOC + EXTENT1(1)
ENDDO

EXTENT1 = NLOC
CALL H5SCREATE_SIMPLE_F(1,EXTENT1,MEMSPACE,ERROR)
IF (NLOC==0_HSIZE_T) CALL H5SSELECT_NONE_F(MEMSPACE,ERROR)
IF (PRESENT(BDATA)) THEN
   CALL H5DWRITE_F(DSET_ID,H5T_STD_U8LE,BDATA,EXTENT1,ERROR,&
                   MEM_SPACE_ID=MEMSPACE,FILE_SPACE_ID=DATASPACE,XFER_PRP=PLIST_ID)
ELSEIF (PRESENT(FDATA)) THEN
   CALL H5DWRITE_F(DSET_ID,H5T_IEEE_F32LE,FDATA,EXTENT1,ERROR,&
                   MEM_SPACE_ID=MEMSPACE,FILE_SPACE_ID=DATASPACE,XFER_PRP=PLIST_ID)
ELSE
   CALL H5DWRITE_F(DSET_ID,H5T_STD_I32LE,IDATA,EXTENT1,ERROR,&
                   MEM_SPACE_ID=MEMSPACE,FILE_SPACE_ID=DATASPACE,XFER_PRP=PLIST_ID)
ENDIF

CALL H5SCLOSE_F(MEMSPACE,ERROR)
CALL H5SCLOSE_F(DATASPACE,ERROR)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

END SUBROUTINE PART_APPEND_1D


!> \brief Extend a (3,N) dataset by N_TOTAL points and write this rank's blocks in one call

SUBROUTINE PART_APPEND_2D(GROUP_ID,SNAME,CHUNK,N_TOTAL,N_BLK,BLK_START,BLK_COUNT,CRP_LIST,PLIST_ID,FDATA,NOFILTER)
INTEGER(HID_T), INTENT(IN) :: GROUP_ID
INTEGER(HID_T), INTENT(INOUT) :: CRP_LIST,PLIST_ID
CHARACTER(*), INTENT(IN) :: SNAME
INTEGER, INTENT(IN) :: CHUNK,N_TOTAL,N_BLK
INTEGER, DIMENSION(:), INTENT(IN) :: BLK_START,BLK_COUNT
REAL(FB), DIMENSION(*), INTENT(IN) :: FDATA
LOGICAL, INTENT(IN), OPTIONAL :: NOFILTER
INTEGER(HID_T) :: DSET_ID,MEMSPACE,DATASPACE
INTEGER(HSIZE_T), DIMENSION(2) :: CDIM,MDIM,SIZE2,START2,EXTENT2
INTEGER(HSIZE_T) :: BASE,NLOC
INTEGER :: I,ERROR

CDIM = (/3_HSIZE_T,INT(MAX(1,CHUNK),HSIZE_T)/)
MDIM = (/3_HSIZE_T,H5S_UNLIMITED_F/)
START2 = (/3_HSIZE_T,0_HSIZE_T/)
CALL PARALLEL_INIT_F32(GROUP_ID,SNAME,CRP_LIST,2,CDIM,START2,DSET_ID,PLIST_ID,MDIM,NOFILTER)

CALL H5DGET_SPACE_F(DSET_ID,DATASPACE,ERROR)
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DATASPACE,CDIM,MDIM,ERROR)
CALL H5SCLOSE_F(DATASPACE,ERROR)
BASE = CDIM(2)
SIZE2 = (/3_HSIZE_T,BASE+INT(N_TOTAL,HSIZE_T)/)
CALL H5DSET_EXTENT_F(DSET_ID,SIZE2,ERROR)

CALL H5DGET_SPACE_F(DSET_ID,DATASPACE,ERROR)
CALL H5SSELECT_NONE_F(DATASPACE,ERROR)
NLOC = 0_HSIZE_T
DO I=1,N_BLK
   IF (BLK_COUNT(I)<=0) CYCLE
   START2  = (/0_HSIZE_T,BASE+INT(BLK_START(I),HSIZE_T)/)
   EXTENT2 = (/3_HSIZE_T,INT(BLK_COUNT(I),HSIZE_T)/)
   CALL H5SSELECT_HYPERSLAB_F(DATASPACE,H5S_SELECT_OR_F,START2,EXTENT2,ERROR)
   NLOC = NLOC + EXTENT2(2)
ENDDO

EXTENT2 = (/3_HSIZE_T,NLOC/)
CALL H5SCREATE_SIMPLE_F(2,EXTENT2,MEMSPACE,ERROR)
IF (NLOC==0_HSIZE_T) CALL H5SSELECT_NONE_F(MEMSPACE,ERROR)
CALL H5DWRITE_F(DSET_ID,H5T_IEEE_F32LE,FDATA,EXTENT2,ERROR,&
                MEM_SPACE_ID=MEMSPACE,FILE_SPACE_ID=DATASPACE,XFER_PRP=PLIST_ID)

CALL H5SCLOSE_F(MEMSPACE,ERROR)
CALL H5SCLOSE_F(DATASPACE,ERROR)
CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

END SUBROUTINE PART_APPEND_2D


!> \brief Create one VTKHDF particle file per Lagrangian particle class

!> \brief Record a particle class's colour on the file's VTKHDF group
!>
!> \param N Particle class index
!>
!> Three floats, written once when the file is created, so the file still says what
!> colour the class is without the generated ParaView script alongside it.
!>
!> An HDF5 attribute rather than a FieldData array: FieldData is temporal in a VTKHDF
!> time series, so an array there needs a Steps/FieldDataOffsets entry for every step,
!> and without one the reader fails the entire read request rather than skipping it.

SUBROUTINE WRITE_PART_COLOR(N)
INTEGER, INTENT(IN) :: N
INTEGER(HSIZE_T), DIMENSION(1) :: ADIMS
REAL(FB), DIMENSION(3) :: RGB
RGB = GET_PART_CLASS_COLOR(N)
ADIMS = (/3_HSIZE_T/)
CALL ADD_ATTRIBUTE_F32(HDF_PART_G1(N),ADIMS,1,'COLOR',RGB)
END SUBROUTINE WRITE_PART_COLOR


!> \brief Create one VTKHDF file per particle class


!>


!> Unlike the other outputs the particle geometry changes at every output time, so


!> only the file and its groups are made here; the points arrive with the data.



SUBROUTINE INITIALIZE_VTKHDF_PART()

CHARACTER(FN_LENGTH) :: FILENAME
INTEGER :: N

IF (N_LAGRANGIAN_CLASSES<1) RETURN

ALLOCATE(HDF_PART_FILE_ID(N_LAGRANGIAN_CLASSES))
ALLOCATE(HDF_PART_PLIST_ID(N_LAGRANGIAN_CLASSES))
ALLOCATE(HDF_PART_CRP_LIST(N_LAGRANGIAN_CLASSES))
ALLOCATE(HDF_PART_G1(N_LAGRANGIAN_CLASSES)) ; ALLOCATE(HDF_PART_G2(N_LAGRANGIAN_CLASSES))
ALLOCATE(HDF_PART_G3(N_LAGRANGIAN_CLASSES)) ; ALLOCATE(HDF_PART_G4(N_LAGRANGIAN_CLASSES))
ALLOCATE(HDF_PART_G5(N_LAGRANGIAN_CLASSES)) ; ALLOCATE(HDF_PART_G6(N_LAGRANGIAN_CLASSES))
ALLOCATE(HDF_PART_G7(N_LAGRANGIAN_CLASSES))
ALLOCATE(HDF_PART_COUNTER(N_LAGRANGIAN_CLASSES))
HDF_PART_COUNTER = 0

DO N=1,N_LAGRANGIAN_CLASSES
   LPC => LAGRANGIAN_PARTICLE_CLASS(N)
   WRITE(FILENAME,'(A,A,A,A)') TRIM(VTK_DIR)//TRIM(CHID),'_PART_',TRIM(LPC%ID),'.vtkhdf'
   CALL CREATE_OPEN_VTKHDF_SERIES(FILENAME,HDF_PART_FILE_ID(N),HDF_PART_PLIST_ID(N),&
      HDF_PART_G1(N),HDF_PART_G2(N),HDF_PART_G3(N),HDF_PART_G4(N),&
      HDF_PART_G5(N),HDF_PART_G6(N),HDF_PART_G7(N))
   CALL WRITE_PART_COLOR(N)
ENDDO

IF (.NOT.VTK_KEEPOPEN) CALL CLOSE_VTKHDF_PART()

END SUBROUTINE INITIALIZE_VTKHDF_PART


!> \brief Reopen the particle VTKHDF files, when the files are not held open



SUBROUTINE OPEN_VTKHDF_PART()
CHARACTER(FN_LENGTH) :: FILENAME
INTEGER :: N
DO N=1,N_LAGRANGIAN_CLASSES
   LPC => LAGRANGIAN_PARTICLE_CLASS(N)
   WRITE(FILENAME,'(A,A,A,A)') TRIM(VTK_DIR)//TRIM(CHID),'_PART_',TRIM(LPC%ID),'.vtkhdf'
   CALL OPEN_VTKHDF_SERIES(FILENAME,HDF_PART_FILE_ID(N),HDF_PART_PLIST_ID(N),&
      HDF_PART_G1(N),HDF_PART_G2(N),HDF_PART_G3(N),HDF_PART_G4(N),&
      HDF_PART_G5(N),HDF_PART_G6(N),HDF_PART_G7(N))
ENDDO
END SUBROUTINE OPEN_VTKHDF_PART


!> \brief Close the particle VTKHDF files



SUBROUTINE CLOSE_VTKHDF_PART()
INTEGER :: N
IF (.NOT.ALLOCATED(HDF_PART_FILE_ID)) RETURN
DO N=1,N_LAGRANGIAN_CLASSES
   CALL CLOSE_VTKHDF_SERIES(HDF_PART_FILE_ID(N),&
      HDF_PART_G1(N),HDF_PART_G2(N),HDF_PART_G3(N),HDF_PART_G4(N),&
      HDF_PART_G5(N),HDF_PART_G6(N),HDF_PART_G7(N))
ENDDO
END SUBROUTINE CLOSE_VTKHDF_PART


!> \brief Append one time step of Lagrangian particle data to the VTKHDF particle files
!> \param T Current simulation time (s)

SUBROUTINE DUMP_PART_VTKHDF(T)

USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), INTENT(IN) :: T
INTEGER :: N,NN,NM,IP,NPP,IZERO,IERR,IOFF,NQ,NPTS_STEP,NOFF_STEP,N_LOCAL,N_MESH_LOCAL,NLOC_PTS,NP_LOCAL,NB
INTEGER(HSIZE_T) :: POINT_BASE,PART_BASE
INTEGER, ALLOCATABLE, DIMENSION(:) :: NPOINTS,POINT_OFFSET,OFFS_OFFSET,TA,CONN,OFFS,&
                                      BLK_START,BLK_COUNT,OBLK_START,OBLK_COUNT
REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: VERTICES
REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: QP
INTEGER(IB8), ALLOCATABLE, DIMENSION(:) :: TYPES
INTEGER(HID_T) :: DSET_ID
INTEGER(HSIZE_T), DIMENSION(1) :: CDIM1,DDIM1,MDIM
REAL(FB), DIMENSION(1) :: VTK_T
INTEGER :: ERROR
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC

IF (N_LAGRANGIAN_CLASSES<1) RETURN

ALLOCATE(NPOINTS(NMESHES)) ; ALLOCATE(POINT_OFFSET(NMESHES)) ; ALLOCATE(OFFS_OFFSET(NMESHES))
N_MESH_LOCAL = MAX(1,UPPER_MESH_INDEX-LOWER_MESH_INDEX+1)
ALLOCATE(BLK_START(N_MESH_LOCAL))  ; ALLOCATE(BLK_COUNT(N_MESH_LOCAL))
ALLOCATE(OBLK_START(N_MESH_LOCAL)) ; ALLOCATE(OBLK_COUNT(N_MESH_LOCAL))

CDIM1 = (/1_HSIZE_T/) ; DDIM1 = (/0_HSIZE_T/) ; MDIM = (/H5S_UNLIMITED_F/)

CLASS_LOOP: DO N=1,N_LAGRANGIAN_CLASSES

   LPC => LAGRANGIAN_PARTICLE_CLASS(N)
   NQ = LPC%N_QUANTITIES

   ! Count this rank's particles of class N in each mesh it owns, then make the counts
   ! known to every rank with a single reduction.

   NPOINTS = 0
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      NPP = 0
      DO IP=1,NLP
         LP=>LAGRANGIAN_PARTICLE(IP)
         IF (LP%SHOW .AND. LP%CLASS_INDEX==N) NPP = NPP + 1
      ENDDO
      NPOINTS(NM) = NPP
   ENDDO
   IF (N_MPI_PROCESSES>1) &
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,NPOINTS,NMESHES,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)

   NPTS_STEP = 0
   NOFF_STEP = 0
   DO NM=1,NMESHES
      POINT_OFFSET(NM) = NPTS_STEP
      OFFS_OFFSET(NM)  = NOFF_STEP
      NPTS_STEP = NPTS_STEP + NPOINTS(NM)
      NOFF_STEP = NOFF_STEP + NPOINTS(NM) + 1   ! NCELLS+1 offsets per part
   ENDDO

   ! Where this step starts within the appended arrays

   CALL GET_DATASET_LENGTH(HDF_PART_G1(N),'NumberOfPoints',PART_BASE)
   CALL GET_DATASET_LENGTH(HDF_PART_G1(N),'Types',POINT_BASE)

   ! Steps metadata.  One cell per particle, so the cell and connectivity-id offsets are
   ! the same as the point offset.

   VTK_T = REAL(T,FB)
   CALL PARALLEL_INIT_F32(HDF_PART_G5(N),'Values',HDF_PART_CRP_LIST(N),1,CDIM1,DDIM1,&
                          DSET_ID,HDF_PART_PLIST_ID(N),MDIM)
   CALL APPEND_RANK1_DATASET_F32(DSET_ID,HDF_PART_PLIST_ID(N),VTK_T,1)
   CALL VTK_DSET_RELEASE(DSET_ID,ERROR)

   CALL APPEND_STEP_INDEX(HDF_PART_G5(N),'PartOffsets',INT(PART_BASE,IB32),&
                          HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N))
   CALL APPEND_STEP_INDEX(HDF_PART_G5(N),'NumberOfParts',INT(NMESHES,IB32),&
                          HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N))
   CALL APPEND_STEP_INDEX(HDF_PART_G5(N),'PointOffsets',INT(POINT_BASE,IB32),&
                          HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N))
   CALL APPEND_STEP_INDEX(HDF_PART_G5(N),'CellOffsets',INT(POINT_BASE,IB32),&
                          HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N))
   CALL APPEND_STEP_INDEX(HDF_PART_G5(N),'ConnectivityIdOffsets',INT(POINT_BASE,IB32),&
                          HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N))

   HDF_PART_COUNTER(N) = HDF_PART_COUNTER(N) + 1
   CALL ADD_ATTRIBUTE_INT(HDF_PART_G5(N),CDIM1,'NSteps',HDF_PART_COUNTER(N))

   ! Per-part counts.  Every rank knows them for every mesh, so rank 0 writes the whole
   ! block and the others contribute nothing; the call itself stays collective.

   NB = 0
   IF (MY_RANK==0) NB = 1
   BLK_START(1) = 0 ; BLK_COUNT(1) = NMESHES
   CALL PART_APPEND_1D(HDF_PART_G1(N),'NumberOfPoints',NMESHES,NMESHES,NB,BLK_START,BLK_COUNT,&
                       HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),IDATA=NPOINTS)
   CALL PART_APPEND_1D(HDF_PART_G1(N),'NumberOfCells',NMESHES,NMESHES,NB,BLK_START,BLK_COUNT,&
                       HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),IDATA=NPOINTS)
   CALL PART_APPEND_1D(HDF_PART_G1(N),'NumberOfConnectivityIds',NMESHES,NMESHES,NB,BLK_START,BLK_COUNT,&
                       HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),IDATA=NPOINTS)

   ! Gather this rank's particles

   NLOC_PTS = SUM(NPOINTS(LOWER_MESH_INDEX:UPPER_MESH_INDEX))
   ALLOCATE(VERTICES(3,MAX(1,NLOC_PTS)))
   ALLOCATE(TA(MAX(1,NLOC_PTS)),STAT=IZERO)    ; CALL ChkMemErr('DUMP','TA',IZERO)
   ALLOCATE(CONN(MAX(1,NLOC_PTS)))
   ALLOCATE(TYPES(MAX(1,NLOC_PTS)))
   ALLOCATE(OFFS(MAX(1,NLOC_PTS+N_MESH_LOCAL)))
   ALLOCATE(QP(MAX(1,NLOC_PTS),MAX(1,NQ)),STAT=IZERO) ; CALL ChkMemErr('DUMP','QP',IZERO)

   N_LOCAL = 0
   NPP = 0
   IOFF = 0
   LOCAL_MESH_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      N_LOCAL = N_LOCAL + 1
      BLK_START(N_LOCAL)  = POINT_OFFSET(NM)
      BLK_COUNT(N_LOCAL)  = NPOINTS(NM)
      OBLK_START(N_LOCAL) = OFFS_OFFSET(NM)
      OBLK_COUNT(N_LOCAL) = NPOINTS(NM) + 1
      ! Connectivity and Offsets are part-local, as they are for the other VTKHDF outputs
      NP_LOCAL = 0
      IOFF = IOFF + 1
      OFFS(IOFF) = 0
      LOAD_LOOP: DO IP=1,NLP
         LP=>LAGRANGIAN_PARTICLE(IP)
         IF (.NOT.LP%SHOW .OR. LP%CLASS_INDEX/=N) CYCLE LOAD_LOOP
         IF (NP_LOCAL>=NPOINTS(NM)) EXIT LOAD_LOOP
         BC=>BOUNDARY_COORD(LP%BC_INDEX)
         NPP = NPP + 1
         NP_LOCAL = NP_LOCAL + 1
         TA(NPP) = LP%TAG
         TYPES(NPP) = 1_IB8            ! VTK_VERTEX
         CONN(NPP) = NP_LOCAL - 1
         IOFF = IOFF + 1
         OFFS(IOFF) = NP_LOCAL
         VERTICES(1:3,NPP) = REAL((/BC%X,BC%Y,BC%Z/),FB)
         DO NN=1,NQ
            QP(NPP,NN) = REAL(PARTICLE_OUTPUT(NM,T,LPC%QUANTITIES_INDEX(NN),IP,&
               Y_INDEX=LPC%QUANTITIES_Y_INDEX(NN),Z_INDEX=LPC%QUANTITIES_Z_INDEX(NN)),FB)
         ENDDO
      ENDDO LOAD_LOOP
   ENDDO LOCAL_MESH_LOOP

   ! Append the geometry and the point data, one collective write per array

   ! Particle coordinates are scattered, so deflate returns about 1.2x on them for the full
   ! CPU price.  The arrays around them return 4x to 130x and stay compressed.

   CALL PART_APPEND_2D(HDF_PART_G1(N),'Points',VTK_PART_CHUNK,NPTS_STEP,N_LOCAL,BLK_START,BLK_COUNT,&
                       HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),VERTICES,NOFILTER=.TRUE.)
   CALL PART_APPEND_1D(HDF_PART_G1(N),'Connectivity',VTK_PART_CHUNK,NPTS_STEP,N_LOCAL,BLK_START,BLK_COUNT,&
                       HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),IDATA=CONN)
   CALL PART_APPEND_1D(HDF_PART_G1(N),'Offsets',VTK_PART_CHUNK,NOFF_STEP,N_LOCAL,OBLK_START,OBLK_COUNT,&
                       HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),IDATA=OFFS)
   CALL PART_APPEND_1D(HDF_PART_G1(N),'Types',VTK_PART_CHUNK,NPTS_STEP,N_LOCAL,BLK_START,BLK_COUNT,&
                       HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),BDATA=TYPES)

   CALL PART_APPEND_1D(HDF_PART_G4(N),'TAG',VTK_PART_CHUNK,NPTS_STEP,N_LOCAL,BLK_START,BLK_COUNT,&
                       HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),IDATA=TA)
   CALL APPEND_STEP_INDEX(HDF_PART_G7(N),'TAG',INT(POINT_BASE,IB32),&
                          HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N))
   DO NN=1,NQ
      CALL PART_APPEND_1D(HDF_PART_G4(N),TRIM(LPC%SMOKEVIEW_LABEL(NN)),VTK_PART_CHUNK,NPTS_STEP,&
                          N_LOCAL,BLK_START,BLK_COUNT,HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N),&
                          FDATA=QP(:,NN))
      CALL APPEND_STEP_INDEX(HDF_PART_G7(N),TRIM(LPC%SMOKEVIEW_LABEL(NN)),INT(POINT_BASE,IB32),&
                             HDF_PART_CRP_LIST(N),HDF_PART_PLIST_ID(N))
   ENDDO

   DEALLOCATE(VERTICES) ; DEALLOCATE(TA)
   DEALLOCATE(CONN) ; DEALLOCATE(TYPES) ; DEALLOCATE(OFFS) ; DEALLOCATE(QP)

ENDDO CLASS_LOOP

DEALLOCATE(NPOINTS) ; DEALLOCATE(POINT_OFFSET) ; DEALLOCATE(OFFS_OFFSET)
DEALLOCATE(BLK_START) ; DEALLOCATE(BLK_COUNT) ; DEALLOCATE(OBLK_START) ; DEALLOCATE(OBLK_COUNT)

END SUBROUTINE DUMP_PART_VTKHDF

#endif



END MODULE VTK_FDS_INTERFACE



