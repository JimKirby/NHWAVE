        PROGRAM MAKE_NESTING_FILES
        USE GLOBAL

        USE INPUT_UTIL
        IMPLICIT NONE
        INTEGER :: I,J,K,line
        CHARACTER(LEN=80) File_Name
        INTEGER :: NumData,nwest,nsouth,neast,nnorth

        REAL,DIMENSION(:),ALLOCATABLE :: TIME_DATA
        REAL :: TRamp_DATA = 0.0
        REAL :: Ramp

! input data
         NumData = 1674   ! 5 hours 
! get data from input.txt

     ! read from input.txt
       FILE_NAME='input.txt'

     ! dimension                                             
       CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
       CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
       CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line)   

        N_COUPLING_EAST = 560
        N_COUPLING_WEST = 0
        N_COUPLING_SOUTH =0
        N_COUPLING_NORTH =672

        J_START_WEST = 1
        J_START_EAST = 1
        I_START_SOUTH = 1
        I_START_NORTH = 1

		OPEN(11,FILE='coupling.txt')
		OPEN(12,FILE='TIME_DATA.txt',status='old')
		OPEN(13,FILE='ETAE.txt',status='old')
		OPEN(14,FILE='UE.txt',status='old')
		OPEN(15,FILE='VE.txt',status='old')
		OPEN(16,FILE='ETAN.txt',status='old')
		OPEN(17,FILE='UN.txt',status='old')
		OPEN(18,FILE='VN.txt',status='old')
! Allocate

         ALLOCATE (TIME_DATA(NumData))
		 
         do K=1,NumData
           read(12,*) TIME_DATA(K)
		 ENDDO

         IF(N_COUPLING_EAST.GT.0)THEN
           ALLOCATE(U_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData),&
               V_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData),&
               W_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData),&
               Z_COUPLING_EAST(N_COUPLING_EAST,NumData), &
               P_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData) &
               )

           DO K=1,NumData
             read(13,*) (Z_COUPLING_EAST(I,K),I=1,N_COUPLING_EAST)
             read(14,*) ((U_COUPLING_EAST(I,J,K),I=1,N_COUPLING_EAST),J=1,Kglob)
             read(15,*) ((V_COUPLING_EAST(I,J,K),I=1,N_COUPLING_EAST),J=1,Kglob)
           ENDDO
		   W_COUPLING_EAST=0.0
           P_COUPLING_EAST=0.0
         ENDIF

         IF(N_COUPLING_WEST.GT.0)THEN
           ALLOCATE(U_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData),&
               V_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData),&
               W_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData),&
               Z_COUPLING_WEST(N_COUPLING_WEST,NumData), &
               P_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData) &
               )
           DO K=1,NumData
           DO I=1,N_COUPLING_WEST
             Z_COUPLING_WEST(I,K)=0.0
           DO J=1,Kglob
             U_COUPLING_WEST(I,J,K)=0.0
             V_COUPLING_WEST(I,J,K)=0.0
             W_COUPLING_WEST(I,J,K)=0.0
             P_COUPLING_WEST(I,J,K)=0.0
           ENDDO
           ENDDO
           ENDDO
         ENDIF

         IF(N_COUPLING_SOUTH.GT.0)THEN
           ALLOCATE(U_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData),&
               V_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData),&
               W_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData),&
               Z_COUPLING_SOUTH(N_COUPLING_SOUTH,NumData), &
               P_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData)&
               )
           DO K=1,NumData
           DO I=1,N_COUPLING_SOUTH
             Z_COUPLING_SOUTH(I,K)=0.0
           DO J=1,Kglob
             U_COUPLING_SOUTH(I,J,K)=0.0
             V_COUPLING_SOUTH(I,J,K)=0.0
             W_COUPLING_SOUTH(I,J,K)=0.0
             P_COUPLING_SOUTH(I,J,K)=0.0
           ENDDO
           ENDDO
           ENDDO
         ENDIF

         IF(N_COUPLING_NORTH.GT.0)THEN
           ALLOCATE(U_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData),&
               V_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData),&
               W_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData),&
               Z_COUPLING_NORTH(N_COUPLING_NORTH,NumData), &
               P_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData) &
               )
           DO K=1,NumData
             read(16,*) (Z_COUPLING_NORTH(I,K),I=1,N_COUPLING_NORTH)
             read(17,*) ((U_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             read(18,*) ((V_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
           ENDDO
		   W_COUPLING_NORTH=0.0
		   P_COUPLING_NORTH=0.0
         ENDIF


         WRITE(11,*)  'coupling data'
         WRITE(11,*)  'boundary info'
! boundary basic info including point number of coupling, start point, etc
! east
         WRITE(11,*)  'N_COUPLING_EAST,J_START_EAST'
         WRITE(11,'(2I8)') N_COUPLING_EAST,J_START_EAST
! west 
         WRITE(11,*)  'N_COUPLING_WEST,J_START_WEST'
         WRITE(11,'(2I8)') N_COUPLING_WEST,J_START_WEST
! south 
         WRITE(11,*)  'N_COUPLING_SOUTH,I_START_SOUTH'
         WRITE(11,'(2I8)') N_COUPLING_SOUTH,I_START_SOUTH
! north 
         WRITE(11,*)  'N_COUPLING_NORTH,I_START_NORTH'
         WRITE(11,'(2I8)') N_COUPLING_NORTH,I_START_NORTH

! read time and variable 
119      FORMAT(5E16.6)

       DO K=1,NumData

         Ramp=TANH(TIME_DATA(K)/TRAMP_DATA)

         WRITE(11,*) 'TIME_COUPLING'
         WRITE(11,'(F16.8)') TIME_DATA(K)
      
! east
         IF(N_COUPLING_EAST.GT.0)THEN
             WRITE(11,*)   'east'
             WRITE(11,119)(Z_COUPLING_EAST(I,K),I=1,N_COUPLING_EAST)
             WRITE(11,119)((U_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
             WRITE(11,119)((V_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
             WRITE(11,119)((W_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
             WRITE(11,119)((P_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
         ELSE
             WRITE(11,*) 'east'
         ENDIF

! west
         IF(N_COUPLING_WEST.GT.0)THEN
             WRITE(11,*)   'west'
             WRITE(11,119)(Z_COUPLING_WEST(I,K),I=1,N_COUPLING_WEST)
             WRITE(11,119)((U_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
             WRITE(11,119)((V_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
             WRITE(11,119)((W_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
             WRITE(11,119)((P_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
         ELSE
             WRITE(11,*) 'west'
         ENDIF          

! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
             WRITE(11,*)   'south'
             WRITE(11,119)(Z_COUPLING_SOUTH(I,K),I=1,N_COUPLING_SOUTH)
             WRITE(11,119)((U_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             WRITE(11,119)((V_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             WRITE(11,119)((W_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             WRITE(11,119)((P_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
         ELSE
             WRITE(11,*) 'south'
         ENDIF

! north
         IF(N_COUPLING_NORTH.GT.0)THEN
             WRITE(11,*)   'north'
             WRITE(11,119)(Z_COUPLING_NORTH(I,K),I=1,N_COUPLING_NORTH)
             WRITE(11,119)((U_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             WRITE(11,119)((V_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             WRITE(11,119)((W_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             WRITE(11,119)((P_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
         ELSE
             WRITE(11,*) 'north'
         ENDIF

        ENDDO ! end numdata
        close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)
		close(17)
		close(18)
        END
        








