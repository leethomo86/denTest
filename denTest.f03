      program denTest
!
!     This program tests the symmetry of a density matrix.
!
      use mqc_gaussian
!
!****x* Main/denTest
!*    NAME
!*      Density Test
!*
!*    SYNOPSIS
!*      Determine symmetry of a density matrix stored on a matrix file.
!
      implicit none
      character(len=:),allocatable::command,fileName,fileName2,help_path
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      integer::iOut=output_unit,iPrint=0,i,j,flag,nAlpha,nBeta,nBasis
      type(mqc_scf_integral)::density,overlap,mo,density2,overlap2,mo2
      logical::foundMO,foundDen
      real(kind=real64),parameter::thresh=1.0e-14,one=1.0d0
      type(mqc_scalar)::mo_overlap,density_dist
!
!*    USAGE
!*      denTest [-f <Matrix_file>] [--distance <matrix_file>] [--print-level <print_level>] 
!*              [--help]
!*
!*    OPTIONS
!
      j = 1
      do i=1,command_argument_count()
        if(i.ne.j) cycle
        call mqc_get_command_argument(i,command)
        if(command.eq.'-f') then
!
!*      -f matrix_file                   Input matrix file containing density.
!*
          call mqc_get_command_argument(i+1,fileName)
          j = i+2
        elseIf(command.eq.'--distance') then
!
!*      --distance matrix_file           Additional matrix file containing MO coefficients
!*                                       and/or density from which distance of first matrix
!*                                       file will be computed.
!*
          call mqc_get_command_argument(i+1,fileName2)
          j = i+2
        elseIf(command.eq.'--print-level') then
!
!*      --print-level print_level        Output print level value from 0 (minimum) to 3 
!*                                       (maximum). Default is 0.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I2)') iPrint
          j = i+2
        elseIf(command.eq.'--help') then 
!
!*      --help                           Output help documentation to terminal. 
!*
          if(command_argument_count().gt.1) call mqc_error_I('Help output requested with multiple arguments',6, &
            'command_argument_count()',command_argument_count()) 

          call mqc_get_command_argument(0,help_path)
          help_path = 'less ' // trim(help_path(1:scan(help_path,'/',.true.))) // 'doc/denTest.txt'
          call execute_command_line(help_path,exitstat=flag)
          if(flag.ne.0) call mqc_error('Help output command failed')
          stop
        else
          call mqc_error_A('Unrecognised input flag',6,'command',command)
        endIf
        deallocate(command)
      endDo
!
      call fileInfo%load(filename)

      nBasis = fileInfo%getVal('nBasis')
      nAlpha = fileInfo%getVal('nAlpha')
      nBeta = fileInfo%getVal('nBeta')

      call fileInfo%getESTObj('overlap',est_integral=overlap)
      call fileInfo%getESTObj('mo coefficients',est_integral=mo,foundObj=foundMO)
      call fileInfo%getESTObj('density',est_integral=density,foundObj=foundDen)
      if(.not.foundDen) then
        if(.not.foundMO) call mqc_error('Density or MOs must be on matrix file') 
        density = matmul(mo%orbitals('occupied',[nAlpha],[nBeta]), &
          dagger(mo%orbitals('occupied',[nAlpha],[nBeta])))
      endIf

      write(iOut,'(A,A)') ' Testing symmetry of density on file ',trim(filename) 
      call density_test(iOut,iPrint,density,overlap)

      if(allocated(fileName2)) then
        call fileInfo%load(filename2)
        call fileInfo%getESTObj('overlap',est_integral=overlap2)
        if(abs(mqc_integral_norm(overlap-overlap2)).gt.thresh) &
          call mqc_error('Distance cannot be computed between different AO basis')
        call fileInfo%getESTObj('mo coefficients',est_integral=mo2)
        call density_distance(iOut,iPrint,nAlpha,nBeta,nBasis,overlap,MO,MO2,thresh,&
          mo_overlap,density_dist)
        write(iOut,'(A,A,A,A)') ' Computing distance between solutions on files ',trim(filename),&
          ' and ',trim(filename2) 
        call mqc_print(one-mo_overlap,iOut,'     Normed MO distance (1-S)')
        call mqc_print(one-abs(mo_overlap),iOut,'     Normed density distance (1-|S|)')
        call mqc_print(density_dist,iOut,'     Density distance metric (N-<PSPS>)')
      endIf

!
      contains
!
      subroutine density_test(iOut,iPrint,density_integral,overlap)

      implicit none
      type(mqc_scf_integral),intent(in)::density_integral,overlap
      integer,intent(in)::iOut,iPrint
      integer::k,density_T_test,density_tao_test
      type(mqc_scalar)::temp_T,temp_tao,half,halfim,threshzero
      type(mqc_matrix)::P_den,Mx,My,Mz,Tmat_inc,TaoMat_inc
      type(mqc_vector)::Tvals_inc,TaoVals_inc

      threshzero = 1.0e-8
      half = 0.5
      halfim=-1*(0.0,0.5)

      P_den=half*(density_integral%getBlock('alpha')+density_integral%getBlock('beta'))
      Mx=half*(density_integral%getBlock('alpha-beta')+density_integral%getBlock('beta-alpha'))
      if(iPrint.ge.2) call Mx%print(6,'Mx')
      My=halfim*(density_integral%getBlock('alpha-beta')-density_integral%getBlock('beta-alpha'))
      if(iPrint.ge.2) call My%print(6,'My')
      Mz=half*(density_integral%getBlock('alpha')-density_integral%getBlock('beta'))
      if(iPrint.ge.2) call Mz%print(6,'Mz')
      call Tmat_inc%init(3,3)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(Mx,Mx),overlap%getBlock('alpha'))),1,1)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(Mx,My),overlap%getBlock('alpha'))),1,2)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(Mx,Mz),overlap%getBlock('alpha'))),1,3)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(My,Mx),overlap%getBlock('alpha'))),2,1)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(My,My),overlap%getBlock('alpha'))),2,2)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(My,Mz),overlap%getBlock('alpha'))),2,3)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(Mz,Mx),overlap%getBlock('alpha'))),3,1)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(Mz,My),overlap%getBlock('alpha'))),3,2)
      call Tmat_inc%put(mqc_matrix_trace(matmul(matmul(Mz,Mz),overlap%getBlock('alpha'))),3,3)
      If(iPrint.ge.3) call Tmat_inc%print(iOut,'T matrix')
      call Tmat_inc%diag(Tvals_inc)
      if(iPrint.ge.2) call Tvals_inc%print(6,'T test eigenvalues')
      
      density_T_test = 0
      do k = 1,3
        temp_T=Tvals_inc%at(k)
        if(temp_T%abs().lt.threshzero) density_T_test = density_T_test + 1
      endDo
      if(iPrint.ge.1) write(iOut,'(A,I2)') 'Number of zero T eigenvalues:',density_T_test
      density_tao_test=0
      if(density_T_test.le.1) then
        Mx = MQC_Matrix_Cast_Real(Mx)
        My = MQC_Matrix_Cast_Real(My)
        Mz = MQC_Matrix_Cast_Real(Mz)
        call TaoMat_inc%init(3,3)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(Mx,Mx),overlap%getBlock('alpha'))),1,1)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(Mx,My),overlap%getBlock('alpha'))),1,2)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(Mx,Mz),overlap%getBlock('alpha'))),1,3)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(My,Mx),overlap%getBlock('alpha'))),2,1)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(My,My),overlap%getBlock('alpha'))),2,2)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(My,Mz),overlap%getBlock('alpha'))),2,3)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(Mz,Mx),overlap%getBlock('alpha'))),3,1)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(Mz,My),overlap%getBlock('alpha'))),3,2)
        call TaoMat_inc%put(mqc_matrix_trace(matmul(matmul(Mz,Mz),overlap%getBlock('alpha'))),3,3)
        If(iPrint.ge.3) call TaoMat_inc%print(iOut,'Tao matrix')
        call TaoMat_inc%diag(TaoVals_inc)
        if(iPrint.ge.2) call TaoVals_inc%print(6,'Tao test eigenvalues')
        do k = 1,3
          temp_tao=TaoVals_inc%at(k)
          if(temp_tao%abs().lt.threshzero) density_tao_test = density_tao_test + 1
        endDo
        if(iPrint.ge.1) write(iOut,'(A,I2)')'Number of zero tao eigenvalues: ',density_tao_test 
      endIf

      select case (density_T_test)
      case (3)
        write(iOut,'(A)') ' Solution is fully symmetric based on T-test'
      case (2)
        write(iOut,'(A)') ' Solution is collinear based on T-test'
      case (0:1)
        write(iOut,'(A)') ' Solution is noncolinear based on T-test'
        select case (density_tao_test)
        case (1:3)
          write(iOut,'(A)') ' Solution is coplanar based on tao-test'
        case (0)
          write(iOut,'(A)') ' Solution is noncoplanar based on tao-test'
        case default
          call mqc_error_I('Unexpected number in density_tao_test',6,'density_tao_test',density_tao_test)
        end select
      case default
        call mqc_error_I('Unexpected number in density_T_test',6,'density_T_test',density_T_test)
      end select

      end subroutine density_test  
!
      subroutine density_distance(iOut,iPrint,nAlpha,nBeta,nBasis,overlap,MO1,MO2,thresh,&
        mo_overlap,density_dist)

      implicit none
      type(mqc_scf_integral),intent(in)::overlap
      type(mqc_scf_integral),intent(inOut)::MO1,MO2
      integer,intent(in)::iOut,iPrint,nAlpha,nBeta,nBasis
      real(kind=real64),intent(in)::thresh
      type(mqc_scalar),intent(inOut)::mo_overlap,density_dist
      logical::have_coplanar,have_noncoplanar
      type(mqc_scf_integral)::density1,density2

      If(mqc_matrix_norm(MO1%getBlock('alpha-beta')).gt.thresh.or. &
        mqc_matrix_norm(MO1%getBlock('beta-alpha')).gt.thresh.or. &
        mqc_matrix_norm(MO2%getBlock('alpha-beta')).gt.thresh.or. &
        mqc_matrix_norm(MO2%getBlock('beta-alpha')).gt.thresh) then
        If(mqc_matrix_norm(aimag(MO1%getBlock('full'))).gt.thresh.or. &
          mqc_matrix_norm(aimag(MO2%getBlock('full'))).gt.thresh) then
          have_coplanar = .false.
          have_noncoplanar = .true.
        else
          have_coplanar = .true.
          have_noncoplanar = .false.
        endIf
      else
        have_coplanar = .false.
        have_noncoplanar = .false.
      endIf

      MO1 = MO1%orbitals('occupied',[nAlpha],[nBeta])
      MO2 = MO2%orbitals('occupied',[nAlpha],[nBeta])

      if(have_coplanar) then
        call coplanar_rotate(iOut,iPrint,nBasis,MO1,MO2,overlap)
      elseIf(have_noncoplanar) then
        call noncoplanar_rotate(iOut,iPrint,nBasis,MO1,MO2,overlap)
      endIf

      mo_overlap = mqc_scf_integral_determinant(matmul(dagger(mo1),matmul(overlap,mo2)))

      density1 = matmul(mo1,dagger(mo1))
      density2 = matmul(mo2,dagger(mo2))

      density_dist = nAlpha + nBeta - &
        mqc_scf_integral_trace(matmul(matmul(density1,overlap),matmul(density2,overlap))) 

      end subroutine density_distance
!
      subroutine coplanar_rotate(iOut,iPrint,nBasis,MOs_rotate,MOs_fixed,overlap)

      implicit none
      type(mqc_scf_integral),intent(inOut)::MOs_rotate
      type(mqc_scf_integral),intent(in)::MOs_fixed,overlap
      integer,intent(in)::iOut,iPrint,nBasis
      type(mqc_vector)::GHFrot_thetas,GHFrot_norms
      type(mqc_matrix)::GHFrot_R
      type(mqc_scalar)::threshzero,GHFrot_angle,GHFrot_step,GHFrot_err
      integer::k,GHFrot_iters=150,GHFrot_count

      threshzero = 1.0e-8 

      GHFrot_thetas = [0.0,pi/2,2*pi]
      call GHFrot_norms%init(3,(100.0,0.0))
      do k=1,3
        GHFrot_angle = GHFrot_thetas%at(k)
        if(iPrint.ge.3) call GHFrot_angle%print(iout,'GHFrot_theta')
        call GHF_rotation(GHFrot_R,GHFrot_angle,nBasis)
        if(iPrint.ge.3) call GHFrot_R%print(iOut,'Rotation matrix')
        if(iPrint.ge.3) call mqc_print(matmul(GHFrot_R,MOs_rotate),iout,'Rotated MOs_rotate')
        call GHFrot_norms%put(abs(mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)), &
          overlap),MOs_fixed))),k)
      endDo
      GHFrot_step = ((3-sqrt(5.0))/2)
      GHFrot_count = 1
      do while (GHFrot_count.le.GHFrot_iters)
        if(iPrint.ge.3) then
          call GHFrot_thetas%print(iout,'GHF theta list')
          call GHFrot_norms%print(iout,'GHF norm list')
        endIf
        if(abs(GHFrot_thetas%at(2)-GHFrot_thetas%at(1)).ge. &
          abs(GHFrot_thetas%at(2)-GHFrot_thetas%at(3))) then
          GHFrot_angle = GHFrot_thetas%at(2) + &
            GHFrot_step*(GHFrot_thetas%at(2)-GHFrot_thetas%at(1))
          if(iPrint.ge.3) call GHFrot_angle%print(iout,'new theta')
          call GHF_rotation(GHFrot_R,GHFrot_angle,nBasis)
!          GHFrot_err = mqc_integral_norm(MOs_fixed - matmul(GHFrot_R,MOs_rotate))
          GHFrot_err = abs(mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)),overlap), &
            MOs_fixed)))
          if(iPrint.ge.3) call GHFrot_err%print(iOut,'new norm')
          if(GHFrot_err.gt.GHFrot_norms%at(2)) then
            call GHFrot_thetas%put(GHFrot_thetas%at(2),3)
            call GHFrot_thetas%put(GHFrot_angle,2)
            call GHFrot_norms%put(GHFrot_norms%at(2),3)
            call GHFrot_norms%put(GHFrot_err,2)
          else
            call GHFrot_thetas%put(GHFrot_angle,1)
            call GHFrot_norms%put(GHFrot_err,1)
          endIf
        else
          GHFrot_angle = GHFrot_thetas%at(2) + &
            GHFrot_step*(GHFrot_thetas%at(3)-GHFrot_thetas%at(2))
          if(iPrint.ge.3) call GHFrot_angle%print(iout,'new theta')
          call GHF_rotation(GHFrot_R,GHFrot_angle,nBasis)
!          GHFrot_err = mqc_integral_norm(MOs_fixed - matmul(GHFrot_R,MOs_rotate))
          GHFrot_err = abs(mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)), &
            overlap),MOs_fixed)))
          if(iPrint.ge.3) call GHFrot_err%print(iOut,'new norm')
          if(GHFrot_err.gt.GHFrot_norms%at(2)) then
            call GHFrot_thetas%put(GHFrot_thetas%at(2),1)
            call GHFrot_thetas%put(GHFrot_angle,2)
            call GHFrot_norms%put(GHFrot_norms%at(2),1)
            call GHFrot_norms%put(GHFrot_err,2)
          else
            call GHFrot_thetas%put(GHFrot_angle,3)
            call GHFrot_norms%put(GHFrot_err,3)
          endIf
        endIf
        GHFrot_count = GHFrot_count + 1
        if(((abs(GHFrot_norms%at(1)-GHFrot_norms%at(2)).lt.threshzero).and. &
          (abs(GHFrot_norms%at(2)-GHFrot_norms%at(3)).lt.threshzero)).or. &
          (GHFrot_norms%at(2).lt.threshzero)) then
          if(iPrint.ge.1) write(iOut,'(A)') 'GHF orbital rotation converged'
          if(iPrint.ge.1) write(iOut,'(A)') 'Setting MOs to common orientation'
          if(iPrint.ge.2)call MOs_rotate%print(iout,'Initial MOs_rotate')
          MOs_rotate = matmul(GHFrot_R,MOs_rotate)
          if(iPrint.ge.2) call MOs_rotate%print(iout,'Rotated MOs_rotate')
          if(iPrint.ge.1) call MOs_fixed%print(iout,'Compare to MOs_fixed')
          exit
        endIf
        if(GHFrot_count.eq.GHFrot_iters) then
          if(iPrint.ge.1) write(iOut,'(A)') 'GHF orbital rotation did not converge'
          if(iPrint.ge.1) write(iOut,'(A)') 'Continuing with input MO orientation'
          exit
        endIf
      endDo

      end subroutine coplanar_rotate 
!
!
      subroutine noncoplanar_rotate(iOut,iPrint,nBasis,MOs_rotate,MOs_fixed,overlap)

      implicit none
      type(mqc_scf_integral),intent(inOut)::MOs_rotate
      type(mqc_scf_integral),intent(in)::MOs_fixed,overlap
      integer,intent(in)::iOut,iPrint,nBasis
      type(mqc_vector)::GHFrot_thetas,GHFrot_norms,GHFrot_thetas_2,GHFrot_thetas_3,idx
      type(mqc_matrix)::GHFrot_R,GHFrot_R_int
      type(mqc_scalar)::threshzero,GHFrot_angle,two,three,mid_i,mid_j,mid_k,r_i,r_j,r_k,r_norm, &
        e_i,e_j,e_k,e_norm,c_i,c_j,c_k,c_norm
      integer::k,GHFrot_iters=200,GHFrot_count

      threshzero = 1.0e-6 
      two = 2.0
      three = 3.0

      GHFrot_thetas = [pi,3*pi/2,pi,pi]
      GHFrot_thetas_2 = [pi,pi,3*pi/2,pi]
      GHFrot_thetas_3 = [pi,pi,pi,3*pi/2]
      call GHFrot_norms%init(4,(100.0,0.0))
      do k=1,4
        call GHFrot_R%identity(nBasis*2,nBasis*2)
        if(iPrint.ge.3) call GHFrot_R%print(iOut,'Rotation matrix at beginning')
        GHFrot_angle = GHFrot_thetas%at(k)
        if(iPrint.ge.3) call GHFrot_angle%print(iout,'GHFrot_theta')
        call complex_rotation(GHFrot_R_int,GHFrot_angle,nBasis)
        if(iPrint.ge.3) call GHFrot_R_int%print(iOut,'Intermediate rotation matrix')
        GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
        if(iPrint.ge.3) call GHFrot_R%print(iOut,'Rotation matrix at 1')
        GHFrot_angle = GHFrot_thetas_2%at(k)
        if(iPrint.ge.3) call GHFrot_angle%print(iout,'GHFrot_theta_2')
        call GHF_rotation(GHFrot_R_int,GHFrot_angle,nBasis)
        if(iPrint.ge.3) call GHFrot_R_int%print(iOut,'Intermediate rotation matrix')
        GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
        if(iPrint.ge.3) call GHFrot_R%print(iOut,'Rotation matrix at 2')
        GHFrot_angle = GHFrot_thetas_3%at(k)
        if(iPrint.ge.3) call GHFrot_angle%print(iout,'GHFrot_theta_3')
        call complex_rotation(GHFrot_R_int,GHFrot_angle,nBasis)
        if(iPrint.ge.3) call GHFrot_R_int%print(iOut,'Intermediate rotation matrix')
        GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
        if(iPrint.ge.3) call GHFrot_R%print(iOut,'Final rotation matrix')
        if(iPrint.ge.3) call mqc_print(matmul(GHFrot_R,MOs_rotate),iout,'Rotated MOs_rotate')
        call GHFrot_norms%put(mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)), &
          overlap),MOs_fixed)),k)
      endDo
      GHFrot_count = 1
      do while (GHFrot_count.le.GHFrot_iters)
        if(iPrint.ge.3) then
          call GHFrot_thetas%print(iout,'GHF theta 1 list')
          call GHFrot_thetas_2%print(iout,'GHF theta 2 list')
          call GHFrot_thetas_3%print(iout,'GHF theta 3 list')
          call GHFrot_norms%print(iout,'GHF norm list')
        endIf
        idx = GHFrot_norms%argsort()
        if(iPrint.ge.3) call idx%print(iout,'GHF norm idx')
        call GHFrot_norms%sort(idx)
        if(iPrint.ge.3) call GHFrot_norms%print(iout,'GHF norm sorted')
        call GHFrot_thetas%sort(idx)
        if(iPrint.ge.3) call GHFrot_thetas%print(iout,'GHF thetas sorted')
        call GHFrot_thetas_2%sort(idx)
        if(iPrint.ge.3) call GHFrot_thetas_2%print(iout,'GHF thetas 2 sorted')
        call GHFrot_thetas_3%sort(idx)
        if(iPrint.ge.3) call GHFrot_thetas_3%print(iout,'GHF thetas 3 sorted')
        mid_i = (GHFrot_thetas%at(2)+GHFrot_thetas%at(3)+GHFrot_thetas%at(4))/three
        mid_j = (GHFrot_thetas_2%at(2)+GHFrot_thetas_2%at(3)+GHFrot_thetas_2%at(4))/three
        mid_k = (GHFrot_thetas_3%at(2)+GHFrot_thetas_3%at(3)+GHFrot_thetas_3%at(4))/three
        if(iPrint.ge.3) call mid_i%print(iout,'mid i')
        if(iPrint.ge.3) call mid_j%print(iout,'mid j')
        if(iPrint.ge.3) call mid_k%print(iout,'mid k')
        r_i = two*mid_i - GHFrot_thetas%at(1)
        r_j = two*mid_j - GHFrot_thetas_2%at(1)
        r_k = two*mid_k - GHFrot_thetas_3%at(1)
        if(iPrint.ge.3) call r_i%print(iout,'r i')
        if(iPrint.ge.3) call r_j%print(iout,'r j')
        if(iPrint.ge.3) call r_k%print(iout,'r k')
        call GHFrot_R%identity(nBasis*2,nBasis*2)
        call complex_rotation(GHFrot_R_int,r_i,nBasis)
        GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
        call GHF_rotation(GHFrot_R_int,r_j,nBasis)
        GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
        call complex_rotation(GHFrot_R_int,r_k,nBasis)
        GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
        r_norm = mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)), &
          overlap),MOs_fixed))
        if(iPrint.ge.3) call r_norm%print(iout,'r norm')
        if(r_norm.gt.GHFrot_norms%at(2).and.r_norm.lt.GHFrot_norms%at(3)) then
          call GHFrot_thetas%put(r_i,1)
          call GHFrot_thetas_2%put(r_j,1)
          call GHFrot_thetas_3%put(r_k,1)
          call GHFrot_norms%put(r_norm,1)
        elseIf(r_norm.gt.GHFrot_norms%at(3)) then
          e_i = two*r_i - mid_i
          e_j = two*r_j - mid_j
          e_k = two*r_k - mid_k
          if(iPrint.ge.3) call e_i%print(iout,'e i')
          if(iPrint.ge.3) call e_j%print(iout,'e j')
          if(iPrint.ge.3) call e_j%print(iout,'e k')
          call GHFrot_R%identity(nBasis*2,nBasis*2)
          call complex_rotation(GHFrot_R_int,e_i,nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          call GHF_rotation(GHFrot_R_int,e_j,nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          call complex_rotation(GHFrot_R_int,e_k,nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          e_norm = mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)), &
            overlap),MOs_fixed))
          if(iPrint.ge.3) call e_norm%print(iout,'e norm')
          if(e_norm.gt.r_norm) then
            call GHFrot_thetas%put(e_i,1)
            call GHFrot_thetas_2%put(e_j,1)
            call GHFrot_thetas_3%put(e_k,1)
            call GHFrot_norms%put(e_norm,1)
          else
            call GHFrot_thetas%put(r_i,1)
            call GHFrot_thetas_2%put(r_j,1)
            call GHFrot_thetas_3%put(r_k,1)
            call GHFrot_norms%put(r_norm,1)
          endIf
        else
          c_i = (mid_i+GHFrot_thetas%at(1))/two
          c_j = (mid_j+GHFrot_thetas_2%at(1))/two
          c_k = (mid_k+GHFrot_thetas_3%at(1))/two
          if(iPrint.ge.3) call c_i%print(iout,'c i')
          if(iPrint.ge.3) call c_j%print(iout,'c j')
          if(iPrint.ge.3) call c_k%print(iout,'c k')
          call GHFrot_R%identity(nBasis*2,nBasis*2)
          call complex_rotation(GHFrot_R_int,c_i,nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          call GHF_rotation(GHFrot_R_int,c_j,nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          call complex_rotation(GHFrot_R_int,c_k,nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          c_norm = mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)), &
            overlap),MOs_fixed))
          if(iPrint.ge.3) call c_norm%print(iout,'c norm')
          if(c_norm.gt.GHFrot_norms%at(1)) then
            call GHFrot_thetas%put(c_i,1)
            call GHFrot_thetas_2%put(c_j,1)
            call GHFrot_thetas_3%put(c_k,1)
            call GHFrot_norms%put(c_norm,1)
          else
            call GHFrot_thetas%put((GHFrot_thetas%at(3)+GHFrot_thetas%at(1))/two,1)
            call GHFrot_thetas_2%put((GHFrot_thetas_2%at(3)+GHFrot_thetas_2%at(1))/two,1)
            call GHFrot_thetas_3%put((GHFrot_thetas_3%at(3)+GHFrot_thetas_3%at(1))/two,1)
            call GHFrot_R%identity(nBasis*2,nBasis*2)
            call complex_rotation(GHFrot_R_int,GHFrot_thetas%at(1),nBasis)
            GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
            call GHF_rotation(GHFrot_R_int,GHFrot_thetas_2%at(1),nBasis)
            GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
            call complex_rotation(GHFrot_R_int,GHFrot_thetas_3%at(1),nBasis)
            GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
            call GHFrot_norms%put(mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)), &
            overlap),MOs_fixed)),1)
            call GHFrot_thetas%put((GHFrot_thetas%at(3)+GHFrot_thetas%at(2))/two,2)
            call GHFrot_thetas_2%put((GHFrot_thetas_2%at(3)+GHFrot_thetas_2%at(2))/two,2)
            call GHFrot_thetas_3%put((GHFrot_thetas_3%at(3)+GHFrot_thetas_3%at(2))/two,2)
            call GHFrot_R%identity(nBasis*2,nBasis*2)
            call complex_rotation(GHFrot_R_int,GHFrot_thetas%at(2),nBasis)
            GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
            call GHF_rotation(GHFrot_R_int,GHFrot_thetas_2%at(2),nBasis)
            GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
            call complex_rotation(GHFrot_R_int,GHFrot_thetas_3%at(2),nBasis)
            GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
            call GHFrot_norms%put(mqc_scf_integral_determinant(matmul(matmul(dagger(matmul(GHFrot_R,MOs_rotate)), &
            overlap),MOs_fixed)),2)
          endIf
        endIf
        GHFrot_count = GHFrot_count + 1
        if((abs(GHFrot_thetas%at(1)-GHFrot_thetas%at(2)).lt.threshzero).and. &
          (abs(GHFrot_thetas%at(1)-GHFrot_thetas%at(3)).lt.threshzero).and. &
          (abs(GHFrot_thetas%at(1)-GHFrot_thetas%at(4)).lt.threshzero).and. &
          (abs(GHFrot_thetas_2%at(1)-GHFrot_thetas_2%at(2)).lt.threshzero).and. &
          (abs(GHFrot_thetas_2%at(1)-GHFrot_thetas_2%at(3)).lt.threshzero).and. &
          (abs(GHFrot_thetas_2%at(1)-GHFrot_thetas_2%at(4)).lt.threshzero).and. & 
          (abs(GHFrot_thetas_3%at(1)-GHFrot_thetas_3%at(2)).lt.threshzero).and. &
          (abs(GHFrot_thetas_3%at(1)-GHFrot_thetas_3%at(3)).lt.threshzero).and. &
          (abs(GHFrot_thetas_3%at(1)-GHFrot_thetas_3%at(4)).lt.threshzero)) then
          if(iPrint.ge.1) write(iOut,'(A)') 'GHF orbital rotation converged'
          if(iPrint.ge.1) write(iOut,'(A)') 'Setting MOs to common orientation'
          if(iPrint.ge.2) call MOs_rotate%print(iout,'Initial MOS_rotate')
          call GHFrot_R%identity(nBasis*2,nBasis*2)
          call complex_rotation(GHFrot_R_int,GHFrot_thetas%at(3),nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          call GHF_rotation(GHFrot_R_int,GHFrot_thetas_2%at(3),nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          call complex_rotation(GHFrot_R_int,GHFrot_thetas_3%at(3),nBasis)
          GHFrot_R = matmul(GHFrot_R_int,GHFrot_R)
          MOs_rotate = matmul(GHFrot_R,MOs_rotate)
          if(iPrint.ge.2) call MOs_rotate%print(iout,'Rotated MOS_rotate')
          if(iPrint.ge.1) call MOs_fixed%print(iout,'Compare to MOs_fixed')
          exit
        endIf
        if(GHFrot_count.eq.GHFrot_iters) then
          if(iPrint.ge.1) write(iOut,'(A)') 'GHF/complex orbital rotation did not converge'
          if(iPrint.ge.1) write(iOut,'(A)') 'Continuing with input MO orientation'
          exit
        endIf
      endDo
      
      end subroutine noncoplanar_rotate 
!
      subroutine GHF_rotation(rotation_matrix,theta,nBasis)
!
      implicit none
      type(mqc_matrix)::rotation_matrix
      type(mqc_scalar)::theta,cangle,sangle,nsangle
      integer::nBasis,i
!
      cangle = cos(theta%rval())
      sangle = sin(theta%rval())
      nsangle = -sin(theta%rval())
      call rotation_matrix%identity(nBasis*2,nBasis*2)
      do i = 1, nBasis
        call rotation_matrix%put(cangle,i,i)
        call rotation_matrix%put(cangle,i+nBasis,i+nBasis)
        call rotation_matrix%put(sangle,i,i+nBasis)
        call rotation_matrix%put(nsangle,i+nBasis,i)
      endDo
!
      return
      end subroutine GHF_rotation
!
      subroutine complex_rotation(rotation_matrix,theta,nBasis)
!
      implicit none
      type(mqc_matrix)::rotation_matrix
      type(mqc_scalar)::theta,eangle,neangle
      integer::nBasis,i
!
      eangle = exp(cmplx(0,theta%rval()))
      neangle = exp(cmplx(0,-theta%rval()))
      call rotation_matrix%identity(nBasis*2,nBasis*2)
      do i = 1, nBasis
        call rotation_matrix%put(eangle,i,i)
        call rotation_matrix%put(neangle,i+nBasis,i+nBasis)
      endDo
!
      return
      end subroutine complex_rotation
!
!
!*    NOTES
!*      Compilation of this program requires the MQC library (https://github.com/MQCPack/mqcPack) 
!*      and the gauopen utility (http://gaussian.com/g16/gauopen.zip) and compilation with the f08
!*      standard. 
!*
!*      Compilation tested using: pgfortran, gfortran. e.g.
!*        gfortran -fdefault-integer-8 -fdefault-real-8 -std=f2008 -I/path/to/mqc/src -o guessGen.exe
!*        /path/to/mqc/src/*.o guessGen.f03 -llapack -lblas /path/to/mqc/src/libmqc.a
!*      
!*      Note that subroutine Wr_LCBuf needs modifying in gauopen/qcmatrix.F as follows:
!*        line 58:       LenBX = (LenBuf/(2*abs(NR)))*abs(NR)
!*        line 60:       Call Wr_CBuf(IU,NTot*abs(NR),LenBX,X)
!*
!*      Documentation generated with robodoc. To update documentation edit robodoc.rc to determine
!*      documentation output type and then run robodoc at the command line in the guessGen directory.
!*
!*    AUTHOR
!*      Lee Thompson, University of Louisville, lee.thompson.1@lousiville.edu
!*
!*    COPYRIGHT 
!*      (c) 2022 by Lee Thompson distributed under terms of the MIT license.
!*
!****
!     
 
 999  end program denTest
