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
      character(len=:),allocatable::command,fileName,help_path
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      integer::iOut=output_unit,iPrint=0,i,j,flag
      type(mqc_scf_integral)::density,overlap
!
!*    USAGE
!*      denTest [-f <Matrix_file>] [--print-level <print_level>] [--help]
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
      call fileInfo%getESTObj('overlap',est_integral=overlap)
      call fileInfo%getESTObj('density',est_integral=density)
      call density_test(iOut,iPrint,density,overlap)
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
        write(iOut,'(A)') 'Solution is fully symmetric based on T-test'
      case (2)
        write(iOut,'(A)') 'Solution is collinear based on T-test'
      case (0:1)
        write(iOut,'(A)') 'Solution is noncolinear based on T-test'
        select case (density_tao_test)
        case (1:3)
          write(iOut,'(A)') 'Solution is coplanar based on tao-test'
        case (0)
          write(iOut,'(A)') 'Solution is noncoplanar based on tao-test'
        case default
          call mqc_error_I('Unexpected number in density_tao_test',6,'density_tao_test',density_tao_test)
        end select
      case default
        call mqc_error_I('Unexpected number in density_T_test',6,'density_T_test',density_T_test)
      end select

      end subroutine density_test  
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
