module constants
    implicit none
    real(8), parameter :: pi = 3.14159265358979323846264338327d0
end module constants

module fftwa3_mpi
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        ! FFTW3 MPI functions
        function fftw_mpi_local_size_2d_transposed(n0, n1, comm, lH, lh0, & 
lW, lw0) bind(C, name="fftw_mpi_local_size_2d_transposed")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: fftw_mpi_local_size_2d_transposed
            integer(c_int), value :: n0, n1
            integer(c_int), value :: comm
            integer(c_ptrdiff_t) :: lH, lh0, lW, lw0
        end function fftw_mpi_local_size_2d_transposed

       function fftw_mpi_plan_dft_r2c_2d(n0, n1, in, out, comm, flags) &
 bind(C, name="fftw_mpi_plan_dft_r2c_2d")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: fftw_mpi_plan_dft_r2c_2d
            integer(c_int), value :: n0, n1
            type(c_ptr), value :: in, out
            integer(c_int), value :: comm, flags
        end function fftw_mpi_plan_dft_r2c_2d

        function fftw_mpi_plan_dft_c2r_2d(n0, n1, in, out, comm, flags) & 
bind(C, name="fftw_mpi_plan_dft_c2r_2d")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: fftw_mpi_plan_dft_c2r_2d
            integer(c_int), value :: n0, n1
            type(c_ptr), value :: in, out
            integer(c_int), value :: comm, flags
        end function fftw_mpi_plan_dft_c2r_2d

        function fftw_malloc(size) bind(C, name='fftw_malloc')
            use, intrinsic :: iso_c_binding
            implicit none
            integer(C_SIZE_T), value :: size
            type(c_ptr) :: fftw_malloc
        end function fftw_malloc
    end interface

end module fftwa3_mpi

module data_types
    use, intrinsic :: iso_c_binding
    implicit none

    type :: Output
        character(len=128) :: name       ! run name
        integer :: T_print               ! interval for printing output
        integer :: T_write               ! interval for saving state
    end type Output

    type :: Arrays
        integer :: W                     ! system width
        integer :: H                     ! system height
        integer(kind=c_ptrdiff_t) :: lH  ! local system height
        integer(kind=c_ptrdiff_t) :: lW  ! local system width
        integer(kind=c_ptrdiff_t) :: lh0 ! local vertical start index
        integer(kind=c_ptrdiff_t) :: lw0 ! local horizontal start index

        real(8) :: u_C_min               
! minimum of smoothed density (graphene)
        real(8) :: u_C_max               ! maximum
! minimum of smoothed density (HBN)
        real(8) :: u_BN_min              
        real(8) :: u_BN_max              ! maximum

        real(8), pointer :: A_C(:)       
! operator for linear part, e^{-k^2 \hat{\mathcal{L}} \Delta t} (graphene)
        real(8), pointer :: B_C(:)       ! operator for nonlinear part
        real(8), pointer :: p_arr_C(:)       ! array for \psi_C
        real(8), pointer :: q_arr_C(:)       ! another array
        real(8), pointer :: u_arr_C(:)       ! smoothed \psi_C
        real(8), pointer :: A_BN(:)      ! operator for linear part (HBN)
        real(8), pointer :: B_BN(:)      ! operator for nonlinear part
        real(8), pointer :: p_arr_B(:)       ! ...
        real(8), pointer :: q_arr_B(:)
        real(8), pointer :: p_arr_N(:)
        real(8), pointer :: q_arr_N(:)
        real(8), pointer :: u_arr_BN(:)

        integer :: p_P_C                 ! FFTW plan for F(p) (graphene)
        integer :: qa_Q_C                 ! F(q)
        integer :: Q_q_C                 ! F^-1(Q)
        integer :: p_P_B                 ! ... (boron)
        integer :: qa_Q_B                 ! ...
        integer :: Q_q_B
        integer :: p_P_N
        integer :: qa_Q_N
        integer :: Q_q_N
        integer :: ua_U_C
        integer :: U_u_C
        integer :: ua_U_BN
        integer :: U_u_BN

        type(c_ptr) :: A_C_cptr

        type(c_ptr) :: p_arr_C_cptr
        type(c_ptr) :: q_arr_C_cptr
        type(c_ptr) :: p_arr_B_cptr
        type(c_ptr) :: q_arr_B_cptr
        type(c_ptr) :: p_arr_N_cptr
        type(c_ptr) :: q_arr_N_cptr
        type(c_ptr) :: u_arr_C_cptr
        type(c_ptr) :: u_arr_BN_cptr
        type(c_ptr) :: B_C_cptr
        type(c_ptr) :: A_BN_cptr
        type(c_ptr) :: B_BN_cptr


    end type Arrays

    type :: Model
! see documentation and/or input file
        real(8) :: alpha_C               
        real(8) :: alpha_BN
        real(8) :: alpha_C_B
        real(8) :: alpha_B_C
        real(8) :: alpha_C_N
        real(8) :: alpha_N_C
        real(8) :: alpha_B_N
        real(8) :: beta_C
        real(8) :: beta_BN
        real(8) :: beta_B_N
        real(8) :: gamma_C_l
        real(8) :: gamma_C_s
        real(8) :: gamma_BN_l
        real(8) :: gamma_BN_s
        real(8) :: gamma_B_N
        real(8) :: gamma_N_B
        real(8) :: gamma_C_B
        real(8) :: gamma_B_C
        real(8) :: gamma_C_N
        real(8) :: gamma_N_C
        real(8) :: delta_C
        real(8) :: delta_BN
        real(8) :: l_BN
        real(8) :: sigma_u
        real(8) :: sigma_mask
    end type Model

    type :: Relaxation
        integer(kind=c_long) :: t0       ! start time
        integer :: iid                    ! rank of process
        integer :: ID                    ! number of processes
        integer :: tt                     ! time step count
        real(8) :: f_C                   ! free energy density (graphene)
        real(8) :: f_BN                  ! ... (HBN)
        real(8) :: p_rex_C                   ! average density (graphene)
        real(8) :: p_rex_BN
        real(8) :: d                     
! sampling step size for calculation box optimization

        integer :: T                     ! total number of iterations
        real(8) :: dx                    ! x-discretization
        real(8) :: dy                    ! y-discretization
        real(8) :: dt                    ! time step
        integer :: T_optimize            
! interval for calculation box optimization
    end type Relaxation

end module data_types

module utilities
    use, intrinsic :: iso_c_binding
    use constants
    use data_types
    use fftwa3_mpi
    !use mpi 
    !use mpi_f08
    
    implicit none
  !  include 'mpif.h'
    contains

    subroutine seed_rngs(relax, input)
        use, intrinsic :: iso_fortran_env, only: int64
        use iso_c_binding, only: c_f_pointer
        implicit none
        type(Relaxation), intent(inout) :: relax
        integer(c_int), intent(in) :: input
        integer :: MPI_INTEGER, MPI_COMM_WORLD, MPI_STATUS_SIZE

        integer, allocatable :: seed(:)
         integer :: ierr, status(MPI_STATUS_SIZE)
        character(len=20) :: str_seed

        ! Read seed from the input file
        read(input, *) seed

        if (relax%iid == 0) then
            if (seed == 0) then
                call random_seed()
            else
                call random_seed(put=seed)
            end if
            call random_number(seed)
            seed = int(seed * huge(seed))
        else
            call MPI_Recv(seed, 1, MPI_INTEGER, relax%iid-1, 0, &
MPI_COMM_WORLD, status, ierr)
            call random_seed(put=seed)
            call random_number(seed)
            seed = int(seed * huge(seed))
        end if

        if (relax%iid /= relax%ID-1) then
            call MPI_Send(seed, 1, MPI_INTEGER, relax%iid+1, 0, &
MPI_COMM_WORLD, ierr)
        end if
    end subroutine seed_rngs

    subroutine configure_arrays(arra, input)
        type(Arrays), intent(inout) :: arra
        integer(c_int), intent(in) :: input
        integer :: lWHh, lWHp, MPI_COMM_WORLD

        ! Read system width and height
        read(input, *) arra%W, arra%H

        ! Determine local data size
        lWHh = fftw_mpi_local_size_2d_transposed(arra%H, arra%W/2+1, & 
MPI_COMM_WORLD, arra%lH, arra%lh0, arra%lW, arra%lw0)
        lWHp = 2 * lWHh

        ! Allocate arrays
        arra%A_C_cptr = fftw_malloc(lWHh * c_sizeof(0.0d0))
        call c_f_pointer(arra%A_C_cptr, arra%A_C, [lWHh])

        arra%B_C_cptr = fftw_malloc(lWHh * sizeof(0.0d0))
        call c_f_pointer(arra% B_C_cptr, arra%B_C, [lWHh])

        arra%A_BN_cptr = fftw_malloc(lWHh * sizeof(0.0d0))
        call c_f_pointer(arra%A_BN_cptr, arra%A_BN, [lWHh])

        arra%B_BN_cptr = fftw_malloc(lWHh * sizeof(0.0d0))
        call c_f_pointer(arra%B_BN_cptr, arra%B_BN, [lWHh])

        arra%p_arr_C_cptr = fftw_malloc(lWHp * sizeof(0.0d0))
        call c_f_pointer(arra%p_arr_C_cptr, arra%p_arr_C, [lWHp])

        arra%q_arr_C_cptr = fftw_malloc(lWHp * sizeof(0.0d0))
        call c_f_pointer(arra%q_arr_C_cptr, arra%q_arr_C, [lWHp])

        arra%p_arr_B_cptr = fftw_malloc(lWHp * sizeof(0.0d0))
        call c_f_pointer(arra%p_arr_B_cptr, arra%p_arr_B, [lWHp])

        arra%q_arr_B_cptr = fftw_malloc(lWHp * sizeof(0.0d0))
        call c_f_pointer(arra%q_arr_B_cptr, arra%q_arr_B, [lWHp])

        arra%p_arr_N_cptr = fftw_malloc(lWHp * sizeof(0.0d0))
        call c_f_pointer(arra%p_arr_N_cptr, arra%p_arr_N, [lWHp])

        arra%q_arr_N_cptr = fftw_malloc(lWHp * sizeof(0.0d0))
        call c_f_pointer(arra%q_arr_N_cptr, arra%q_arr_N, [lWHp])

        arra%u_arr_C_cptr = fftw_malloc(lWHp * sizeof(0.0d0))
        call c_f_pointer(arra%u_arr_C_cptr, arra%u_arr_C, [lWHp])

        arra%u_arr_BN_cptr = fftw_malloc(lWHp * sizeof(0.0d0))
        call c_f_pointer(arra%u_arr_BN_cptr, arra%u_arr_BN, [lWHp])

        ! Set up FFTW plans
        arra%p_P_C = fftw_mpi_plan_dft_r2c_2d(arra%H, arra%W, &
arra%p_arr_C_cptr, arra%p_arr_C_cptr, MPI_COMM_WORLD, IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT))
        arra%qa_Q_C = fftw_mpi_plan_dft_r2c_2d(arra%H, arra%W, arra%q_arr_C_cptr, arra%q_arr_C_cptr, &
 MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_OUT)
        arra%Q_q_C = fftw_mpi_plan_dft_c2r_2d(arra%H, arra%W, arra%q_arr_C_cptr, arra%q_arr_C, &
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_IN)
        arra%p_P_B = fftw_mpi_plan_dft_r2c_2d(arra%H, arra%W, arra%p_arr_B_cptr, arra%p_arr_B_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_OUT)
        arra%qa_Q_B = fftw_mpi_plan_dft_r2c_2d(arra%H, arra%W, arra%q_arr_B_cptr, arra%q_arr_B_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_OUT)
        arra%Q_q_B = fftw_mpi_plan_dft_c2r_2d(arra%H, arra%W, arra%q_arr_B_cptr, arra%q_arr_B_cptr, &
 MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_IN)
        arra%p_P_N = fftw_mpi_plan_dft_r2c_2d(arra%H, arra%W, arra%p_arr_N_cptr, arra%p_arr_N_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_OUT)
        arra%q_Q_N = fftw_mpi_plan_dft_r2c_2d(arra%H, arra%W, arra%q_arr_N_cptr, arra%q_arr_N_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_OUT)
        arra%Q_q_N = fftw_mpi_plan_dft_c2r_2d(arra%H, arra%W, arra%q_arr_N_cptr, arra%q_arr_N_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_IN)

        arra%ua_U_C = fftw_mpi_plan_dft_r2c_2d(arra%H, arra%W, arra%u_arr_C_cptr, arra%u_arr_C_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_OUT)
        arra%U_u_C = fftw_mpi_plan_dft_c2r_2d(arra%H, arra%W, arra%u_arr_C, arra%u_arr_C_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_IN)
        arra%ua_U_BN = fftw_mpi_plan_dft_r2c_2d(arra%H, arra%W, arra%u_arr_BN_cptr, arra%u_arr_BN_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_OUT)
        arra%U_u_BN = fftw_mpi_plan_dft_c2r_2d(arra%H, arra%W, arra%u_arr_BN_cptr, arra%u_arr_BN_cptr, & 
MPI_COMM_WORLD, FFTW_MEASURE.or.FFTW_MPI_TRANSPOSED_IN)
    end subroutine configure_arrays

    subroutine configure_output(outp, input)
        type(Output), intent(inout) :: outp
        integer(c_int), intent(in) :: input

        ! Read intervals for printing and saving from input file
        read(input, *) outp%T_print, outp%T_write
    end subroutine configure_output

    subroutine configure_model(mmode, input)
        type(Model), intent(inout) :: mmode
        integer(c_int), intent(in) :: input

        ! Read model parameters from input file
        read(input, *) mmode%alpha_C, mmode%beta_C, mmode%gamma_C_l, & 
mmode%gamma_C_s, mmode%delta_C
        read(input, *) mmode%alpha_BN, mmode%beta_BN, mmode%gamma_BN_l, & 
mmode%gamma_BN_s, mmode%delta_BN
        read(input, *) mmode%alpha_B_N, mmode%beta_B_N, mmode%gamma_B_N
        read(input, *) mmode%alpha_C_B, mmode%gamma_C_B
        read(input, *) mmode%alpha_C_N, mmode%gamma_C_N
        read(input, *) mmode%alpha_B_C, mmode%gamma_B_C
        read(input, *) mmode%alpha_N_C, mmode%gamma_N_C
        read(input, *) mmode%l_BN
        read(input, *) mmode%sigma_u, mmode%sigma_mask
    end subroutine configure_model

    subroutine configure_relaxation(relax, input)
        type(Relaxation), intent(inout) :: relax
        integer(c_int), intent(in) :: input

        ! Read from input file
        read(input, *) relax%T, relax%dx, relax%dy, relax%dt, & 
relax%T_optimize
    end subroutine configure_relaxation

    function OMA(x0, y0, u, v, l0, theta) result(oma_a)
        real(8), intent(in) :: x0, y0, u, v, l0, theta
        real(8) :: qx, qy, cosine, sine, x, y
        real(8) :: oma_a

        qx = 2.0d0 * pi / l0       ! wave numbers
        qy = qx / sqrt(3.0d0)
        cosine = cos(theta)
        sine = sin(theta)
        x = cosine * x0 - sine * y0 + u       ! rotation matrix applied
        y = sine * x0 + cosine * y0 + v
        oma_a = cos(qx * x) * cos(qy * y) + 0.5d0 * cos(2.0d0 * qy * y)
    end function OMA

    subroutine embedded_crystallite(arra, dx, dy, l0, l_BN, R, phases, & 
p_C_l, p_C_s, A_C, theta_C, u_arr_C, v_C, p_BN_l, p_BN_s, A_BN, theta_BN, u_arr_BN, v_BN)
        use, intrinsic :: iso_c_binding
        type(Arrays), intent(inout) :: arra
        real(8), intent(in) :: dx, dy, l0, l_BN, R
        integer, intent(in) :: phases
        real(8), intent(in) :: p_C_l, p_C_s, A_C, u_arr_C, v_C
        real(8), intent(inout) :: theta_C, theta_BN
        real(8), intent(in) :: p_BN_l, p_BN_s, A_BN, u_arr_BN, v_BN

        integer :: Wp, w, h, gh, k
        real(8) :: deg2rad, R2, y2, x

        deg2rad = pi / 180.0d0
        theta_C = theta_C * deg2rad
        theta_BN = theta_BN * deg2rad
        R2 = R * R

        Wp = 2 * (arra%W / 2 + 1)   ! padded width

        do h = 1, arra%lH
            gh = arra%lh0 + h       ! global vertical index
            y2 = (gh - 0.5d0 * arra%H) * dy
            y2 = y2 * y2
            k = Wp * (h - 1)          ! local row start index

            do w = 1, arra%W
                x = (w - 0.5d0 * arra%W) * dx
                if (x * x + y2 < R2) then    ! crystallite
                    if (phases == 0) then    ! graphene
                        arra%q_arr_C(k + w) = p_C_s + A_C * OMA((w - 0.5d0 * arra%W) * dx, & 
(gh - 0.5d0 * arra%H) * dy, u_arr_C, v_C, l0, theta_C)
                        arra%q_arr_B(k + w) = p_BN_l
                        arra%q_arr_N(k + w) = p_BN_l
                    else                    ! HBN
                        arra%q_arr_C(k + w) = p_C_l
                        arra%q_arr_B(k + w) = p_BN_s + A_BN * OMA((w - 0.5d0 * arra%W) * dx, & 
(gh - 0.5d0 * arra%H) * dy, 0.5d0 * l_BN * l0 + u_arr_BN, -0.5d0 * l_BN * & 
l0 / sqrt(3.0d0) + v_BN, l_BN * l0, theta_BN)
                        arra%q_arr_N(k + w) = p_BN_s + A_BN * OMA((w - 0.5d0 * arra%W) * dx, & 
(gh - 0.5d0 * arra%H) * dy, 0.5d0 * l_BN * l0 + u_arr_BN, 0.5d0 * l_BN * l0 / sqrt(3.0d0) & 
+ v_BN, l_BN * l0, theta_BN)
                    end if
                else                        ! matrix
                    if (phases == 1) then    ! HBN
                        arra%q_arr_C(k + w) = p_C_l
                        arra%q_arr_B(k + w) = p_BN_s + A_BN * OMA((w - 0.5d0 * arra%W) * dx, & 
(gh - 0.5d0 * arra%H) * dy, 0.5d0 * l_BN * l0 + u_arr_BN, -0.5d0 * l_BN * l0 & 
/ sqrt(3.0d0) + v_BN, l_BN * l0, theta_BN)
                        arra%q_arr_N(k + w) = p_BN_s + A_BN * OMA((w - 0.5d0 * arra%W) * dx, & 
(gh - 0.5d0 * arra%H) * dy, 0.5d0 * l_BN * l0 + u_arr_BN, 0.5d0 * l_BN * l0 & 
/ sqrt(3.0d0) + v_BN, l_BN * l0, theta_BN)
                    else                    ! graphene
                        arra%q_arr_C(k + w) = p_C_s + A_C * OMA((w - 0.5d0 * arra%W) * dx, & 
(gh - 0.5d0 * arra%H) * dy, u_arr_C, v_C, l0, theta_C)
                        arra%q_arr_B(k + w) = p_BN_l
                        arra%q_arr_N(k + w) = p_BN_l
                    end if
                end if
            end do
        end do
    end subroutine embedded_crystallite

    ! Map x (min <= x <= max) linearly between 0 and 1
    real(8) function map(x, min, max)
        real(8), intent(in) :: x, min, max
        if (min == max) then
            map = min    ! avoid divide by zero
        else
            map = (x - min) / (max - min)
        end if
    end function map

    ! Save state into a file
    subroutine write_state(arra, relax, outp)
        type(Arrays), intent(in) :: arra
        type(Relaxation), intent(in) :: relax
        type(Output), intent(in) :: outp

        character(len=128) :: filename
        integer :: Wp, i, w, h, k
        real(8) :: v_C, v_BN, dv
        integer :: file, ierr, MPI_STATUS_SIZE
        integer :: status(MPI_STATUS_SIZE)

        integer :: MPI_COMM_WORLD 

        Wp = 2 * (arra%W / 2 + 1)
        write(filename, '(A, "_t:", I0, ".dat")') trim(outp%name), & 
relax%tt

        do i = 0, relax%ID - 1
            call MPI_Barrier(MPI_COMM_WORLD, ierr)    
! synchronize processes
            if (relax%iid == i) then
                if (relax%iid == 0) then
                    open(file, file=filename, status='unknown', &
action='write', iostat=ierr)
                else
                    open(file, file=filename, status='unknown', & 
position='append', iostat=ierr)
                end if
                do h = 1, arra%lH
                    k = Wp * (h - 1)
                    do w = 1, arra%W
                        v_C = map(arra%u_arr_BN(k + w), arra%u_BN_max, & 
arra%u_BN_min)
                        v_BN = map(arra%u_arr_C(k + w), arra%u_C_min, & 
arra%u_C_max)
                        dv = (v_C - v_BN) / (v_C + v_BN)
                        write(file, '(6E20.12)') arra%q_arr_C(k + w), & 
arra%q_arr_B(k + w), &
 arra%q_arr_N(k + w), arra%u_arr_C(k + w), arra%u_arr_BN(k + w), dv
                    end do
                end do
                close(file)
            end if
        end do
    end subroutine write_state

    ! Print output
    subroutine print_output(relax, outp)
        type(Relaxation), intent(in) :: relax
        type(Output), intent(in) :: outp

        character(len=128) :: filename
        integer :: file, ierr

        if (relax%iid == 0) then
            print *, relax%tt, int(time() - relax%t0), relax%dx, relax%dy, & 
                relax%f_C + relax%f_BN, relax%f_C, relax%f_BN, & 
relax%p_rex_C, relax%p_rex_BN
            write(filename, '(A, ".out")') trim(outp%name)
            open(file, file=filename, status='unknown', position='append', & 
iostat=ierr)
            write(file, '(I0, I0, 6E20.12)') relax%tt, int(time() - & 
relax%t0), relax%dx, relax%dy, &
                relax%f_C + relax%f_BN, relax%f_C, relax%f_BN, & 
relax%p_rex_C, relax%p_rex_BN
            close(file)
        end if
    end subroutine print_output

    ! Update operators for linear and nonlinear parts
    subroutine update_AB(arra, mmode, relax)
        type(Arrays), intent(inout) :: arra
        type(Model), intent(in) :: mmode
        type(Relaxation), intent(in) :: relax

        integer :: w, h, gw, k
        real(8) :: ky, kx2, k2, divWH, dkx, dky, d2, l, expl, q_BN2

        divWH = 1.0d0 / (arra%W * arra%H)
        dkx = 2.0d0 * pi / (relax%dx * arra%W)
        dky = 2.0d0 * pi / (relax%dy * arra%H)
        q_BN2 = 1.0d0 / (mmode%l_BN * mmode%l_BN)

        do w = 1, arra%lW
            gw = arra%lw0 + w
            kx2 = gw * dkx
            kx2 = kx2 * kx2
            k = arra%H * (w - 1)
            do h = 1, arra%H
                if (h <= arra%H / 2) then
                    ky = (h - 1) * dky
                else
                    ky = (h - 1 - arra%H) * dky
                end if
                k2 = kx2 + ky * ky
                ! graphene
                d2 = (1.0d0 - k2)
                d2 = d2 * d2
                l = mmode%alpha_C + mmode%beta_C * d2
                expl = exp(-k2 * l * relax%dt)
                arra%A_C(k + h) = expl
                if (l == 0.0d0) then
                    arra%B_C(k + h) = -k2 * relax%dt
                else
                    arra%B_C(k + h) = (expl - 1.0d0) / l
                end if
                arra%B_C(k + h) = arra%B_C(k + h) * divWH
                ! HBN
                d2 = (q_BN2 - k2)
                d2 = d2 * d2
                l = mmode%alpha_BN + mmode%beta_BN * d2
                expl = exp(-k2 * l * relax%dt)
                arra%A_BN(k + h) = expl
                if (l == 0.0d0) then
                    arra%B_BN(k + h) = -k2 * relax%dt
                else
                    arra%B_BN(k + h) = (expl - 1.0d0) / l
                end if
                arra%B_BN(k + h) = arra%B_BN(k + h) * divWH
            end do
        end do
    end subroutine update_AB

    ! Scale p arrays (in k-space) by 1/(W*H)
    subroutine scale_Ps(arra)
        type(Arrays), intent(inout) :: arra

        real(8), pointer :: P_scale_Ps_C(:), & 
P_scale_Ps_B(:), P_scale_Ps_N(:)
        real(8) :: divWH
        integer :: i, lA

        divWH = 1.0d0 / (arra%W * arra%H)
        P_scale_Ps_C => arra%p_arr_C
        P_scale_Ps_B => arra%p_arr_B
        P_scale_Ps_N => arra%p_arr_N
        lA = arra%lW * arra%H

        do i = 1, lA
            P_scale_Ps_C(i) = P_scale_Ps_C(i) * divWH
            P_scale_Ps_B(i) = P_scale_Ps_B(i) * divWH
            P_scale_Ps_N(i) = P_scale_Ps_N(i) * divWH
        end do
    end subroutine scale_Ps

    ! Scale q arrays (in k-space) by 1/(W*H)
    subroutine scale_Qs(arra)
        type(Arrays), intent(inout) :: arra

        real(8), pointer :: Q_C(:), Q_B(:), Q_N(:)
        real(8) :: divWH
        integer :: i, lA

        divWH = 1.0d0 / (arra%W * arra%H)
        Q_C => arra%q_arr_C
        Q_B => arra%q_arr_B
        Q_N => arra%q_arr_N
        lA = arra%lW * arra%H

        do i = 1, lA
            Q_C(i) = Q_C(i) * divWH
            Q_B(i) = Q_B(i) * divWH
            Q_N(i) = Q_N(i) * divWH
        end do
    end subroutine scale_Qs

    ! Compute average free energy densities and average densities
        subroutine fp(arra, mmode, relax)
        type(Arrays), intent(inout) :: arra
        type(Model), intent(in) :: mmode
        type(Relaxation), intent(inout) :: relax

        integer:: MPI_IN_PLACE, MPI_DOUBLE, MPI_MIN, MPI_MAX, MPI_SUM
        integer :: ierr, MPI_COMM_WORLD
 
        integer :: Wp, w, h, k, gw
        real(8) :: dkx, dky, kx2, ky, k2, d2
        real(8) :: q2, div3, a, divWH
        real(8) :: p_fp_C, p_fp_B, p_fp_N, & 
q_fp_C, q_fp_B, q_fp_N
        real(8) :: p_C2, p_B2, p_N2, v_C, v_BN, gamma_C, gamma_BN, dv, maskk
        real(8) :: alpha_C, beta_C, delta_C, alpha_C_B, gamma_C_B
        real(8) :: alpha_C_N, gamma_C_N, alpha_BN, beta_BN, delta_BN
        real(8) :: alpha_B_N, beta_B_N, gamma_B_N, alpha_B_C, gamma_B_C
        real(8) :: alpha_N_C, gamma_N_C
        real(8), pointer :: Qq_fp_C(:), Qq_fp_B(:), Qq_fp_N(:)

        Wp = 2 * (arra%W / 2 + 1)
        dkx = 2.0d0 * pi / (relax%dx * arra%W)
        dky = 2.0d0 * pi / (relax%dy * arra%H)
        div3 = 1.0d0 / 3.0d0
        a = -1.0d0 / (2.0d0 * mmode%sigma_mask * mmode%sigma_mask)
        q2 = 1.0d0 / (mmode%l_BN * mmode%l_BN)

        ! Save current state to p arrays
        arra%p_arr_C(:) = arra%q_arr_C(:)
        arra%p_arr_B(:) = arra%q_arr_B(:)
        arra%p_arr_N(:) = arra%q_arr_N(:)

        ! Execute FFTW plans
        call fftw_execute(arra%qa_Q_C)
        call fftw_execute(arra%qa_Q_B)
        call fftw_execute(arra%qa_Q_N)

        ! Scale Qs
        call scale_Qs(arra)

        Qq_fp_C => arra%q_arr_C
        Qq_fp_B => arra%q_arr_B
        Qq_fp_N => arra%q_arr_N

        ! Compute in k-space
        do w = 1, arra%lW
            gw = arra%lw0 + w - 1
            kx2 = gw * dkx
            kx2 = kx2 * kx2
            k = arra%H * (w - 1)
            do h = 1, arra%H
                if (h <= arra%H / 2) then
                    ky = (h - 1) * dky
                else
                    ky = (h - 1 - arra%H) * dky
                end if
                k2 = kx2 + ky * ky

                d2 = (1.0d0 - k2)
                d2 = d2 * d2
                Qq_fp_C(k + h) = Qq_fp_C(k + h) * d2

                d2 = (q2 - k2)
                d2 = d2 * d2
                Qq_fp_B(k + h) = Qq_fp_B(k + h) * d2
                Qq_fp_N(k + h) = Qq_fp_N(k + h) * d2
            end do
        end do

        ! Execute FFTW plans to restore q arrays
        call fftw_execute(arra%Q_q_C)
        call fftw_execute(arra%Q_q_B)
        call fftw_execute(arra%Q_q_N)

        ! Reset variables for average (free energy) density
        relax%f_C = 0.0d0
        relax%f_BN = 0.0d0
        relax%p_rex_C = 0.0d0
        relax%p_rex_BN = 0.0d0

        ! Model parameters
        alpha_C = mmode%alpha_C
        beta_C = mmode%beta_C
        delta_C = mmode%delta_C
        alpha_C_B = mmode%alpha_C_B
        gamma_C_B = mmode%gamma_C_B
        alpha_C_N = mmode%alpha_C_N
        gamma_C_N = mmode%gamma_C_N
        alpha_BN = mmode%alpha_BN
        beta_BN = mmode%beta_BN
        delta_BN = mmode%delta_BN
        alpha_B_N = mmode%alpha_B_N
        beta_B_N = mmode%beta_B_N
        gamma_B_N = mmode%gamma_B_N
        alpha_B_C = mmode%alpha_B_C
        gamma_B_C = mmode%gamma_B_C
        alpha_N_C = mmode%alpha_N_C
        gamma_N_C = mmode%gamma_N_C

        do h = 1, arra%lH
            k = Wp * (h - 1)
            do w = 1, arra%W
                p_fp_C = arra%p_arr_C(k + w)
                p_fp_B = arra%p_arr_B(k + w)
                p_fp_N = arra%p_arr_N(k + w)
                q_fp_C = arra%q_arr_C(k + w)
                q_fp_B = arra%q_arr_B(k + w)
                q_fp_N = arra%q_arr_N(k + w)
                p_C2 = p_fp_C * p_fp_C
                p_B2 = p_fp_B * p_fp_B
                p_N2 = p_fp_N * p_fp_N
                v_C = map(arra%u_arr_BN(k + w), arra%u_BN_max, & 
arra%u_BN_min)
                v_BN = map(arra%u_arr_C(k + w), arra%u_C_min, & 
arra%u_C_max)
                gamma_C = v_C * (mmode%gamma_C_s - mmode%gamma_C_l) + & 
mmode%gamma_C_l
                gamma_BN = v_BN * (mmode%gamma_BN_s - mmode%gamma_BN_l) + & 
mmode%gamma_BN_l
                dv = (v_C - v_BN) / (v_C + v_BN)
                maskk = exp(a * dv * dv)

                ! Average free energy density (graphene)
                relax%f_C = relax%f_C + &
       0.5d0 * alpha_C * p_C2 + 0.5d0 * beta_C * p_fp_C * q_fp_C + &
 div3 * gamma_C * p_fp_C * p_C2 + 0.25d0 * delta_C * p_C2 * p_C2 + &
                    maskk * (alpha_C_B * p_fp_C * p_fp_B + 0.5d0 * & 
gamma_C_B * (p_C2 * p_fp_B + p_fp_C * p_B2) & 
+ alpha_C_N * p_fp_C * p_fp_N + 0.5d0 * gamma_C_N * (p_C2 * p_fp_N + & 
p_fp_C * p_N2))

                ! Average free energy density (HBN)
                relax%f_BN = relax%f_BN + &
               0.5d0 * alpha_BN * p_B2 + 0.5d0 * beta_BN * p_fp_B * q_fp_B & 
+ div3 * gamma_BN * p_fp_B * p_B2 + 0.25d0 * delta_BN * p_B2 * p_B2 + &
                    0.5d0 * alpha_BN * p_N2 + 0.5d0 * beta_BN * p_fp_N * &
q_fp_N & 
+ div3 * gamma_BN * p_fp_N * p_N2 + 0.25d0 * delta_BN * p_N2 * p_N2 + &
                    alpha_B_N * p_fp_B * p_fp_N + beta_B_N * p_fp_B * &
q_fp_N & 
+ 0.5d0 * gamma_B_N * (p_B2 * p_fp_N + p_fp_B * p_N2) + &
                    maskk * (alpha_B_C * p_fp_B * p_fp_C + 0.5d0 * &
gamma_B_C * (p_B2 * p_fp_C + p_fp_B * p_C2) & 
+ alpha_N_C * p_fp_N * p_fp_C + 0.5d0 * gamma_N_C * (p_N2 * p_fp_C + &
p_fp_N * p_C2))

                ! Restore q arrays
                arra%q_arr_C(k + w) = p_fp_C
                arra%q_arr_B(k + w) = p_fp_B
                arra%q_arr_N(k + w) = p_fp_N

                ! Average densities
                relax%p_rex_C = relax%p_rex_C + p_fp_C
             relax%p_rex_BN = relax%p_rex_BN + 0.5d0 * (p_fp_B + p_fp_N)
            end do
        end do

        ! Communication between processes
        call MPI_Allreduce(MPI_IN_PLACE, relax%f_C, 1, MPI_DOUBLE, & 
MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, relax%f_BN, 1, MPI_DOUBLE, & 
MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, relax%p_rex_C, 1, MPI_DOUBLE, & 
MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, relax%p_rex_BN, 1, MPI_DOUBLE, & 
MPI_SUM, MPI_COMM_WORLD, ierr)

        divWH = 1.0d0 / (arra%W * arra%H)
        relax%f_C = relax%f_C * divWH
        relax%f_BN = relax%f_BN * divWH
        relax%p_rex_C = relax%p_rex_C * divWH
        relax%p_rex_BN = relax%p_rex_BN * divWH

        ! Execute FFTW plans to restore p arrays in k-space
        call fftw_execute(arra%p_P_C)
        call fftw_execute(arra%p_P_B)
        call fftw_execute(arra%p_P_N)

        ! Scale Ps
        call scale_Ps(arra)
    end subroutine fp

    ! Update smoothed density fields
    subroutine update_us(arra, mmode, relax)
        type(Arrays), intent(inout) :: arra
        type(Model), intent(in) :: mmode
        type(Relaxation), intent(in) :: relax

        integer :: MPI_DOUBLE, MPI_MAX, MPI_MIN, MPI_COMM_WORLD  

        integer :: Wp, w, gw, h, k
        integer :: MPI_IN_PLACE
        real(8) :: divWH, ky, kx2, k2, a, dkx, dky
        real(8), pointer :: Uu_C(:), Uu_BN(:)

        Wp = 2 * (arra%W / 2 + 1)
        divWH = 1.0d0 / (arra%W * arra%H)
        a = -1.0d0 / (2.0d0 * mmode%sigma_u * mmode%sigma_u)
        dkx = 2.0d0 * pi / (relax%dx * arra%W)
        dky = 2.0d0 * pi / (relax%dy * arra%H)

        do h = 1, arra%lH
            k = Wp * (h - 1)
            do w = 1, arra%W
                arra%u_arr_C(k + w) = arra%q_arr_C(k + w)
                arra%u_arr_BN(k + w) = 0.5d0 * (arra%q_arr_B(k + w) + & 
arra%q_arr_N(k + w))
            end do
        end do

        call fftw_execute(arra%u_U_C)
        call fftw_execute(arra%u_U_BN)

        Uu_C => arra%u_arr_C
        Uu_BN => arra%u_arr_BN

        do w = 1, arra%lW
            gw = arra%lw0 + w - 1
            kx2 = gw * dkx
            kx2 = kx2 * kx2
            k = arra%H * (w - 1)
            do h = 1, arra%H
                if (h <= arra%H / 2) then
                    ky = (h - 1) * dky
                else
                    ky = (h - 1 - arra%H) * dky
                end if
                k2 = kx2 + ky * ky
                Uu_C(k + h) = Uu_C(k + h) * exp(a * k2) * divWH
                Uu_BN(k + h) = Uu_BN(k + h) * exp(a * k2) * divWH
            end do
        end do

        call fftw_execute(arra%U_u_C)
        call fftw_execute(arra%U_u_BN)

        ! Reset variables for density field extrema
        arra%u_C_min = 1.0d100
        arra%u_C_max = -1.0d100
        arra%u_BN_min = 1.0d100
        arra%u_BN_max = -1.0d100

        do h = 1, arra%lH
            k = Wp * (h - 1)
            do w = 1, arra%W
                if (arra%u_arr_C(k + w) < arra%u_C_min) then 
                   arra%u_C_min = arra%u_arr_C(k + w)
                if (arra%u_arr_C(k + w) > arra%u_C_max) then 
                   arra%u_C_max = arra%u_arr_C(k + w)
                if (arra%u_arr_BN(k + w) < arra%u_BN_min) then 
                   arra%u_BN_min = arra%u_arr_BN(k + w)
                if (arra%u_arr_BN(k + w) > arra%u_BN_max) then 
                   arra%u_BN_max = arra%u_arr_BN(k + w)
                end if
                 end if
                 end if
                 end if  
            end do
        end do

        call MPI_Allreduce(MPI_IN_PLACE, arra%u_C_min, 1, & 
MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD)
        call MPI_Allreduce(MPI_IN_PLACE, arra%u_C_max, 1, & 
MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)
        call MPI_Allreduce(MPI_IN_PLACE, arra%u_BN_min, 1, & 
MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD)
        call MPI_Allreduce(MPI_IN_PLACE, arra%u_BN_max, 1, & 
MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)
    end subroutine update_us

    ! Perform one iteration of the semi-implicit spectral method
    subroutine step(arra, mmode, relax)
        type(Arrays), intent(inout) :: arra
        type(Model), intent(in) :: mmode
        type(Relaxation), intent(in) :: relax

        integer :: Wp, w, h, k, gw
        real(8) :: v_C, v_BN, gamma_C, gamma_BN, dv, mask, a
        real(8) :: q_step_C, q_step_B, q_step_N, q_C2, q_B2, q_N2
        real(8) :: delta_C, alpha_C_B, gamma_C_B, alpha_C_N, gamma_C_N
        real(8) :: delta_BN, alpha_B_N, beta_B_N, gamma_B_N, alpha_B_C, & 
gamma_B_C
        real(8) :: alpha_N_C, gamma_N_C
        real(8) :: WH, dkx, dky, kx2, ky, k2, d2, q2
        complex(8) :: P_B_prev
        real(8), pointer :: P_C(:), P_B(:), P_N(:), Q_C(:), Q_B(:), Q_N(:)

        Wp = 2 * (arra%W / 2 + 1)
        a = -1.0d0 / (2.0d0 * mmode%sigma_mask * mmode%sigma_mask)
        WH = arra%W * arra%H
        dkx = 2.0d0 * pi / (relax%dx * arra%W)
        dky = 2.0d0 * pi / (relax%dy * arra%H)
        q2 = 1.0d0 / (mmode%l_BN * mmode%l_BN)

        ! Avoid numerical instability
        if (mod(relax%tt, 100) == 0) then
            arra%p_arr_C(:) = arra%q_arr_C(:)
            arra%p_arr_B(:) = arra%q_arr_B(:)
            arra%p_arr_N(:) = arra%q_arr_N(:)
            call fftw_execute(arra%p_P_C)
            call fftw_execute(arra%p_P_B)
            call fftw_execute(arra%p_P_N)
            call scale_Ps(arra)
        end if

        ! Nonlinear part
        delta_C = mmode%delta_C
        alpha_C_B = mmode%alpha_C_B
        gamma_C_B = mmode%gamma_C_B
        alpha_C_N = mmode%alpha_C_N
        gamma_C_N = mmode%gamma_C_N
        delta_BN = mmode%delta_BN
        alpha_B_N = mmode%alpha_B_N
        beta_B_N = mmode%beta_B_N
        gamma_B_N = mmode%gamma_B_N
        alpha_B_C = mmode%alpha_B_C
        gamma_B_C = mmode%gamma_B_C
        alpha_N_C = mmode%alpha_N_C
        gamma_N_C = mmode%gamma_N_C

        do h = 1, arra%lH
            k = Wp * (h - 1)
            do w = 1, arra%W
                v_C = map(arra%u_arr_BN(k + w), arra%u_BN_max, & 
arra%u_BN_min)
                v_BN = map(arra%u_arr_C(k + w), arra%u_C_min, & 
arra%u_C_max)
                gamma_C = v_C * (mmode%gamma_C_s - mmode%gamma_C_l) + & 
mmode%gamma_C_l
                gamma_BN = v_BN * (mmode%gamma_BN_s - mmode%gamma_BN_l) + &  
mmode%gamma_BN_l
                dv = (v_C - v_BN) / (v_C + v_BN)
                mask = exp(a * dv * dv)
                q_step_C = arra%q_arr_C(k + w)
                q_step_B = arra%q_arr_B(k + w)
                q_step_N = arra%q_arr_N(k + w)
                q_C2 = q_step_C * q_step_C
                q_B2 = q_step_B * q_step_B
                q_N2 = q_step_N * q_step_N
                ! compute nonlinear part (C)
                arra%q_arr_C(k + w) = gamma_C * q_C2 + delta_C * q_step_C & 
* q_C2 + &
                    mask * (alpha_C_B * q_step_B + gamma_C_B * (q_step_C * & 
q_step_B + 0.5d0 * q_B2) & 
+ alpha_C_N * q_step_N + gamma_C_N * (q_step_C * q_step_N + 0.5d0 * q_N2))
                ! (B)
                arra%q_arr_B(k + w) = gamma_BN * q_B2 + delta_BN * & 
q_step_B * q_B2 & 
+ alpha_B_N * q_step_N + gamma_B_N * (q_step_B * q_step_N + 0.5d0 * q_N2) & 
+ &
                    mask * (alpha_B_C * q_step_C + gamma_B_C * (q_step_B * & 
q_step_C + 0.5d0 * q_C2))
                ! (N)
                arra%q_arr_N(k + w) = gamma_BN * q_N2 + delta_BN * & 
q_step_N * q_N2 & 
+ alpha_B_N * q_step_B + gamma_B_N * (q_step_N * q_step_B + 0.5d0 * q_B2) & 
+ &
                    mask * (alpha_N_C * q_step_C + gamma_N_C * (q_step_N * & 
q_step_C + 0.5d0 * q_C2))
            end do
        end do

        call fftw_execute(arra%qa_Q_C)
        call fftw_execute(arra%qa_Q_B)
        call fftw_execute(arra%qa_Q_N)

        P_C => arra%p_arr_C
        P_B => arra%p_arr_B
        P_N => arra%p_arr_N
        Q_C => arra%q_arr_C
        Q_B => arra%q_arr_B
        Q_N => arra%q_arr_N

        ! Update \psi_C, \psi_B, and \psi_N in k-space
        do w = 1, arra%lW
            gw = arra%lw0 + w - 1
            kx2 = gw * dkx
            kx2 = kx2 * kx2
            k = arra%H * (w - 1)
            do h = 1, arra%H
                P_B_prev = P_B(k + h)    
! P_B(k + h) will be overwritten 
! so value saved for computing P_N(k + h)
                if (h <= arra%H / 2) then
                    ky = (h - 1) * dky
                else
                    ky = (h - 1 - arra%H) * dky
                end if
                k2 = kx2 + ky * ky
                ! Update \psi_C in k-space
                P_C(k + h) = arra%A_C(k + h) * P_C(k + h) + arra%B_C(k + & 
h) * Q_C(k + h)
                ! \psi_B and \psi_N
                d2 = q2 - k2
                d2 = d2 * d2
                P_B(k + h) = arra%A_BN(k + h) * P_B(k + h) & 
+ arra%B_BN(k + h) * (Q_B(k + h) + beta_B_N * d2 * P_N(k + h) * WH)
                P_N(k + h) = arra%A_BN(k + h) * P_N(k + h) & 
+ arra%B_BN(k + h) * (Q_N(k + h) + beta_B_N * d2 * P_B_prev * WH)
                ! Copy p's to q's
                Q_C(k + h) = P_C(k + h)
                Q_B(k + h) = P_B(k + h)
                Q_N(k + h) = P_N(k + h)
            end do
        end do

        call fftw_execute(arra%Q_q_C)
        call fftw_execute(arra%Q_q_B)
        call fftw_execute(arra%Q_q_N)

        ! Update smoothed density fields
        call update_us(arra, mmode, relax)
    end subroutine step

    ! Optimize calculation box size with system size
    subroutine optimize(arra, mmode, relax)
        type(Arrays), intent(inout) :: arra
        type(Model), intent(in) :: mmode
        type(Relaxation), intent(inout) :: relax

        real(8) :: dxs(5), dys(5), fs(5)
        real(8) :: l0, dw, dh, dr, x, ddx, ddy, ddr
        integer :: i

        dxs = [relax%dx, relax%dx-relax%d, relax%dx+relax%d, relax%dx, & 
relax%dx]
        dys = [relax%dy, relax%dy, relax%dy, relax%dy-relax%d, & 
relax%dy+relax%d]

        do i = 1, 5
            relax%dx = dxs(i)
            relax%dy = dys(i)
            call fp(arra, mmode, relax)
            fs(i) = relax%f_C + relax%f_BN
        end do

        relax%dx = (relax%d * (fs(2) - fs(3)) + 2.0d0 * dxs(1) * (-2.0d0 * fs(1) & 
+ fs(2) + fs(3))) / (2.0d0 * (fs(2) - 2.0d0 * fs(1) + fs(3)))
        relax%dy = (relax%d * (fs(4) - fs(5)) + 2.0d0 * dys(1) * (-2.0d0 * fs(1) & 
+ fs(4) + fs(5))) / (2.0d0 * (fs(4) - 2.0d0 * fs(1) + fs(5)))

        l0 = 4.0d0 * pi / sqrt(3.0d0)
        dw = arra%W * (relax%dx - dxs(1))
        dh = arra%H * (relax%dy - dys(1))
        dr = sqrt(dw * dw + dh * dh)
        x = 0.25d0 * l0

        if (dr > x) then
            x = x / dr
            relax%dx = x * relax%dx + (1.0d0 - x) * dxs(1)
            relax%dy = x * relax%dy + (1.0d0 - x) * dys(1)
        end if

        call update_AB(arra, mmode, relax)

        ddx = relax%dx - dxs(1)
        ddy = relax%dy - dys(1)
        ddr = sqrt(ddx * ddx + ddy * ddy)

        if (ddr < relax%d) then
            relax%d = relax%d * 0.5d0
        else
            relax%d = relax%d * 2.0d0
        end if

        if (relax%d < 1.0d-6) relax%d = relax%d * 2.0d0
    end subroutine optimize

    ! Relax the system for T time steps
    subroutine relaxx(arra, mmode, relax, outp)
        type(Arrays), intent(inout) :: arra
        type(Model), intent(in) :: mmode
        type(Relaxation), intent(inout) :: relax
        type(Output), intent(in) :: outp
        integer :: ii 

        do ii = 0, relax%T
             relax%tt = ii 
            if (relax%T_optimize > 0 .and. relax%tt > 0 .and. & 
mod(relax%tt, relax%T_optimize) == 0) then  
            call optimize(arra, mmode, relax)
            if (mod(relax%tt, outp%T_print) == 0) then
                call fp(arra, mmode, relax)
                call print_output(relax, outp)
            end if
            if (mod(relax%tt, outp%T_write) == 0) then 
               call write_state(arra, relax, outp)
               end if 
            if (relax%tt < relax%T) then 
               call step(arra, mmode, relax)
               end if 
             end if
        end do
    end subroutine relaxx

    ! Free allocated arrays
    subroutine clear_arrays(arra)
        type(Arrays), intent(inout) :: arra
 
        call fftw_free(arra%p_arr_C_cptr)
        call fftw_free(arra%q_arr_C_cptr)
        call fftw_free(arra%p_arr_B_cptr)
        call fftw_free(arra%q_arr_B_cptr)
        call fftw_free(arra%p_arr_N_cptr)
        call fftw_free(arra%q_arr_N_cptr)
        call fftw_free(arra%u_arr_C_cptr)
        call fftw_free(arra%u_arr_BN_cptr)
        call fftw_free(arra%A_C_cptr)
        call fftw_free(arra%B_C_cptr)
        call fftw_free(arra%A_BN_cptr)
        call fftw_free(arra%B_BN_cptr)
    end subroutine clear_arrays

end module utilities

program heterostructure
  use, intrinsic :: iso_c_binding
 ! use mpi
  use omp_lib
   
    use fftwa3_mpi
    use constants
    use data_types
    use utilities
!    use mpi_f08

    implicit none
!    include 'aslfftw3-mpi.f03'
!   include 'mpi.h'
!   include 'fftw3.f03'

    type(Arrays) :: arra
    type(Model) :: mmode
    type(Relaxation) :: relax
    type(Output) :: outp

    integer :: ierr, input, argc
    character(len=128) :: filename

    call MPI_Init(ierr)
    call fftw_mpi_init()

    call MPI_Comm_rank(MPI_COMM_WORLD, relax%iid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, relax%ID, ierr)
    relax%t0 = time()
    relax%tt = 0
    relax%d = 0.0001d0

    call get_command_argument(1, outp%name)

    write(filename, '(A, ".in")') trim(outp%name)
    open(unit=input, file=filename, status='old', action='read', 
iostat=ierr)
    if (ierr /= 0) then
        print *, 'Error: Input file not found!'
        call MPI_Abort(MPI_COMM_WORLD, ierr)
    end if

    write(filename, '(A, ".out")') trim(outp%name)
    open(unit=input, file=filename, status='unknown', action='write')
    close(input)

    do
        read(input, '(A)', advance='no') filename
        select case (filename)
            case ('S')
                call seed_rngs(relax, input)
            case ('O')
                call configure_output(outp, input)
            case ('A')
                call configure_arrays(arra, input)
            case ('I')
                call initialize_system(arra, input)
            case ('M')
                call configure_model(mmode, input)
            case ('R')
                call configure_relaxation(relax, input)
                call update_AB(arra, mmode, relax)
                call update_us(arra, mmode, relax)
                call relax(arra, mmode, relax, outp)
            case default
                print *, 'Invalid label: ', filename
                call MPI_Abort(MPI_COMM_WORLD, ierr)
        end select
    end do

    close(input)

    call clear_arrays(arra)
    call MPI_Finalize(ierr)

end program heterostructure
