module T
    implicit none

    integer, parameter :: ip_ = 4
    integer, parameter :: ipc_ = 4
    integer, parameter :: rp_ = 8
    integer, parameter :: rpc_ = 8
    REAL ( KIND = rp_ ), PARAMETER :: zero = 0.0_rp_
    REAL ( KIND = rp_ ), PARAMETER :: one = 1.0_rp_
    REAL ( KIND = rp_ ), PARAMETER :: ten = 10.0_rp_

    TYPE :: LANCELOT_problem_type
        INTEGER ( KIND = ip_ ) :: n, ng, nel
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IELING
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISTADG
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IELVAR
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISTAEV
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: INTVAR
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISTADH
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ICNA
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISTADA
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: KNDOFG
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ITYPEE
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISTEPA
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ITYPEG
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISTGPA
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: A
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: B
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: BL
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: BU
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: X
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: C
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: Y
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: GSCALE
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: ESCALE
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: VSCALE
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: EPVALU
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: GPVALU
        LOGICAL, ALLOCATABLE, DIMENSION( : ) :: GXEQX
        LOGICAL, ALLOCATABLE, DIMENSION( : ) :: INTREP
        CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: VNAMES
        CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: GNAMES
    END TYPE LANCELOT_problem_type

    TYPE, BIND( C ) :: scu_inform_type
        INTEGER ( KIND = ipc_ ) :: status
        INTEGER ( KIND = ipc_ ) :: alloc_status
        INTEGER ( KIND = ipc_ ), DIMENSION( 3 ) :: inertia
    END TYPE scu_inform_type

    TYPE :: SILS_ainfo
        INTEGER ( KIND = ip_ ) :: flag = 0
        INTEGER ( KIND = ip_ ) :: more = 0
        INTEGER ( KIND = ip_ ) :: nsteps = 0
        INTEGER ( KIND = ip_ ) :: nrltot = - 1
        INTEGER ( KIND = ip_ ) :: nirtot = - 1
        INTEGER ( KIND = ip_ ) :: nrlnec = - 1
        INTEGER ( KIND = ip_ ) :: nirnec = - 1
        INTEGER ( KIND = ip_ ) :: nrladu = - 1
        INTEGER ( KIND = ip_ ) :: niradu = - 1
        INTEGER ( KIND = ip_ ) :: ncmpa  = 0
        INTEGER ( KIND = ip_ ) :: oor = 0
        INTEGER ( KIND = ip_ ) :: dup = 0
        INTEGER ( KIND = ip_ ) :: maxfrt = - 1
        INTEGER ( KIND = ip_ ) :: stat = 0
        INTEGER ( KIND = ip_ ) :: faulty = 0
        REAL ( KIND = rp_ ) :: opsa = - 1.0_rp_
        REAL ( KIND = rp_ ) :: opse = - 1.0_rp_
    END TYPE SILS_ainfo

    TYPE :: SILS_finfo
        INTEGER ( KIND = ip_ ) :: flag = 0
        INTEGER ( KIND = ip_ ) :: more = 0
        INTEGER ( KIND = ip_ ) :: maxfrt = - 1
        INTEGER ( KIND = ip_ ) :: nebdu  = - 1
        INTEGER ( KIND = ip_ ) :: nrlbdu = - 1
        INTEGER ( KIND = ip_ ) :: nirbdu = - 1
        INTEGER ( KIND = ip_ ) :: nrltot = - 1
        INTEGER ( KIND = ip_ ) :: nirtot = - 1
        INTEGER ( KIND = ip_ ) :: nrlnec = - 1
        INTEGER ( KIND = ip_ ) :: nirnec = - 1
        INTEGER ( KIND = ip_ ) :: ncmpbr = - 1
        INTEGER ( KIND = ip_ ) :: ncmpbi = - 1
        INTEGER ( KIND = ip_ ) :: ntwo = - 1
        INTEGER ( KIND = ip_ ) :: neig = - 1
        INTEGER ( KIND = ip_ ) :: delay = - 1
        INTEGER ( KIND = ip_ ) :: signc = - 1
        INTEGER ( KIND = ip_ ) :: static = - 1
        INTEGER ( KIND = ip_ ) :: modstep = -1
        INTEGER ( KIND = ip_ ) :: rank = - 1
        INTEGER ( KIND = ip_ ) :: stat = - 1
        INTEGER ( KIND = ip_ ) :: faulty = - 1
        INTEGER ( KIND = ip_ ) :: step = - 1
        REAL ( KIND = rp_ ) :: opsa = - 1.0_rp_
        REAL ( KIND = rp_ ) :: opse = - 1.0_rp_
        REAL ( KIND = rp_ ) :: opsb = - 1.0_rp_
        REAL ( KIND = rp_ ) :: maxchange = - 1.0_rp_
        REAL ( KIND = rp_ ) :: smin = 0.0_rp_
        REAL ( KIND = rp_ ) :: smax = 0.0_rp_
    END TYPE SILS_finfo

    TYPE :: SILS_sinfo
        INTEGER ( KIND = ip_ ) :: flag = - 1
        INTEGER ( KIND = ip_ ) :: stat = - 1
        REAL ( KIND = rp_ ) :: cond = - 1.0_rp_
        REAL ( KIND = rp_ ) :: cond2 = - 1.0_rp_
        REAL ( KIND = rp_ ) :: berr = - 1.0_rp_
        REAL ( KIND = rp_ ) :: berr2 = - 1.0_rp_
        REAL ( KIND = rp_ ) :: error = - 1.0_rp_
    END TYPE SILS_sinfo

    TYPE :: LANCELOT_inform_type
        INTEGER ( KIND = ip_ ) :: status = 7777
        INTEGER ( KIND = ip_ ) :: alloc_status = 0
        INTEGER ( KIND = ip_ ) :: iter = - 1
        INTEGER ( KIND = ip_ ) :: itercg = - 1
        INTEGER ( KIND = ip_ ) :: itcgmx = - 1
        INTEGER ( KIND = ip_ ) :: ncalcf = 0
        INTEGER ( KIND = ip_ ) :: ncalcg = 0
        INTEGER ( KIND = ip_ ) :: nvar = 0
        INTEGER ( KIND = ip_ ) :: ngeval = 0
        INTEGER ( KIND = ip_ ) :: iskip = 0
        INTEGER ( KIND = ip_ ) :: ifixed = 0
        INTEGER ( KIND = ip_ ) :: nsemib = 0
        REAL ( KIND = rp_ ) :: aug = HUGE( one )
        REAL ( KIND = rp_ ) :: obj = HUGE( one )
        REAL ( KIND = rp_ ) :: pjgnrm = HUGE( one )
        REAL ( KIND = rp_ ) :: cnorm = zero
        REAL ( KIND = rp_ ) :: ratio = zero
        REAL ( KIND = rp_ ) :: mu = zero
        REAL ( KIND = rp_ ) :: radius = zero
        REAL ( KIND = rp_ ) :: ciccg = zero
        LOGICAL :: newsol = .FALSE.
        CHARACTER ( LEN = 80 ) :: bad_alloc =  REPEAT( ' ', 80 )
        TYPE ( SCU_inform_type ) :: SCU_info
        TYPE ( SILS_ainfo ) :: SILS_infoa
        TYPE ( SILS_finfo ) :: SILS_infof
        TYPE ( SILS_sinfo ) :: SILS_infos
        INTEGER ( KIND = ip_ ) :: ende = 6666
    END TYPE LANCELOT_inform_type

    TYPE :: SILS_control
        INTEGER ( KIND = ip_ ) :: ICNTL( 30 ) =                                 &
            (/ 6, 6, 0, 2139062143, 1, 32639, 32639, 32639, 32639, 14,           &
                9, 8, 8, 9, 10, 32639, 32639, 32639, 32689, 24,                   &
                11, 9, 8, 9, 10, 0, 0, 0, 0, 0 /)
        INTEGER ( KIND = ip_ ) :: lp = 6
        INTEGER ( KIND = ip_ ) :: wp = 6
        INTEGER ( KIND = ip_ ) :: mp = 6
        INTEGER ( KIND = ip_ ) :: sp = - 1
        INTEGER ( KIND = ip_ ) :: ldiag = 0
        INTEGER ( KIND = ip_ ) :: la = 0
        INTEGER ( KIND = ip_ ) :: liw = 0
        INTEGER ( KIND = ip_ ) :: maxla = HUGE( 0 )
        INTEGER ( KIND = ip_ ) :: maxliw = HUGE( 0 )
        INTEGER ( KIND = ip_ ) :: pivoting = 1
        INTEGER ( KIND = ip_ ) :: nemin = 1
        INTEGER ( KIND = ip_ ) :: factorblocking = 16
        INTEGER ( KIND = ip_ ) :: solveblocking = 16
        INTEGER ( KIND = ip_ ) :: thresh = 50
        INTEGER ( KIND = ip_ ) :: ordering = 3
        INTEGER ( KIND = ip_ ) :: scaling = 0
        REAL ( KIND = rp_ ) :: CNTL( 5 ) =                                      &
            (/ 0.1_rp_,  1.0_rp_,  0.0_rp_,  0.0_rp_,  0.0_rp_ /)
        REAL ( KIND = rp_ ) :: multiplier  = 2.0_rp_
        REAL ( KIND = rp_ ) :: reduce  = 2.0_rp_
        REAL ( KIND = rp_ ) :: u = 0.1_rp_
        REAL ( KIND = rp_ ) :: static_tolerance = 0.0_rp_
        REAL ( KIND = rp_ ) :: static_level = 0.0_rp_
        REAL ( KIND = rp_ ) :: tolerance  = 0.0_rp_
        REAL ( KIND = rp_ ) :: convergence = 0.5_rp_
    END TYPE SILS_control

    TYPE :: LANCELOT_control_type
        INTEGER ( KIND = ip_ ) :: error = 6
        INTEGER ( KIND = ip_ ) :: out = 6
        INTEGER ( KIND = ip_ ) :: alive_unit = 60
        CHARACTER ( LEN = 30 ) :: alive_file = 'ALIVE.d                       '
        INTEGER ( KIND = ip_ ) :: print_level = 0
        INTEGER ( KIND = ip_ ) :: maxit = 1000
        INTEGER ( KIND = ip_ ) :: start_print = - 1
        INTEGER ( KIND = ip_ ) :: stop_print = - 1
        INTEGER ( KIND = ip_ ) :: print_gap = 1
        INTEGER ( KIND = ip_ ) :: linear_solver = 8
        INTEGER ( KIND = ip_ ) :: icfact = 5
        INTEGER ( KIND = ip_ ) :: semibandwidth = 5
        INTEGER ( KIND = ip_ ) :: max_sc = 100
        INTEGER ( KIND = ip_ ) :: io_buffer = 75
        INTEGER ( KIND = ip_ ) :: more_toraldo = 0
        INTEGER ( KIND = ip_ ) :: non_monotone = 1
        INTEGER ( KIND = ip_ ) :: first_derivatives = 0
        INTEGER ( KIND = ip_ ) :: second_derivatives = 0
        REAL ( KIND = rp_ ) :: stopc = ten ** ( - 5 )
        REAL ( KIND = rp_ ) :: stopg = ten ** ( - 5 )
        REAL ( KIND = rp_ ) :: min_aug = - ( HUGE( one ) / 8.0_rp_ )
        REAL ( KIND = rp_ ) :: acccg = 0.01_rp_
        REAL ( KIND = rp_ ) :: initial_radius = - one
        REAL ( KIND = rp_ ) :: maximum_radius = ten ** 20
        REAL ( KIND = rp_ ) :: eta_successful = 0.01_rp_
        REAL ( KIND = rp_ ) :: eta_very_successful = 0.9_rp_
        REAL ( KIND = rp_ ) :: eta_extremely_successful = 0.95_rp_
        REAL ( KIND = rp_ ) :: gamma_smallest = 0.0625_rp_
        REAL ( KIND = rp_ ) :: gamma_decrease = 0.25_rp_
        REAL ( KIND = rp_ ) :: gamma_increase = 2.0_rp_
        REAL ( KIND = rp_ ) :: mu_meaningful_model = 0.01_rp_
        REAL ( KIND = rp_ ) :: mu_meaningful_group = 0.1_rp_
        REAL ( KIND = rp_ ) :: initial_mu = 0.1_rp_
        REAL ( KIND = rp_ ) :: mu_decrease = 0.1_rp_
        REAL ( KIND = rp_ ) :: mu_steering_decrease = 0.7_rp_
        REAL ( KIND = rp_ ) :: mu_tol = 0.1_rp_
        REAL ( KIND = rp_ ) :: firstg = 0.1_rp_
        REAL ( KIND = rp_ ) :: firstc = 0.1_rp_
        INTEGER ( KIND = ip_ ) :: num_mudec = HUGE( 1 )
        INTEGER ( KIND = ip_ ) :: num_mudec_per_iteration = HUGE( 1 )
        REAL ( KIND = rp_ ) :: kappa_3 = ten ** ( - 5 )
        REAL ( KIND = rp_ ) :: kappa_t = 0.9_rp_
        REAL ( KIND = rp_ ) :: mu_min = zero
        REAL ( KIND = rp_ ) :: cpu_time_limit = - one
        LOGICAL :: quadratic_problem = .FALSE.
        LOGICAL :: steering = .FALSE.
        LOGICAL :: two_norm_tr = .FALSE.
        LOGICAL :: exact_gcp = .TRUE.
        LOGICAL :: gn_model = .FALSE.
        LOGICAL :: gn_model_after_cauchy = .FALSE.
        LOGICAL :: magical_steps = .FALSE.
        LOGICAL :: accurate_bqp = .FALSE.
        LOGICAL :: structured_tr = .FALSE.
        LOGICAL :: print_max = .FALSE.
        LOGICAL :: full_solution = .TRUE.
        LOGICAL :: space_critical = .FALSE.
        LOGICAL :: deallocate_error_fatal = .FALSE.
        TYPE ( SILS_control ) :: SILS_cntl
    END TYPE LANCELOT_control_type

    TYPE :: EXTEND_save_type
        INTEGER ( KIND = ip_ ) :: lirnh, ljcnh, llink_min, lirnh_min
        INTEGER ( KIND = ip_ ) :: ljcnh_min, lh_min, lh, litran_min
        INTEGER ( KIND = ip_ ) :: lwtran_min, lwtran, litran, l_link_e_u_v
        INTEGER ( KIND = ip_ ) :: llink, lrowst, lpos, lused, lfilled
    END TYPE EXTEND_save_type

    TYPE :: CAUCHY_save_type
        INTEGER ( KIND = ip_ ) :: iterca, iter, itmax, nfreed, nbreak, nzero
        REAL ( KIND = rp_ ) :: tk, gxt, hxt, epstl2, tpttp, tcauch
        REAL ( KIND = rp_ ) :: tbreak, deltat, epsqrt, gxtold, g0tp
        REAL ( KIND = rp_ ) :: t, tamax , ptp, gtp, flxt, tnew
        LOGICAL :: prnter, pronel, recomp
    END TYPE CAUCHY_save_type

    TYPE :: CG_save_type
        INTEGER ( KIND = ip_ ) :: iter, itsle
        REAL ( KIND = rp_ ) :: alpha, oldgns, onepep
        LOGICAL :: prnter, pronel
    END TYPE CG_save_type

    TYPE :: ASMBL_save_type
        LOGICAL :: ptr_status
        INTEGER ( KIND = ip_ ), DIMENSION( 30 ) :: ICNTL
        INTEGER ( KIND = ip_ ), DIMENSION( 20 ) :: INFO
        REAL ( KIND = rp_ ), DIMENSION( 5 ) :: CNTL
    END TYPE ASMBL_save_type

    TYPE :: PRECN_save_type
        INTEGER ( KIND = ip_ ) :: liw, lw, nsemiw, nupdat, liccgg
        INTEGER ( KIND = ip_ ) :: nextra, nz01, iaj
        REAL :: tfactr, t1stsl, tupdat, tsolve
        INTEGER ( KIND = ip_ ) :: ICNTL_iccg( 5 ), KEEP_iccg( 12 )
        INTEGER ( KIND = ip_ ) :: INFO_iccg( 10 )
        REAL ( KIND = rp_ ) :: CNTL_iccg( 3 )
    END TYPE PRECN_save_type

    TYPE :: OTHERS_fdgrad_save_type
        LOGICAL :: backwd
    END TYPE OTHERS_fdgrad_save_type

    TYPE :: SCU_matrix_type
        INTEGER ( KIND = ip_ ) :: n, m, m_max
        INTEGER ( KIND = ip_ ) :: class = 0
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: BD_row
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: BD_col_start
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: CD_col
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: CD_row_start
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: BD_val
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: CD_val
    END TYPE SCU_matrix_type

    TYPE :: SCU_data_type
        INTEGER ( KIND = ip_ ) :: m, m_max, jumpto, jcol, newdia, sign_determinant
        INTEGER ( KIND = ip_ ) :: class = 3
        LOGICAL :: got_factors
        REAL ( KIND = rp_ ) :: dianew
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: R
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: W
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : , : ) :: Q
    END TYPE SCU_data_type

    TYPE :: SMT_type
        INTEGER ( KIND = ip_ ) :: m, n, ne
        CHARACTER, ALLOCATABLE, DIMENSION(:) :: id, type
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION(:) :: row, col, ptr
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION(:) :: val
    END TYPE

    TYPE :: SILS_factors
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( :, : )  :: keep
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : )  :: iw
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : )  :: iw1
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( :, : )  :: iw2
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : )  :: val
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : )  :: w
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : )  :: r
        INTEGER ( KIND = ip_ ) :: n = - 1
        INTEGER ( KIND = ip_ ) :: nrltot = - 1
        INTEGER ( KIND = ip_ ) :: nirtot = - 1
        INTEGER ( KIND = ip_ ) :: nrlnec = - 1
        INTEGER ( KIND = ip_ ) :: nirnec = - 1
        INTEGER ( KIND = ip_ ) :: nsteps = - 1
        INTEGER ( KIND = ip_ ) :: maxfrt = - 1
        INTEGER ( KIND = ip_ ) :: latop = - 1
        INTEGER ( KIND = ip_ ) :: dim_iw1 = - 1
        INTEGER ( KIND = ip_ ) :: pivoting = - 1
        REAL ( KIND = rp_ ) :: ops = - one
    END TYPE SILS_factors

    TYPE :: LANCELOT_save_type
        LOGICAL :: full_solution
        INTEGER ( KIND = ip_ ) :: igetfd
        LOGICAL :: unsucc
        INTEGER ( KIND = ip_ ) :: nobjgr, m, icrit, ncrit, p_type
        REAL ( KIND = rp_ ) :: ocnorm, cnorm_major, etak, eta0, omegak, omega0
        REAL ( KIND = rp_ ) :: tau, tau_steering, gamma1, alphae, betae, alphak
        REAL ( KIND = rp_ ) :: alphao, betao, omega_min, eta_min, epstol, epsgrd
        REAL ( KIND = rp_ ) :: cnorm
        CHARACTER ( LEN = 5 ), DIMENSION( 5 ) :: STATE
        LOGICAL :: itzero, reeval
        INTEGER ( KIND = ip_ ) :: ifactr, ldx, lfxi, lgxi, lhxi, lggfx, nvar2
        INTEGER ( KIND = ip_ ) :: nfreef, nfree, nnonnz, nadd, icfact
        INTEGER ( KIND = ip_ ) :: jumpto, nbprod, infor, number, nfixed, ibqpst
        INTEGER ( KIND = ip_ ) :: nmhist, maxsel, ntotin, nfreec, lnguvl, lnhuvl
        INTEGER ( KIND = ip_ ) :: ntype, nsets, nvargp, l_suc, msweep, nbnd
        INTEGER ( KIND = ip_ ) :: mortor_its, ntotel,inform_status
        INTEGER ( KIND = ip_ ) :: nvrels, nnza, error, out, print_level
        INTEGER ( KIND = ip_ ) :: start_print, stop_print, print_gap
        INTEGER ( KIND = ip_ ) :: n_steering, n_steering_this_iteration
        INTEGER ( KIND = ip_ ) :: first_derivatives, second_derivatives
        REAL ( KIND = rp_ ) :: epstlp, gmodel, vscmax, rad, maximum_radius
        REAL ( KIND = rp_ ) :: epsrcg, fnew, radmin, cgstop, diamin, diamax
        REAL ( KIND = rp_ ) :: ared, prered, rho, fmodel, curv, dxsqr, fcp, f0
        REAL ( KIND = rp_ ) :: stepmx, smallh, resmin, qgnorm, oldrad, epscns
        REAL ( KIND = rp_ ) :: radtol, fill, step, teneps, stpmin, epstln
        REAL ( KIND = rp_ ) :: f_min, f_r, f_c, sigma_r, sigma_c, findmx
        REAL ( KIND = rp_ ) :: f_min_lag, f_r_lag, f_c_lag
        REAL ( KIND = rp_ ) :: f_min_viol, f_r_viol, f_c_viol
        REAL ( KIND = rp_ ) :: violation, delta_qv, delta_qv_steering
        LOGICAL :: alllin, altriv, next, second, print_header, modchl
        LOGICAL :: iprcnd, munks , seprec, densep, calcdi, xactcp, reusec
        LOGICAL :: gmpspr, slvbqp, refact, fdgrad, centrl, dprcnd, strctr
        LOGICAL :: use_band, icfs, mortor, firsup, twonrm, direct, myprec
        LOGICAL :: prcond, firstc, nobnds, getders, save_c
        LOGICAL :: printt, printi, printm, printw, printd, printe, set_printe
        LOGICAL :: set_printt, set_printi, set_printm, set_printw, set_printd
        LOGICAL :: skipg, steering, new_major
        CHARACTER ( LEN = 6 ) :: cgend, lisend
        CHARACTER ( LEN = 1 ) :: cgend1, lisend1
        REAL ( KIND = KIND( 1.0E0 ) ) :: t, time, tmv, tca, tls, tup
        INTEGER ( KIND = ip_ ), DIMENSION( 5 ) :: ISYS
        CHARACTER ( LEN = 6 ), DIMENSION( 6 ) :: CGENDS
        CHARACTER ( LEN = 6 ), DIMENSION( 5 ) :: LSENDS
        CHARACTER ( LEN = 1 ), DIMENSION( 6 ) :: CGENDS1
        CHARACTER ( LEN = 1 ), DIMENSION( 5 ) :: LSENDS1
        TYPE( CAUCHY_save_type ) :: CAUCHY
        TYPE( CG_save_type ) :: CG
        TYPE( ASMBL_save_type ) :: ASMBL
        TYPE( PRECN_save_type ) :: PRECN
        TYPE( OTHERS_fdgrad_save_type ) :: OTHERS
        TYPE( EXTEND_save_type ) :: EXTEND
    END TYPE LANCELOT_save_type

    TYPE :: LANCELOT_data_type
        TYPE( LANCELOT_save_type ) :: S
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ITRANS
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ROW_start
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: POS_in_H
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: USED
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: FILLED
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: LINK_elem_uses_var
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: WTRANS
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISYMMD
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISWKSP
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISTAJC
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISTAGV
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISVGRP
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISLGRP
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IGCOLJ
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IVALJR
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IUSED
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ITYPER
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISSWTR
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISSITR
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISET
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: ISVSET
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: INVSET
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IFREE
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: INDEX
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IFREEC
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: INNONZ
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: LIST_elements
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : , : ) :: ISYMMH
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: FUVALS_temp
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: P
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: X0
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: XCP
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: GX0
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: RADII
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: DELTAX
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: QGRAD
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: GRJAC
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: CDASH
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: C2DASH
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: GV_old
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : , : ) :: BND
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : , : ) :: BND_radius
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IW_asmbl
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: NZ_comp_w
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: W_ws
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: W_el
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: W_in
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: H_el
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: H_in
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : , : ) :: IKEEP
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : , : ) :: IW1
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IW
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: IVUSE
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: H_col_ptr
        INTEGER ( KIND = ip_ ), ALLOCATABLE, DIMENSION( : ) :: L_col_ptr
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: W
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: RHS
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: RHS2
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: P2
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: G
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: DIAG
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: BREAKP
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : ) :: GRAD
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : , : ) :: W1
        REAL ( KIND = rp_ ), ALLOCATABLE, DIMENSION( : , : ) :: OFFDIA
        REAL ( KIND = rp_ ), POINTER, DIMENSION( : ) :: GROUP_SCALING => NULL( )
        LOGICAL, POINTER, DIMENSION( : ) :: GXEQX_AUG => NULL( )
        TYPE ( SCU_matrix_type ) :: SCU_matrix
        TYPE ( SCU_data_type ) :: SCU_data
        TYPE ( SMT_type ) :: matrix
        TYPE ( SILS_factors ) :: SILS_data
    END TYPE LANCELOT_data_type
end

subroutine test()
    use T
    interface
        SUBROUTINE LANCELOT_initialize( data, control )
            use T
            TYPE ( LANCELOT_data_type ), INTENT( INOUT ) :: data
            TYPE ( LANCELOT_control_type ), INTENT( OUT ) :: control
        end subroutine

        SUBROUTINE LANCELOT_solve( prob, RANGE , GVALS, FT, XT, FUVALS, lfuval,   &
                        ICALCF, ICALCG, IVAR, Q, DGRAD, control,       &
                        inform, data, ELFUN, GROUP, ELFUN_flexible,    &
                        ELDERS )
            use T

            TYPE ( LANCELOT_control_type ), INTENT( INOUT ) :: control
            TYPE ( LANCELOT_inform_type ), INTENT( INOUT ) :: inform
            TYPE ( LANCELOT_data_type ), INTENT( INOUT ) :: data
            TYPE ( LANCELOT_problem_type ), INTENT( INOUT ), TARGET :: prob
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: lfuval
            INTEGER ( KIND = ip_ ), INTENT( INOUT ), DIMENSION( prob%n  ) :: IVAR
            INTEGER ( KIND = ip_ ), INTENT( INOUT ), DIMENSION( prob%nel ) :: ICALCF
            INTEGER ( KIND = ip_ ), INTENT( INOUT ), DIMENSION( prob%ng ) :: ICALCG
            REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( prob%ng, 3 ) :: GVALS
            REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( prob%n ) :: Q, XT, DGRAD
            REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( prob%ng ) :: FT
            REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( lfuval ) :: FUVALS

            INTERFACE
                SUBROUTINE RANGE ( ielemn, transp, W1, W2, nelvar, ninvar, ieltyp,      &
                                    lw1, lw2 )
                    USE T

                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ielemn, nelvar, ninvar, ieltyp
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: lw1, lw2
                    LOGICAL, INTENT( IN ) :: transp
                    REAL ( KIND = rp_ ), INTENT( IN ), DIMENSION( lw1 ) :: W1
                    REAL ( KIND = rp_ ), DIMENSION( lw2 ) :: W2
                END SUBROUTINE RANGE

                SUBROUTINE ELFUN ( FUVALS, XVALUE, EPVALU, ncalcf, ITYPEE, ISTAEV,      &
                                    IELVAR, INTVAR, ISTADH, ISTEPA, ICALCF, ltypee,      &
                                    lstaev, lelvar, lntvar, lstadh, lstepa, lcalcf,      &
                                    lfuval, lxvalu, lepvlu, ifflag, ifstat )
                    USE T

                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ncalcf, ifflag, ltypee, lstaev
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: lelvar, lntvar, lstadh, lstepa
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: lcalcf, lfuval, lxvalu, lepvlu
                    INTEGER ( KIND = ip_ ), INTENT( OUT ) :: ifstat
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ITYPEE(ltypee), ISTAEV(lstaev)
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: IELVAR(lelvar), INTVAR(lntvar)
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ISTADH(lstadh), ISTEPA(lstepa)
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ICALCF(lcalcf)
                    REAL ( KIND = rp_ ), INTENT( IN ) :: XVALUE(lxvalu)
                    REAL ( KIND = rp_ ), INTENT( IN ) :: EPVALU(lepvlu)
                    REAL ( KIND = rp_ ), INTENT( INOUT ) :: FUVALS(lfuval)
                END SUBROUTINE ELFUN

                SUBROUTINE ELFUN_flexible(                                              &
                                    FUVALS, XVALUE, EPVALU, ncalcf, ITYPEE, ISTAEV,      &
                                    IELVAR, INTVAR, ISTADH, ISTEPA, ICALCF, ltypee,      &
                                    lstaev, lelvar, lntvar, lstadh, lstepa, lcalcf,      &
                                    lfuval, lxvalu, lepvlu, llders, ifflag, ELDERS,      &
                                    ifstat )
                    USE T

                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ncalcf, ifflag, ltypee, lstaev
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: lelvar, lntvar, lstadh, lstepa
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: lcalcf, lfuval, lxvalu, lepvlu
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: llders
                    INTEGER ( KIND = ip_ ), INTENT( OUT ) :: ifstat
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ITYPEE(ltypee), ISTAEV(lstaev)
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: IELVAR(lelvar), INTVAR(lntvar)
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ISTADH(lstadh), ISTEPA(lstepa)
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: ICALCF(lcalcf), ELDERS(2,llders)
                    REAL ( KIND = rp_ ), INTENT( IN ) :: XVALUE(lxvalu)
                    REAL ( KIND = rp_ ), INTENT( IN ) :: EPVALU(lepvlu)
                    REAL ( KIND = rp_ ), INTENT( INOUT ) :: FUVALS(lfuval)
                END SUBROUTINE ELFUN_flexible

                SUBROUTINE GROUP ( GVALUE, lgvalu, FVALUE, GPVALU, ncalcg,              &
                                    ITYPEG, ISTGPA, ICALCG, ltypeg, lstgpa,              &
                                    lcalcg, lfvalu, lgpvlu, derivs, igstat )
                    USE T

                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: lgvalu, ncalcg, ltypeg, lstgpa
                    INTEGER ( KIND = ip_ ), INTENT( IN ) :: lcalcg, lfvalu, lgpvlu
                    INTEGER ( KIND = ip_ ), INTENT( OUT ) :: igstat
                    LOGICAL, INTENT( IN ) :: derivs
                    INTEGER ( KIND = ip_ ), INTENT( IN ), DIMENSION( ltypeg ) :: ITYPEG
                    INTEGER ( KIND = ip_ ), INTENT( IN ), DIMENSION( lstgpa ) :: ISTGPA
                    INTEGER ( KIND = ip_ ), INTENT( IN ), DIMENSION( lcalcg ) :: ICALCG
                    REAL ( KIND = rp_ ), INTENT( IN ), DIMENSION( lfvalu ) :: FVALUE
                    REAL ( KIND = rp_ ), INTENT( IN ), DIMENSION( lgpvlu ) :: GPVALU
                    REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( lgvalu, 3 ) :: GVALUE
                END SUBROUTINE GROUP
            END INTERFACE

            INTEGER ( KIND = ip_ ), INTENT( INOUT ), OPTIONAL,                        &
                                                    DIMENSION( 2, prob%nel ) :: ELDERS
            OPTIONAL :: ELFUN, ELFUN_flexible, GROUP
        end subroutine
    end interface

    INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
    REAL ( KIND = wp ), PARAMETER :: infinity = 10.0_wp ** 20
    TYPE ( LANCELOT_control_type ) :: control
    TYPE ( LANCELOT_inform_type ) :: info
    TYPE ( LANCELOT_data_type ) :: data
    TYPE ( LANCELOT_problem_type ) :: prob
    INTEGER :: i, lfuval
    INTEGER, PARAMETER :: n = 3, ng = 6, nel = 3, nnza = 4, nvrels = 7
    INTEGER, PARAMETER :: ntotel = 3, ngpvlu = 0, nepvlu = 0
    INTEGER :: IVAR( n ), ICALCF( nel ), ICALCG( ng )
    REAL ( KIND = wp ) :: Q( n ), XT( n ), DGRAD( n ), FVALUE( ng )
    REAL ( KIND = wp ) :: GVALUE( ng, 3 )
    REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: FUVALS
    EXTERNAL RANGE, ELFUN, GROUP

    prob%n = n ; prob%ng = ng ; prob%nel = nel
    ALLOCATE( prob%ISTADG( prob%ng + 1 ), prob%ISTGPA( prob%ng + 1 ) )
    ALLOCATE( prob%ISTADA( prob%ng + 1 ), prob%ISTAEV( prob%nel + 1 ) )
    ALLOCATE( prob%ISTEPA( prob%nel + 1 ), prob%ITYPEG( prob%ng ) )
    ALLOCATE( prob%KNDOFG( prob%ng ), prob%ITYPEE( prob%nel ) )
    ALLOCATE( prob%INTVAR( prob%nel + 1 ) )
    ALLOCATE( prob%IELING( ntotel ), prob%IELVAR( nvrels ), prob%ICNA( nnza ) )
    ALLOCATE( prob%ISTADH( prob%nel + 1 ), prob%A( nnza ) )
    ALLOCATE( prob%B( prob%ng ), prob%BL( prob%n ), prob%BU( prob%n ) )
    ALLOCATE( prob%X( prob%n ), prob%Y( prob%ng ), prob%C( prob%ng ) )
    ALLOCATE( prob%GPVALU( ngpvlu ), prob%EPVALU( nepvlu ) )
    ALLOCATE( prob%ESCALE( ntotel ), prob%GSCALE( prob%ng ) )
    ALLOCATE( prob%VSCALE( prob%n ) )
    ALLOCATE( prob%INTREP( prob%nel ), prob%GXEQX( prob%ng ) )
    ALLOCATE( prob%GNAMES( prob%ng ), prob%VNAMES( prob%n ) )
    prob%ISTADG = (/ 1, 1, 2, 3, 3, 4, 4 /)
    prob%IELVAR = (/ 2, 1, 3, 2, 3, 1, 2 /)
    prob%ISTAEV = (/ 1, 4, 6, 8 /) ; prob%INTVAR( : nel ) = (/ 2, 2, 2 /)
    prob%IELING = (/ 1, 2, 3 /) ; prob%ICNA = (/ 1, 2, 1, 2 /)
    prob%ISTADA = (/ 1, 2, 2, 2, 3, 3, 5 /)
    prob%KNDOFG = (/ 1, 1, 1, 1, 1, 2 /)
    prob%ITYPEG = (/ 1, 0, 2, 0, 1, 3 /) ; prob%ITYPEE = (/ 1, 2, 2 /)
    prob%ISTGPA = (/ 0, 0, 0, 0, 0, 0 /) ; prob%ISTEPA = (/ 0, 0, 0 /)
    prob%A = (/ 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp /)
    prob%B = (/ 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp /)
    prob%BL = (/ - infinity, -1.0_wp, 1.0_wp /)
    prob%BU = (/ infinity, 1.0_wp, 2.0_wp /)
    prob%GSCALE = (/ 1.0_wp, 1.0_wp, 3.0_wp, 1.0_wp, 2.0_wp, 1.0_wp /)
    prob%ESCALE = (/ 1.0_wp, 1.0_wp, 1.0_wp /)
    prob%VSCALE = (/ 1.0_wp, 1.0_wp, 1.0_wp /)
    prob%X = (/ 0.0_wp, 0.0_wp, 1.5_wp /) ; prob%Y( 6 ) = 0.0_wp
    prob%INTREP = (/ .TRUE., .FALSE., .FALSE. /)
    prob%GXEQX = (/ .FALSE., .TRUE., .FALSE., .TRUE., .FALSE., .FALSE. /)
    lfuval = prob%nel + 2 * prob%n
    DO i = 1, prob%nel
        lfuval = lfuval + ( prob%INTVAR( i ) * ( prob%INTVAR( i ) + 3 ) ) / 2
    END DO
    ALLOCATE( FUVALS( lfuval ) )

    CALL LANCELOT_initialize( data, control )
    control%maxit = 100 ; control%out = 6
    control%stopg = 0.00001_wp ; control%stopc = 0.00001_wp
    control%linear_solver = 1
    control%exact_gcp = .FALSE.
    info%status = 0

    CALL LANCELOT_solve( &
    prob, RANGE, GVALUE, FVALUE, XT, FUVALS, lfuval, ICALCF, ICALCG, &
    IVAR, Q, DGRAD, control, info, data, ELFUN = ELFUN , GROUP = GROUP )
end subroutine