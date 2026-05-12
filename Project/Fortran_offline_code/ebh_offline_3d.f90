
!=============================================================================
! MODULE 1 – kinds
!=============================================================================
module ebh_kinds
  implicit none
  integer, parameter :: rk = selected_real_kind(15,300)
  integer, parameter :: NV = 6    ! Voigt size (strain/stress components)
  integer, parameter :: NG = 8    ! Gauss points per hex8 (2x2x2)
end module ebh_kinds

!=============================================================================
! MODULE 2 – Gmsh mesh reader  (format 2 ASCII)
!
! The reader fills:
!   nodes(1:n_nodes, 1:3)     – coordinates (x,y,z)
!   hex_conn(1:n_hex, 1:8)    – connectivity (1-based node indices)
!   hex_tag(1:n_hex)          – physical-group tag of each hex element
!   face_tag(1:n_face, 1:4)   – quad face connectivity (4 nodes)
!   face_phys(1:n_face)       – physical tag of each quad face
!
! Expected physical tags (user sets these in Gmsh):
!   Volumes : tag = phase number  (e.g. 1=matrix, 2=fiber)
!   Surfaces: +x=1, -x=2, +y=3, -y=4, +z=5, -z=6
!             (slaves = -x,-y,-z i.e. tags 2,4,6)
!=============================================================================
module ebh_mesh
  use ebh_kinds
  implicit none
  private
  public :: read_gmsh2, MeshData

  type MeshData
    integer                    :: n_nodes, n_hex, n_face
    real(rk), allocatable      :: coords(:,:)    ! (n_nodes, 3)
    integer,  allocatable      :: hex_conn(:,:)  ! (n_hex,   8)
    integer,  allocatable      :: hex_tag(:)     ! (n_hex)      phase tag
    integer,  allocatable      :: face_conn(:,:) ! (n_face,  4)
    integer,  allocatable      :: face_phys(:)   ! (n_face)     surface tag
  end type MeshData

contains

  subroutine read_gmsh2(filename, m)
    character(len=*), intent(in) :: filename
    type(MeshData),   intent(out) :: m

    integer  :: u, ios, i, j
    integer  :: n_nodes_file, n_elems_file
    integer  :: elem_num, elem_type, n_tags, phys_tag, geom_tag
    integer  :: n_hex_tmp, n_face_tmp
    character(len=256) :: line

    ! ---- temporary dynamic storage ----
    integer, parameter :: MAXEL = 500000
    integer, allocatable :: tmp_hex(:,:), tmp_htag(:)
    integer, allocatable :: tmp_face(:,:), tmp_ftag(:)
    integer :: conn8(8), conn4(4), extra_tag

    allocate(tmp_hex (MAXEL, 8), tmp_htag(MAXEL))
    allocate(tmp_face(MAXEL, 4), tmp_ftag(MAXEL))
    n_hex_tmp  = 0
    n_face_tmp = 0

    open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'read_gmsh2: cannot open file'

    ! ---- scan until $Nodes ----
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) stop 'read_gmsh2: $Nodes not found'
      if (trim(adjustl(line)) == '$Nodes') exit
    end do

    read(u,*) n_nodes_file
    m%n_nodes = n_nodes_file
    allocate(m%coords(n_nodes_file, 3))

    do i = 1, n_nodes_file
      read(u,*,iostat=ios) j, m%coords(i,1), m%coords(i,2), m%coords(i,3)
      if (ios /= 0) stop 'read_gmsh2: error reading nodes'
    end do

    ! skip $EndNodes, scan until $Elements
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) stop 'read_gmsh2: $Elements not found'
      if (trim(adjustl(line)) == '$Elements') exit
    end do

    read(u,*) n_elems_file

    do i = 1, n_elems_file
      read(u,*,iostat=ios) elem_num, elem_type, n_tags, phys_tag, geom_tag, &
                           (extra_tag, j=1, n_tags-2)
      ! rewind and re-read properly with connectivity
      backspace(u)

      if (elem_type == 5) then
        ! Hex8: type 5
        read(u,*) elem_num, elem_type, n_tags, phys_tag, &
                  (extra_tag, j=1, n_tags-1), conn8
        n_hex_tmp = n_hex_tmp + 1
        if (n_hex_tmp > MAXEL) stop 'read_gmsh2: too many hexahedra'
        tmp_hex(n_hex_tmp,:) = conn8
        tmp_htag(n_hex_tmp)  = phys_tag

      else if (elem_type == 3) then
        ! Quad4: type 3
        read(u,*) elem_num, elem_type, n_tags, phys_tag, &
                  (extra_tag, j=1, n_tags-1), conn4
        n_face_tmp = n_face_tmp + 1
        if (n_face_tmp > MAXEL) stop 'read_gmsh2: too many quads'
        tmp_face(n_face_tmp,:) = conn4
        tmp_ftag(n_face_tmp)   = phys_tag

      else
        ! skip other element types (lines, points, etc.)
        read(u,'(A)') line
      end if
    end do

    close(u)

    m%n_hex  = n_hex_tmp
    m%n_face = n_face_tmp

    allocate(m%hex_conn (n_hex_tmp,  8))
    allocate(m%hex_tag  (n_hex_tmp))
    allocate(m%face_conn(n_face_tmp, 4))
    allocate(m%face_phys(n_face_tmp))

    m%hex_conn  = tmp_hex (1:n_hex_tmp,  :)
    m%hex_tag   = tmp_htag(1:n_hex_tmp)
    m%face_conn = tmp_face(1:n_face_tmp, :)
    m%face_phys = tmp_ftag(1:n_face_tmp)

    deallocate(tmp_hex, tmp_htag, tmp_face, tmp_ftag)

    write(*,'(A,I7,A,I7,A,I7)') &
      '  Nodes: ', m%n_nodes, '   Hex8: ', m%n_hex, '   Quad4 faces: ', m%n_face

  end subroutine read_gmsh2

end module ebh_mesh

!=============================================================================
! MODULE 3 – Hex8 element routines
!
! Voigt order: [eps11, eps22, eps33, eps12, eps23, eps13]
!              (engineering shear in positions 4-6)
!=============================================================================
module ebh_hex8
  use ebh_kinds
  implicit none
  private
  public :: gauss_pts_hex8, shape_hex8, Bmat_hex8, Kelem_hex8, &
            Felem_eigenstrain, vol_hex8, isotropic_D


contains

  !---------------------------------------------------------------------------
  ! 2-point Gauss rule for hex8 (2x2x2 = 8 points)
  !---------------------------------------------------------------------------
  subroutine gauss_pts_hex8(xi_g, w_g)
    real(rk), intent(out) :: xi_g(NG,3), w_g(NG)
    real(rk), parameter   :: a = 1.0_rk/sqrt(3.0_rk), w1 = 1.0_rk
    integer :: i, j, k, n
    n = 0
    do k = -1, 1, 2
    do j = -1, 1, 2
    do i = -1, 1, 2
      n = n + 1
      xi_g(n,1) = real(i,rk)*a
      xi_g(n,2) = real(j,rk)*a
      xi_g(n,3) = real(k,rk)*a
      w_g(n)    = w1
    end do; end do; end do
  end subroutine gauss_pts_hex8

  !---------------------------------------------------------------------------
  ! Shape functions and their natural-coordinate derivatives for hex8
  !  N(1:8)          – shape function values at (xi,eta,zeta)
  !  dN(1:8, 1:3)    – dN/d(xi,eta,zeta)
  !  Node ordering: Gmsh hex8 numbering (right-hand rule)
  !    1=(-1,-1,-1), 2=(+1,-1,-1), 3=(+1,+1,-1), 4=(-1,+1,-1)
  !    5=(-1,-1,+1), 6=(+1,-1,+1), 7=(+1,+1,+1), 8=(-1,+1,+1)
  !---------------------------------------------------------------------------
  subroutine shape_hex8(xi, eta, zeta, N, dN)
    real(rk), intent(in)  :: xi, eta, zeta
    real(rk), intent(out) :: N(8), dN(8,3)
    real(rk) :: xm,xp, ym,yp, zm,zp
    xm=1.0_rk-xi; xp=1.0_rk+xi
    ym=1.0_rk-eta;yp=1.0_rk+eta
    zm=1.0_rk-zeta;zp=1.0_rk+zeta

    N(1)=0.125_rk*xm*ym*zm; N(2)=0.125_rk*xp*ym*zm
    N(3)=0.125_rk*xp*yp*zm; N(4)=0.125_rk*xm*yp*zm
    N(5)=0.125_rk*xm*ym*zp; N(6)=0.125_rk*xp*ym*zp
    N(7)=0.125_rk*xp*yp*zp; N(8)=0.125_rk*xm*yp*zp

    ! dN/dxi
    dN(1,1)=-0.125_rk*ym*zm; dN(2,1)= 0.125_rk*ym*zm
    dN(3,1)= 0.125_rk*yp*zm; dN(4,1)=-0.125_rk*yp*zm
    dN(5,1)=-0.125_rk*ym*zp; dN(6,1)= 0.125_rk*ym*zp
    dN(7,1)= 0.125_rk*yp*zp; dN(8,1)=-0.125_rk*yp*zp
    ! dN/deta
    dN(1,2)=-0.125_rk*xm*zm; dN(2,2)=-0.125_rk*xp*zm
    dN(3,2)= 0.125_rk*xp*zm; dN(4,2)= 0.125_rk*xm*zm
    dN(5,2)=-0.125_rk*xm*zp; dN(6,2)=-0.125_rk*xp*zp
    dN(7,2)= 0.125_rk*xp*zp; dN(8,2)= 0.125_rk*xm*zp
    ! dN/dzeta
    dN(1,3)=-0.125_rk*xm*ym; dN(2,3)=-0.125_rk*xp*ym
    dN(3,3)=-0.125_rk*xp*yp; dN(4,3)=-0.125_rk*xm*yp
    dN(5,3)= 0.125_rk*xm*ym; dN(6,3)= 0.125_rk*xp*ym
    dN(7,3)= 0.125_rk*xp*yp; dN(8,3)= 0.125_rk*xm*yp
  end subroutine shape_hex8

  !---------------------------------------------------------------------------
  ! B matrix (6x24) and Jacobian determinant at a single Gauss point
  !  coords(8,3)  – nodal coordinates of element
  !---------------------------------------------------------------------------
  subroutine Bmat_hex8(xi, eta, zeta, coords, B, detJ)
    real(rk), intent(in)  :: xi, eta, zeta, coords(8,3)
    real(rk), intent(out) :: B(NV,24), detJ

    real(rk) :: N(8), dN(8,3), J(3,3), Jinv(3,3), dNdx(8,3)
    integer  :: a, d

    call shape_hex8(xi, eta, zeta, N, dN)

    ! Jacobian  J_{ij} = sum_a  dN_a/dxi_j * x_{a,i}
    J = 0.0_rk
    do a = 1, 8
      do d = 1, 3
        J(1,d) = J(1,d) + dN(a,d)*coords(a,1)
        J(2,d) = J(2,d) + dN(a,d)*coords(a,2)
        J(3,d) = J(3,d) + dN(a,d)*coords(a,3)
      end do
    end do

    call inv3(J, Jinv, detJ)
    if (detJ <= 0.0_rk) stop 'Bmat_hex8: negative Jacobian – check element orientation'

    ! dN/dx_phys = dN/dxi * J^{-1}
  dNdx = matmul(dN, Jinv)

    B = 0.0_rk
    do a = 1, 8
      ! rows: [eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps13]  (engineering shear)
      B(1, 3*(a-1)+1) = dNdx(a,1)
      B(2, 3*(a-1)+2) = dNdx(a,2)
      B(3, 3*(a-1)+3) = dNdx(a,3)
      B(4, 3*(a-1)+1) = dNdx(a,2);  B(4, 3*(a-1)+2) = dNdx(a,1)
      B(5, 3*(a-1)+2) = dNdx(a,3);  B(5, 3*(a-1)+3) = dNdx(a,2)
      B(6, 3*(a-1)+1) = dNdx(a,3);  B(6, 3*(a-1)+3) = dNdx(a,1)
    end do
  end subroutine Bmat_hex8

  !---------------------------------------------------------------------------
  ! Element stiffness matrix  Ke(24x24)  for isotropic material
  !  E, nu – Young's modulus and Poisson ratio
  !---------------------------------------------------------------------------
  subroutine Kelem_hex8(coords, E_mod, nu, Ke)
    real(rk), intent(in)  :: coords(8,3), E_mod, nu
    real(rk), intent(out) :: Ke(24,24)

    real(rk) :: xi_g(NG,3), w_g(NG), B(NV,24), detJ, D(NV,NV)
    integer  :: ig

    call gauss_pts_hex8(xi_g, w_g)
    call isotropic_D(E_mod, nu, D)
    Ke = 0.0_rk
    do ig = 1, NG
      call Bmat_hex8(xi_g(ig,1),xi_g(ig,2),xi_g(ig,3), coords, B, detJ)
      Ke = Ke + w_g(ig)*detJ * matmul(transpose(B), matmul(D,B))
    end do
  end subroutine Kelem_hex8

  !---------------------------------------------------------------------------
  ! Element force vector due to eigenstrain mu(6) on element
  !  F_e = integral B^T D mu dV     (24-vector)
  !---------------------------------------------------------------------------
  subroutine Felem_eigenstrain(coords, E_mod, nu, mu, Fe)
    real(rk), intent(in)  :: coords(8,3), E_mod, nu, mu(NV)
    real(rk), intent(out) :: Fe(24)

    real(rk) :: xi_g(NG,3), w_g(NG), B(NV,24), detJ, D(NV,NV)
    integer  :: ig

    call gauss_pts_hex8(xi_g, w_g)
    call isotropic_D(E_mod, nu, D)
    Fe = 0.0_rk
    do ig = 1, NG
      call Bmat_hex8(xi_g(ig,1),xi_g(ig,2),xi_g(ig,3), coords, B, detJ)
      Fe = Fe + w_g(ig)*detJ * matmul(transpose(B), matmul(D, mu))
    end do
  end subroutine Felem_eigenstrain

  !---------------------------------------------------------------------------
  ! Volume of a hex8 element  (sum of Gauss-point contributions)
  !---------------------------------------------------------------------------
  real(rk) function vol_hex8(coords)
    real(rk), intent(in) :: coords(8,3)
    real(rk) :: xi_g(NG,3), w_g(NG), B(NV,24), detJ
    integer  :: ig
    call gauss_pts_hex8(xi_g, w_g)
    vol_hex8 = 0.0_rk
    do ig = 1, NG
      call Bmat_hex8(xi_g(ig,1),xi_g(ig,2),xi_g(ig,3), coords, B, detJ)
      vol_hex8 = vol_hex8 + w_g(ig)*detJ
    end do
  end function vol_hex8

  !---------------------------------------------------------------------------
  ! Isotropic constitutive matrix D (6x6) in Voigt notation
  !---------------------------------------------------------------------------
  subroutine isotropic_D(E_mod, nu, D)
    real(rk), intent(in)  :: E_mod, nu
    real(rk), intent(out) :: D(NV,NV)
    real(rk) :: c, s, lam, mu2
    c   = E_mod / ((1.0_rk+nu)*(1.0_rk-2.0_rk*nu))
    lam = c * nu
    mu2 = c * (1.0_rk-2.0_rk*nu) * 0.5_rk   ! = G
    D = 0.0_rk
    D(1,1)=c*(1.0_rk-nu); D(1,2)=lam;          D(1,3)=lam
    D(2,1)=lam;           D(2,2)=c*(1.0_rk-nu); D(2,3)=lam
    D(3,1)=lam;           D(3,2)=lam;           D(3,3)=c*(1.0_rk-nu)
    D(4,4)=mu2
    D(5,5)=mu2
    D(6,6)=mu2
  end subroutine isotropic_D

  !---------------------------------------------------------------------------
  ! 3x3 matrix inverse (Cramer's rule) and determinant
  !---------------------------------------------------------------------------
  subroutine inv3(A, Ainv, detA)
    real(rk), intent(in)  :: A(3,3)
    real(rk), intent(out) :: Ainv(3,3), detA
    detA = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
         - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
         + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    if (abs(detA) < 1.0e-30_rk) stop 'inv3: singular matrix'
    Ainv(1,1)= (A(2,2)*A(3,3)-A(2,3)*A(3,2))/detA
    Ainv(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))/detA
    Ainv(1,3)= (A(1,2)*A(2,3)-A(1,3)*A(2,2))/detA
    Ainv(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))/detA
    Ainv(2,2)= (A(1,1)*A(3,3)-A(1,3)*A(3,1))/detA
    Ainv(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))/detA
    Ainv(3,1)= (A(2,1)*A(3,2)-A(2,2)*A(3,1))/detA
    Ainv(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))/detA
    Ainv(3,3)= (A(1,1)*A(2,2)-A(1,2)*A(2,1))/detA
  end subroutine inv3

end module ebh_hex8

!=============================================================================
! MODULE 4 – Periodic boundary conditions (3D, Lagrange multiplier method)
!
! HIERARCHICAL VERSION:
!   1) 8 corners fixed to zero
!   2) strict edge nodes paired by edge families
!   3) face-interior nodes paired across opposite faces
!
! This avoids mixing categories while enforcing full 3D periodicity.
!
! Face physical tags expected in the mesh:
!   MASTER_X = 1, SLAVE_X = 2
!   MASTER_Y = 3, SLAVE_Y = 4
!   MASTER_Z = 5, SLAVE_Z = 6
!=============================================================================
module ebh_pbc
  use ebh_kinds
  use ebh_mesh, only: MeshData
  implicit none
  private
  public :: build_pbc_constraints, MASTER_X, SLAVE_X, &
            MASTER_Y, SLAVE_Y, MASTER_Z, SLAVE_Z

  integer, parameter :: MASTER_X = 1, SLAVE_X = 2
  integer, parameter :: MASTER_Y = 3, SLAVE_Y = 4
  integer, parameter :: MASTER_Z = 5, SLAVE_Z = 6

contains

  subroutine build_pbc_constraints(m, n_dof_tot, slave_dof, master_dof, &
                                  corner_dof, C_mat, RHS)
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: n_dof_tot
    integer, allocatable, intent(out) :: slave_dof(:), master_dof(:)
    integer, allocatable, intent(out) :: corner_dof(:)
    real(rk), allocatable, intent(out) :: C_mat(:,:), RHS(:)

    integer, allocatable :: corner_nodes(:)
    integer, allocatable :: s_tmp(:), m_tmp(:)
    integer :: n_corner, n_cdofs
    integer :: n_pairs, maxpairs, n_constr
    integer :: i

    maxpairs = 3 * m%n_nodes * 6
    allocate(s_tmp(maxpairs), m_tmp(maxpairs))
    n_pairs = 0

    !--------------------------------------------------
    ! 1) corners
    !--------------------------------------------------
    call collect_corner_nodes_by_bbox(m%coords, m%n_nodes, corner_nodes, n_corner)
    if (n_corner /= 8) then
      write(*,'(A,I0,A)') 'WARNING: expected 8 corners, detected = ', n_corner, ''
      stop 'build_pbc_constraints: corner detection failed'
    end if

    n_cdofs = 3*n_corner
    allocate(corner_dof(n_cdofs))
    call build_corner_dofs(corner_nodes, n_corner, corner_dof)

    !--------------------------------------------------
    ! 2) strict edge nodes
    !--------------------------------------------------
    call add_edge_family_constraints(m%coords, m%n_nodes, 1, s_tmp, m_tmp, n_pairs) ! x-edges
    call add_edge_family_constraints(m%coords, m%n_nodes, 2, s_tmp, m_tmp, n_pairs) ! y-edges
    call add_edge_family_constraints(m%coords, m%n_nodes, 3, s_tmp, m_tmp, n_pairs) ! z-edges

    !--------------------------------------------------
    ! 3) face-interior nodes
    !--------------------------------------------------
    call add_face_constraints(m, MASTER_X, SLAVE_X, 1, s_tmp, m_tmp, n_pairs)
    call add_face_constraints(m, MASTER_Y, SLAVE_Y, 2, s_tmp, m_tmp, n_pairs)
    call add_face_constraints(m, MASTER_Z, SLAVE_Z, 3, s_tmp, m_tmp, n_pairs)
  

    allocate(slave_dof(n_pairs), master_dof(n_pairs))
    slave_dof  = s_tmp(1:n_pairs)
    master_dof = m_tmp(1:n_pairs)

    !--------------------------------------------------
    ! 4) assemble C
    !--------------------------------------------------
    n_constr = n_pairs + n_cdofs
    allocate(C_mat(n_constr, n_dof_tot), RHS(n_constr))
    C_mat = 0.0_rk
    RHS   = 0.0_rk

    do i = 1, n_pairs
      C_mat(i, slave_dof(i))  =  1.0_rk
      C_mat(i, master_dof(i)) = -1.0_rk
    end do

    do i = 1, n_cdofs
      C_mat(n_pairs+i, corner_dof(i)) = 1.0_rk
    end do

    write(*,'(A,I0)') '  Corner nodes detected   = ', n_corner
    write(*,'(A,I0)') '  PBC pairs (DOF-level)   = ', n_pairs
    write(*,'(A,I0)') '  Corner DOFs fixed       = ', n_cdofs
    write(*,'(A,I0)') '  Total constraints       = ', n_constr

    deallocate(corner_nodes, s_tmp, m_tmp)
  end subroutine build_pbc_constraints

  !===========================================================================
  ! Geometry helpers
  !===========================================================================

  subroutine get_bbox(coords, xmin,xmax,ymin,ymax,zmin,zmax,tol)
    real(rk), intent(in)  :: coords(:,:)
    real(rk), intent(out) :: xmin,xmax,ymin,ymax,zmin,zmax,tol
    real(rk) :: Lref

    xmin = minval(coords(:,1)); xmax = maxval(coords(:,1))
    ymin = minval(coords(:,2)); ymax = maxval(coords(:,2))
    zmin = minval(coords(:,3)); zmax = maxval(coords(:,3))

    Lref = max( max(xmax-xmin, ymax-ymin), zmax-zmin )
    tol  = 1.0e-8_rk * max(Lref, 1.0_rk)
    if (tol == 0.0_rk) tol = 1.0e-12_rk
  end subroutine get_bbox

  integer function n_boundary_planes_of_node(x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,tol)
    real(rk), intent(in) :: x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,tol
    n_boundary_planes_of_node = 0
    if (abs(x-xmin) < tol .or. abs(x-xmax) < tol) n_boundary_planes_of_node = n_boundary_planes_of_node + 1
    if (abs(y-ymin) < tol .or. abs(y-ymax) < tol) n_boundary_planes_of_node = n_boundary_planes_of_node + 1
    if (abs(z-zmin) < tol .or. abs(z-zmax) < tol) n_boundary_planes_of_node = n_boundary_planes_of_node + 1
  end function n_boundary_planes_of_node

  subroutine collect_corner_nodes_by_bbox(coords, n_nodes, nlist_corner, n_corner)
    real(rk), intent(in) :: coords(n_nodes,3)
    integer,  intent(in) :: n_nodes
    integer, allocatable, intent(out) :: nlist_corner(:)
    integer, intent(out) :: n_corner

    real(rk) :: xmin,xmax,ymin,ymax,zmin,zmax,tol
    integer  :: i, k, tmp(8)
    logical  :: on_x, on_y, on_z

    call get_bbox(coords, xmin,xmax,ymin,ymax,zmin,zmax,tol)

    k = 0
    tmp = 0
    do i = 1, n_nodes
      on_x = (abs(coords(i,1)-xmin) < tol) .or. (abs(coords(i,1)-xmax) < tol)
      on_y = (abs(coords(i,2)-ymin) < tol) .or. (abs(coords(i,2)-ymax) < tol)
      on_z = (abs(coords(i,3)-zmin) < tol) .or. (abs(coords(i,3)-zmax) < tol)
      if (on_x .and. on_y .and. on_z) then
        k = k + 1
        if (k > 8) stop 'collect_corner_nodes_by_bbox: more than 8 corners detected'
        tmp(k) = i
      end if
    end do

    n_corner = k
    allocate(nlist_corner(n_corner))
    if (n_corner > 0) nlist_corner = tmp(1:n_corner)
  end subroutine collect_corner_nodes_by_bbox

  subroutine build_corner_dofs(clist, nc, cdof)
    integer, intent(in)  :: clist(:), nc
    integer, intent(out) :: cdof(nc*3)
    integer :: i, d
    do i = 1, nc
      do d = 1, 3
        cdof(3*(i-1)+d) = 3*(clist(i)-1) + d
      end do
    end do
  end subroutine build_corner_dofs

  !===========================================================================
  ! Edge constraints
  !===========================================================================

  subroutine add_edge_family_constraints(coords, n_nodes, axis_dir, s_arr, m_arr, p)
    real(rk), intent(in)    :: coords(n_nodes,3)
    integer,  intent(in)    :: n_nodes, axis_dir
    integer,  intent(inout) :: s_arr(:), m_arr(:), p

    real(rk) :: xmin,xmax,ymin,ymax,zmin,zmax,tol
    integer, allocatable :: e1(:), e2(:), e3(:), e4(:)
    integer :: n1,n2,n3,n4

    call get_bbox(coords, xmin,xmax,ymin,ymax,zmin,zmax,tol)

    select case(axis_dir)
    case(1) ! edges parallel to x: fixed (y,z)
      call collect_strict_edge_nodes(coords, n_nodes, 1, ymin, zmin, e1, n1)
      call collect_strict_edge_nodes(coords, n_nodes, 1, ymax, zmin, e2, n2)
      call collect_strict_edge_nodes(coords, n_nodes, 1, ymin, zmax, e3, n3)
      call collect_strict_edge_nodes(coords, n_nodes, 1, ymax, zmax, e4, n4)

    case(2) ! edges parallel to y: fixed (x,z)
      call collect_strict_edge_nodes(coords, n_nodes, 2, xmin, zmin, e1, n1)
      call collect_strict_edge_nodes(coords, n_nodes, 2, xmax, zmin, e2, n2)
      call collect_strict_edge_nodes(coords, n_nodes, 2, xmin, zmax, e3, n3)
      call collect_strict_edge_nodes(coords, n_nodes, 2, xmax, zmax, e4, n4)

    case(3) ! edges parallel to z: fixed (x,y)
      call collect_strict_edge_nodes(coords, n_nodes, 3, xmin, ymin, e1, n1)
      call collect_strict_edge_nodes(coords, n_nodes, 3, xmax, ymin, e2, n2)
      call collect_strict_edge_nodes(coords, n_nodes, 3, xmin, ymax, e3, n3)
      call collect_strict_edge_nodes(coords, n_nodes, 3, xmax, ymax, e4, n4)

    case default
      stop 'add_edge_family_constraints: invalid axis_dir'
    end select

    if (n1 /= n2 .or. n1 /= n3 .or. n1 /= n4) stop 'add_edge_family_constraints: inconsistent number of edge nodes'

    ! Use edge 1 as master reference; slaves = edges 2, 3, 4
    call pair_node_lists_by_axis(coords, e1, n1, e2, n2, axis_dir, s_arr, m_arr, p)
    call pair_node_lists_by_axis(coords, e1, n1, e3, n3, axis_dir, s_arr, m_arr, p)
    call pair_node_lists_by_axis(coords, e1, n1, e4, n4, axis_dir, s_arr, m_arr, p)

    deallocate(e1,e2,e3,e4)
  end subroutine add_edge_family_constraints

  subroutine collect_strict_edge_nodes(coords, n_nodes, axis_dir, v1, v2, nlist, ncount)
    real(rk), intent(in) :: coords(n_nodes,3)
    integer,  intent(in) :: n_nodes, axis_dir
    real(rk), intent(in) :: v1, v2
    integer, allocatable, intent(out) :: nlist(:)
    integer, intent(out) :: ncount

    real(rk) :: xmin,xmax,ymin,ymax,zmin,zmax,tol
    integer :: i, nb, tmp(n_nodes)
    logical :: ok

    call get_bbox(coords, xmin,xmax,ymin,ymax,zmin,zmax,tol)

    ncount = 0
    tmp = 0

    do i = 1, n_nodes
      nb = n_boundary_planes_of_node(coords(i,1),coords(i,2),coords(i,3), xmin,xmax,ymin,ymax,zmin,zmax,tol)
      if (nb /= 2) cycle

      select case(axis_dir)
      case(1)
        ok = (abs(coords(i,2)-v1) < tol) .and. (abs(coords(i,3)-v2) < tol)
      case(2)
        ok = (abs(coords(i,1)-v1) < tol) .and. (abs(coords(i,3)-v2) < tol)
      case(3)
        ok = (abs(coords(i,1)-v1) < tol) .and. (abs(coords(i,2)-v2) < tol)
      case default
        stop 'collect_strict_edge_nodes: invalid axis_dir'
      end select

      if (ok) then
        ncount = ncount + 1
        tmp(ncount) = i
      end if
    end do

    allocate(nlist(ncount))
    if (ncount > 0) nlist = tmp(1:ncount)
  end subroutine collect_strict_edge_nodes

  subroutine pair_node_lists_by_axis(coords, master_nodes, nm, slave_nodes, ns, axis_dir, s_arr, m_arr, p)
    real(rk), intent(in)    :: coords(:,:)
    integer,  intent(in)    :: master_nodes(:), nm, slave_nodes(:), ns, axis_dir
    integer,  intent(inout) :: s_arr(:), m_arr(:), p

    real(rk) :: tol, xmin,xmax,ymin,ymax,zmin,zmax
    integer :: is, im, sl, ms, d
    logical :: matched

    if (nm /= ns) stop 'pair_node_lists_by_axis: nm /= ns'

    call get_bbox(coords, xmin,xmax,ymin,ymax,zmin,zmax,tol)

    do is = 1, ns
      sl = slave_nodes(is)
      matched = .false.
      do im = 1, nm
        ms = master_nodes(im)

        if (abs(coords(sl,axis_dir)-coords(ms,axis_dir)) < tol) then
          do d = 1, 3
            p = p + 1
            s_arr(p) = 3*(sl-1)+d
            m_arr(p) = 3*(ms-1)+d
          end do
          matched = .true.
          exit
        end if
      end do
      if (.not. matched) stop 'pair_node_lists_by_axis: could not match edge node'
    end do
  end subroutine pair_node_lists_by_axis

  !===========================================================================
  ! Face-interior constraints
  !===========================================================================

  subroutine add_face_constraints(m, master_tag, slave_tag, periodic_axis, s_arr, m_arr, p)
    type(MeshData), intent(in)    :: m
    integer,        intent(in)    :: master_tag, slave_tag, periodic_axis
    integer,        intent(inout) :: s_arr(:), m_arr(:), p

    integer, allocatable :: mlist(:), slist(:)
    integer :: nm, ns

    call collect_face_interior_nodes(m, master_tag, mlist, nm)
    call collect_face_interior_nodes(m, slave_tag,  slist, ns)

    if (nm /= ns) stop 'add_face_constraints: inconsistent number of face-interior nodes'

    call pair_face_lists(coords=m%coords, master_nodes=mlist, nm=nm, slave_nodes=slist, ns=ns, &
                         periodic_axis=periodic_axis, s_arr=s_arr, m_arr=m_arr, p=p)

    deallocate(mlist, slist)
  end subroutine add_face_constraints

  subroutine collect_face_interior_nodes(m, tag, nlist, ncount)
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: tag
    integer, allocatable, intent(out) :: nlist(:)
    integer, intent(out) :: ncount

    real(rk) :: xmin,xmax,ymin,ymax,zmin,zmax,tol
    integer  :: i, j, n, nb, tmp(m%n_face*4)
    logical  :: found

    call get_bbox(m%coords, xmin,xmax,ymin,ymax,zmin,zmax,tol)

    ncount = 0
    tmp = 0

    do i = 1, m%n_face
      if (m%face_phys(i) /= tag) cycle
      do j = 1, 4
        n = m%face_conn(i,j)

        nb = n_boundary_planes_of_node(m%coords(n,1),m%coords(n,2),m%coords(n,3), &
                                       xmin,xmax,ymin,ymax,zmin,zmax,tol)

        if (nb /= 1) cycle

        found = .false.
        block
          integer :: k
          do k = 1, ncount
            if (tmp(k) == n) then
              found = .true.
              exit
            end if
          end do
        end block

        if (.not. found) then
          ncount = ncount + 1
          tmp(ncount) = n
        end if
      end do
    end do

    allocate(nlist(ncount))
    if (ncount > 0) nlist = tmp(1:ncount)
  end subroutine collect_face_interior_nodes

  subroutine pair_face_lists(coords, master_nodes, nm, slave_nodes, ns, periodic_axis, s_arr, m_arr, p)
    real(rk), intent(in)    :: coords(:,:)
    integer,  intent(in)    :: master_nodes(:), nm, slave_nodes(:), ns, periodic_axis
    integer,  intent(inout) :: s_arr(:), m_arr(:), p

    real(rk) :: tol, xmin,xmax,ymin,ymax,zmin,zmax
    integer  :: is, im, sl, ms, d, a1, a2
    logical  :: matched

    if (nm /= ns) stop 'pair_face_lists: nm /= ns'

    call get_bbox(coords, xmin,xmax,ymin,ymax,zmin,zmax,tol)

    select case(periodic_axis)
    case(1)
      a1 = 2; a2 = 3
    case(2)
      a1 = 1; a2 = 3
    case(3)
      a1 = 1; a2 = 2
    case default
      stop 'pair_face_lists: invalid periodic_axis'
    end select

    do is = 1, ns
      sl = slave_nodes(is)
      matched = .false.
      do im = 1, nm
        ms = master_nodes(im)

        if (abs(coords(sl,a1)-coords(ms,a1)) < tol .and. &
            abs(coords(sl,a2)-coords(ms,a2)) < tol) then
          do d = 1, 3
            p = p + 1
            s_arr(p) = 3*(sl-1)+d
            m_arr(p) = 3*(ms-1)+d
          end do
          matched = .true.
          exit
        end if
      end do
      if (.not. matched) stop 'pair_face_lists: could not match face-interior node'
    end do
  end subroutine pair_face_lists

end module ebh_pbc

!=============================================================================
! MODULE 5 – Global assembly and solver
!=============================================================================
module ebh_assembly
  use ebh_kinds
  use ebh_mesh,  only: MeshData
  use ebh_hex8,  only: Kelem_hex8, Felem_eigenstrain, Bmat_hex8, &
                        gauss_pts_hex8, vol_hex8
  implicit none
  private
  public :: assemble_K, assemble_F_macro, assemble_F_eigenstrain, &
            solve_augmented, average_strain_partition, &
            get_elem_coords, get_elem_dofs

  integer, parameter :: NDOF_E = 24   ! DOFs per hex8 element

contains

  !---------------------------------------------------------------------------
  ! get_elem_coords: extract 8x3 nodal coordinates for element ie
  !---------------------------------------------------------------------------
  subroutine get_elem_coords(m, ie, coords_e)
    use ebh_mesh, only: MeshData
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: ie
    real(rk),       intent(out) :: coords_e(8,3)
    integer :: a, nd
    do a = 1, 8
      nd = m%hex_conn(ie,a)
      coords_e(a,:) = m%coords(nd,:)
    end do
  end subroutine get_elem_coords

  !---------------------------------------------------------------------------
  ! get_elem_dofs: map local DOF index (1..24) to global DOF (1..3*n_nodes)
  !---------------------------------------------------------------------------
  subroutine get_elem_dofs(m, ie, gdofs)
    use ebh_mesh, only: MeshData
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: ie
    integer,        intent(out) :: gdofs(NDOF_E)
    integer :: a, nd
    do a = 1, 8
      nd = m%hex_conn(ie,a)
      gdofs(3*(a-1)+1) = 3*(nd-1)+1
      gdofs(3*(a-1)+2) = 3*(nd-1)+2
      gdofs(3*(a-1)+3) = 3*(nd-1)+3
    end do
  end subroutine get_elem_dofs

  !---------------------------------------------------------------------------
  ! assemble_K: global stiffness matrix K(n_dof,n_dof)
  !   phase_E(tag), phase_nu(tag) – material properties indexed by physical tag
  !---------------------------------------------------------------------------
  subroutine assemble_K(m, n_phases, phase_E, phase_nu, K)
    use ebh_mesh, only: MeshData
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: n_phases
    real(rk),       intent(in)  :: phase_E(n_phases), phase_nu(n_phases)
    real(rk),       intent(out) :: K(:,:)

    real(rk) :: coords_e(8,3), Ke(NDOF_E,NDOF_E)
    integer  :: gdofs(NDOF_E)
    integer  :: ie, j, kc, tag

    K = 0.0_rk
    do ie = 1, m%n_hex
      tag = m%hex_tag(ie)
      if (tag < 1 .or. tag > n_phases) stop 'assemble_K: unknown phase tag'
      call get_elem_coords(m, ie, coords_e)
      call Kelem_hex8(coords_e, phase_E(tag), phase_nu(tag), Ke)
      call get_elem_dofs(m, ie, gdofs)
      do j = 1, NDOF_E
        do kc = 1, NDOF_E
          K(gdofs(j), gdofs(kc)) = K(gdofs(j), gdofs(kc)) + Ke(j,kc)
        end do
      end do
    end do
  end subroutine assemble_K

  !---------------------------------------------------------------------------
  ! assemble_F_macro:
  !   For each of the 6 unit macro-strain components eps_bar_kl,
  !   assemble the corresponding load vector F_kl(n_dof).
  !
  !   From TFA: the macro-strain enters as a body force equivalent:
  !     F_macro = -K_0 * H_macro  (see Fish 2013, eq. used to solve for H^kl)
  !   but more directly, the RVE problem with eps_bar prescribed gives:
  !     f^{kl}_i = -integral_Theta  B^T D eps_bar^{kl} dV   (per element)
  !   where eps_bar^{kl} has only component kl = 1 (Voigt unit vector).
  !
  !   Returns F_macro(n_dof, 6) – one column per macro-strain component.
  !---------------------------------------------------------------------------
  subroutine assemble_F_macro(m, n_phases, phase_E, phase_nu, F_macro)
    use ebh_mesh, only: MeshData
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: n_phases
    real(rk),       intent(in)  :: phase_E(n_phases), phase_nu(n_phases)
    real(rk),       intent(out) :: F_macro(:,:)   ! (n_dof, 6)

    real(rk) :: coords_e(8,3), Fe(NDOF_E), mu_unit(NV)
    integer  :: gdofs(NDOF_E), ie, kl, tag, j, i

    F_macro = 0.0_rk
    do kl = 1, NV
      mu_unit    = 0.0_rk
      mu_unit(kl) = 1.0_rk
      do ie = 1, m%n_hex
        tag = m%hex_tag(ie)
        call get_elem_coords(m, ie, coords_e)
        call Felem_eigenstrain(coords_e, phase_E(tag), phase_nu(tag), mu_unit, Fe)
        call get_elem_dofs(m, ie, gdofs)
        do j = 1, NDOF_E
          F_macro(gdofs(j), kl) = F_macro(gdofs(j), kl) + Fe(j)
        end do
      end do
    end do


    ! Note: the RVE problem is K u = -F_macro  (the minus comes from moving
    ! the initial-strain term to the RHS).  We store F_macro here; the sign
    ! is applied in solve_augmented.
  end subroutine assemble_F_macro

  !---------------------------------------------------------------------------
  ! assemble_F_eigenstrain:
  !   Assemble load vector for a UNIT eigenstrain on partition A, component kl.
  !   Returns F_eig(n_dof).
  !---------------------------------------------------------------------------
  subroutine assemble_F_eigenstrain(m, n_phases, phase_E, phase_nu, &
                                    part_id, part_A, kl, F_eig)
    use ebh_mesh, only: MeshData
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: n_phases
    real(rk),       intent(in)  :: phase_E(n_phases), phase_nu(n_phases)
    integer,        intent(in)  :: part_id(:)   ! partition ID per element
    integer,        intent(in)  :: part_A       ! target partition
    integer,        intent(in)  :: kl           ! Voigt component (1..6)
    real(rk),       intent(out) :: F_eig(:)     ! (n_dof)

    real(rk) :: coords_e(8,3), Fe(NDOF_E), mu_unit(NV)
    integer  :: gdofs(NDOF_E), ie, tag, j

    F_eig   = 0.0_rk
    mu_unit = 0.0_rk
    mu_unit(kl) = 1.0_rk

    do ie = 1, m%n_hex
      if (part_id(ie) /= part_A) cycle
      tag = m%hex_tag(ie)
      call get_elem_coords(m, ie, coords_e)
      call Felem_eigenstrain(coords_e, phase_E(tag), phase_nu(tag), mu_unit, Fe)
      call get_elem_dofs(m, ie, gdofs)
      do j = 1, NDOF_E
        F_eig(gdofs(j)) = F_eig(gdofs(j)) + Fe(j)
      end do
    end do
  end subroutine assemble_F_eigenstrain

  !---------------------------------------------------------------------------
  ! solve_augmented:
  !   Solves  [K   C^T] [u     ]   [f ]
  !           [C    0 ] [lambda] = [0 ]
  !   using LAPACK DGESV.
  !   Returns u(n_dof).
  !---------------------------------------------------------------------------
  subroutine solve_augmented(K, C_mat, f_rhs, n_dof, n_constr, u_out)
    real(rk), intent(in)  :: K(:,:), C_mat(:,:), f_rhs(:)
    integer,  intent(in)  :: n_dof, n_constr
    real(rk), intent(out) :: u_out(n_dof)

    integer  :: n_aug, info
    real(rk), allocatable :: K_aug(:,:), rhs_aug(:)
    integer,  allocatable :: ipiv(:)

    n_aug = n_dof + n_constr
    allocate(K_aug(n_aug,n_aug), rhs_aug(n_aug), ipiv(n_aug))

    K_aug = 0.0_rk
    K_aug(1:n_dof,       1:n_dof)       = K
    K_aug(n_dof+1:n_aug, 1:n_dof)       = C_mat
    K_aug(1:n_dof,       n_dof+1:n_aug) = transpose(C_mat)

    rhs_aug        = 0.0_rk
    rhs_aug(1:n_dof) = f_rhs

    call dgesv(n_aug, 1, K_aug, n_aug, ipiv, rhs_aug, n_aug, info)
    if (info /= 0) stop 'solve_augmented: DGESV failed'

    u_out = rhs_aug(1:n_dof)
    deallocate(K_aug, rhs_aug, ipiv)
  end subroutine solve_augmented

  !---------------------------------------------------------------------------
  ! average_strain_partition:
  !   Computes volume-averaged strain over partition B.
  !     avg_eps(6) = (1/V_B) * sum_{e in B} integral_e B * u_e dV
  !---------------------------------------------------------------------------
  subroutine average_strain_partition(m, part_id, part_B, u, avg_eps)
    use ebh_mesh, only: MeshData
    type(MeshData), intent(in) :: m
    integer,        intent(in) :: part_id(:), part_B
    real(rk),       intent(in) :: u(:)           ! full displacement vector
    real(rk),       intent(out):: avg_eps(NV)

    real(rk) :: coords_e(8,3), B(NV,NDOF_E), detJ, ue(NDOF_E)
    real(rk) :: eps_gp(NV), xi_g(NG,3), w_g(NG)
    real(rk) :: vol_B
    integer  :: gdofs(NDOF_E), ie, ig, j

    call gauss_pts_hex8(xi_g, w_g)
    avg_eps = 0.0_rk
    vol_B   = 0.0_rk

    do ie = 1, m%n_hex
      if (part_id(ie) /= part_B) cycle
      call get_elem_coords(m, ie, coords_e)
      call get_elem_dofs  (m, ie, gdofs)
      do j = 1, NDOF_E
        ue(j) = u(gdofs(j))
      end do
      do ig = 1, NG
        call Bmat_hex8(xi_g(ig,1),xi_g(ig,2),xi_g(ig,3), coords_e, B, detJ)
        eps_gp  = matmul(B, ue)
        avg_eps = avg_eps + w_g(ig)*detJ * eps_gp
        vol_B   = vol_B   + w_g(ig)*detJ
      end do
    end do

    if (vol_B > 0.0_rk) avg_eps = avg_eps / vol_B
  end subroutine average_strain_partition

end module ebh_assembly

!=============================================================================
! MODULE 8 – I/O utilities
!=============================================================================
module ebh_io
  use ebh_kinds
  implicit none
  private
  public :: write_E_inf, write_P_inf, write_vtk_hex


contains

  !---------------------------------------------------------------------------
  ! Write E_inf to text file
  !  Format: n_parts*NV rows x NV columns
  !  Row  (B-1)*NV + ij  corresponds to  E^{kl=column, B}_ij=row
  !---------------------------------------------------------------------------
  subroutine write_E_inf(filename, E_inf, n_parts)
    character(len=*), intent(in) :: filename
    real(rk),         intent(in) :: E_inf(:,:,:)   ! (NV, NV, n_parts)
    integer,          intent(in) :: n_parts

    integer :: u, ios, B_part, ij, kl

    open(newunit=u, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) stop 'write_E_inf: cannot open file'

    write(u,'(A)') '# E_inf(ij, kl, B)  -- rows: (B-1)*6+ij, cols: kl'
    do B_part = 1, n_parts
      do ij = 1, NV
        write(u,'(6(1X,ES22.14))') (E_inf(ij,kl,B_part), kl=1,NV)
      end do
    end do
    close(u)
    write(*,'(A,A)') '  Written: ', trim(filename)
  end subroutine write_E_inf

  !---------------------------------------------------------------------------
  ! Write P_inf to text file
  !  Format: n_parts*NV rows x n_parts*NV columns
  !  Row  (B-1)*NV+ij,  Col  (A-1)*NV+kl  -->  P^{kl,B,A}_ij
  !---------------------------------------------------------------------------
  subroutine write_P_inf(filename, P_inf, n_parts)
    character(len=*), intent(in) :: filename
    real(rk),         intent(in) :: P_inf(:,:,:,:)  ! (NV,NV,n_parts,n_parts)
    integer,          intent(in) :: n_parts

    integer  :: u, ios, B_part, ij, A, kl
    character(len=32) :: fmt_str

    open(newunit=u, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) stop 'write_P_inf: cannot open file'

    write(fmt_str,'(A,I0,A)') '(', NV*n_parts, '(1X,ES22.14))'
    write(u,'(A)') '# P_inf(ij,kl,B,A) -- rows: (B-1)*6+ij, cols: (A-1)*6+kl'
    do B_part = 1, n_parts
      do ij = 1, NV
        write(u,fmt_str) ((P_inf(ij,kl,B_part,A), kl=1,NV), A=1,n_parts)
      end do
    end do
    close(u)
    write(*,'(A,A)') '  Written: ', trim(filename)
  end subroutine write_P_inf

  !---------------------------------------------------------------------------
  ! Write displacement field as VTK UnstructuredGrid for hex8
  !---------------------------------------------------------------------------
  subroutine write_vtk_hex(filename, coords, hex_conn, u, n_nodes, n_hex)
    character(len=*), intent(in) :: filename
    real(rk),         intent(in) :: coords(n_nodes,3)
    integer,          intent(in) :: hex_conn(n_hex,8)
    real(rk),         intent(in) :: u(n_nodes*3)
    integer,          intent(in) :: n_nodes, n_hex

    integer :: fid, i, ios

    open(newunit=fid, file=trim(filename)//'.vtu', status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) return

    write(fid,'(A)') '<?xml version="1.0"?>'
    write(fid,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1">'
    write(fid,'(A)') '  <UnstructuredGrid>'
    write(fid,'(A,I0,A,I0,A)') &
      '    <Piece NumberOfPoints="',n_nodes,'" NumberOfCells="',n_hex,'">'
    write(fid,'(A)') '      <Points>'
    write(fid,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    do i = 1, n_nodes
      write(fid,'(3(1X,ES16.8))') coords(i,1), coords(i,2), coords(i,3)
    end do
    write(fid,'(A)') '        </DataArray>'
    write(fid,'(A)') '      </Points>'
    write(fid,'(A)') '      <PointData Vectors="Displacement">'
    write(fid,'(A)') '        <DataArray type="Float64" Name="Displacement"'// &
                     ' NumberOfComponents="3" format="ascii">'
    do i = 1, n_nodes
      write(fid,'(3(1X,ES16.8))') u(3*(i-1)+1), u(3*(i-1)+2), u(3*(i-1)+3)
    end do
    write(fid,'(A)') '        </DataArray>'
    write(fid,'(A)') '      </PointData>'
    write(fid,'(A)') '      <Cells>'
    write(fid,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do i = 1, n_hex
      write(fid,'(8(1X,I0))') hex_conn(i,:) - 1   ! VTK is 0-based
    end do
    write(fid,'(A)') '        </DataArray>'
    write(fid,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
    do i = 1, n_hex
      write(fid,'(1X,I0)') 8*i
    end do
    write(fid,'(A)') '        </DataArray>'
    write(fid,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
    do i = 1, n_hex
      write(fid,'(A)') ' 12'   ! VTK_HEXAHEDRON = 12
    end do
    write(fid,'(A)') '        </DataArray>'
    write(fid,'(A)') '      </Cells>'
    write(fid,'(A)') '    </Piece>'
    write(fid,'(A)') '  </UnstructuredGrid>'
    write(fid,'(A)') '</VTKFile>'
    close(fid)
  end subroutine write_vtk_hex

end module ebh_io

!=============================================================================
! MODULE 7 – Influence function computation (the actual offline stage)
!
!  E_inf(NV, NV, n_parts)           E^{kl,B}_{ij}
!  P_inf(NV, NV, n_parts, n_parts)  P^{kl,B,A}_{ij}
!
!  Paper equations (8)-(10):
!    eps^B = E^{kl,B} eps_bar_kl  +  P^{kl,B,A} mu^A_kl
!
!  Two sets of RVE problems:
!   (I)  6 RHS:          unit macro strain, mu=0    -->  E^{kl,B}
!   (II) 6*n_parts RHS:  eps_bar=0, unit mu on A   -->  P^{kl,B,A}
!=============================================================================
module ebh_influence
  use ebh_kinds
  use ebh_io
  use ebh_mesh,     only: MeshData
  use ebh_assembly, only: assemble_K, assemble_F_macro, assemble_F_eigenstrain, &
                           solve_augmented, average_strain_partition
  use ebh_pbc,      only: build_pbc_constraints
  implicit none
  private
  public :: compute_influence_functions


contains

  subroutine compute_influence_functions(m, n_phases, phase_E, phase_nu, &
                                         part_id, n_parts, E_inf, P_inf)
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: n_phases, n_parts
    real(rk),       intent(in)  :: phase_E(n_phases), phase_nu(n_phases)
    integer,        intent(in)  :: part_id(:)    ! (n_hex) partition per element
    real(rk),       allocatable, intent(out) :: E_inf(:,:,:)   ! (NV, NV, n_parts)
    real(rk),       allocatable, intent(out) :: P_inf(:,:,:,:) ! (NV, NV, n_parts, n_parts)

    integer  :: n_dof, n_constr, kl, A, B_part
    real(rk), allocatable :: K(:,:), C_mat(:,:), RHS_c(:)
    real(rk), allocatable :: F_macro(:,:), F_eig(:), u_sol(:)
    integer,  allocatable :: slave_dof(:), master_dof(:), corner_dof(:)
    character(len=64) :: fname
    integer :: i, j
    character(len=128) :: fname2
    n_dof = 3 * m%n_nodes

    allocate(E_inf(NV, NV, n_parts))
    allocate(P_inf(NV, NV, n_parts, n_parts))
    E_inf = 0.0_rk
    P_inf = 0.0_rk

    ! -------- Stiffness matrix (assembled once) --------
    write(*,'(A)') '  Assembling global stiffness K ...'
    allocate(K(n_dof, n_dof))
    call assemble_K(m, n_phases, phase_E, phase_nu, K)

    ! -------- Periodic BC constraints --------
    write(*,'(A)') '  Building periodic BC constraints ...'
    call build_pbc_constraints(m, n_dof, slave_dof, master_dof, &
                                corner_dof, C_mat, RHS_c)
    n_constr = size(C_mat, 1)
    write(*,'(A,I6)') '    n_constraints = ', n_constr

    ! -------- Macro-strain load vectors --------
    write(*,'(A)') '  Assembling macro-strain force vectors ...'
    allocate(F_macro(n_dof, NV))
    call assemble_F_macro(m, n_phases, phase_E, phase_nu, F_macro)
 

    allocate(u_sol(n_dof), F_eig(n_dof))

    ! ================================================================
    ! PROBLEM SET I: 6 macro-strain problems  -->  E^{kl,B}
    ! ================================================================
    write(*,'(A)') '  Solving macro-strain problems (6 RHS) ...'
    do kl = 1, NV
      write(*,'(A,I1,A,I1)') '    kl = ', kl, ' / ', NV
      ! RHS = -F_macro(:,kl)   (minus sign: K u = -F_macro)
      call solve_augmented(K, C_mat, -F_macro(:,kl), n_dof, n_constr, u_sol)
      write(fname,'(A,I1)') 'mode_macro_', kl
      call write_vtk_hex(trim(fname), m%coords, m%hex_conn, u_sol, m%n_nodes, m%n_hex)
      ! Average strain over each partition B
      do B_part = 1, n_parts
        call average_strain_partition(m, part_id, B_part, u_sol, E_inf(:,kl,B_part))
        ! Add identity term:  eps^B = I*eps_bar + (fluctuation)
        ! The fluctuation is what average_strain_partition returns (gradient of u_fluct).
        ! The full E^{kl,B} = delta_{ij,kl} + <dH^kl/dy>_B
        ! i.e. we add the Kronecker delta for the macro-strain component
        E_inf(kl, kl, B_part) = E_inf(kl, kl, B_part) + 1.0_rk
      end do
    end do
    
    ! ================================================================
    ! PROBLEM SET II: 6 * n_parts eigenstrain problems  -->  P^{kl,B,A}
    ! ================================================================
    write(*,'(A,I4,A,I4,A)') '  Solving eigenstrain problems (', &
                               NV*n_parts, ' = ', NV, ' x n_parts RHS) ...'
    do A = 1, n_parts
      write(*,'(A,I4,A,I4)') '    Partition A = ', A, ' / ', n_parts
      do kl = 1, NV
        call assemble_F_eigenstrain(m, n_phases, phase_E, phase_nu, &
                                    part_id, A, kl, F_eig)

        call solve_augmented(K, C_mat, F_eig, n_dof, n_constr, u_sol)

        write(fname2,'(A,I0,A,I0)') 'vtk_outputs/mode_phase_A', A, '_kl', kl
        call write_vtk_hex(trim(fname2), m%coords, m%hex_conn, u_sol, m%n_nodes, m%n_hex)

        do B_part = 1, n_parts
          call average_strain_partition(m, part_id, B_part, u_sol, P_inf(:,kl,B_part,A))
        end do
      end do
    end do

    deallocate(K, C_mat, RHS_c, F_macro, F_eig, u_sol)
    deallocate(slave_dof, master_dof, corner_dof)
    write(*,'(A)') '  Influence functions computed.'
  end subroutine compute_influence_functions

end module ebh_influence

!=============================================================================
! MODULE 9 – Partition builder: strips lineales O sectores angulares
!=============================================================================
module ebh_partitions_strips
  use ebh_kinds
  use ebh_mesh, only: MeshData
  implicit none
  private
  public :: assign_strip_partitions, assign_angular_partitions, &
            assign_radial_angular_partitions, write_partition_file

contains

  !==========================================================================
  ! MODO 1: strips lineales (comportamiento original)
  !==========================================================================
  subroutine assign_strip_partitions(m, matrix_tag, fiber_tag, M_strips, &
                                      strip_axis, part_id, n_parts)
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: matrix_tag, fiber_tag
    integer,        intent(in)  :: M_strips
    integer,        intent(in)  :: strip_axis   ! 1=x, 2=y, 3=z
    integer, allocatable, intent(out) :: part_id(:)
    integer,              intent(out) :: n_parts

    integer  :: ie, n_mat, i, current_strip, count_in_strip
    integer  :: elems_per_strip, remainder, strip_size
    real(rk) :: cx, cy, cz
    integer,  allocatable :: mat_idx(:), sort_ord(:)
    real(rk), allocatable :: mat_coord(:)

    allocate(part_id(m%n_hex))
    part_id = 0
    n_mat = 0

    do ie = 1, m%n_hex
      if     (m%hex_tag(ie) == fiber_tag)  then; part_id(ie) = 1
      else if(m%hex_tag(ie) == matrix_tag) then; n_mat = n_mat + 1
      else
        write(*,'(A,I4,A,I6)') '  WARNING: unknown tag ', m%hex_tag(ie), &
                                ' on element ', ie
      end if
    end do

    if (n_mat   == 0)       stop 'assign_strip_partitions: no matrix elements'
    if (M_strips < 1)       stop 'assign_strip_partitions: M_strips < 1'
    if (M_strips > n_mat)   stop 'assign_strip_partitions: more strips than elements'

    write(*,'(A,I7,A,I4,A,I1)') &
      '  [strip] Matrix elems: ', n_mat, '  M_strips: ', M_strips, &
      '  axis: ', strip_axis

    allocate(mat_idx(n_mat), mat_coord(n_mat), sort_ord(n_mat))
    n_mat = 0
    do ie = 1, m%n_hex
      if (m%hex_tag(ie) /= matrix_tag) cycle
      n_mat = n_mat + 1
      mat_idx(n_mat) = ie
      call elem_centroid(m, ie, cx, cy, cz)
      select case (strip_axis)
        case(1); mat_coord(n_mat) = cx
        case(2); mat_coord(n_mat) = cy
        case(3); mat_coord(n_mat) = cz
        case default; stop 'assign_strip_partitions: strip_axis must be 1,2,3'
      end select
    end do

    do i = 1, n_mat; sort_ord(i) = i; end do
    call argsort_ins(mat_coord, n_mat, sort_ord)

    elems_per_strip = n_mat / M_strips
    remainder       = mod(n_mat, M_strips)
    current_strip   = 1
    count_in_strip  = 0

    do i = 1, n_mat
      ie = mat_idx(sort_ord(i))
      part_id(ie) = current_strip + 1
      count_in_strip = count_in_strip + 1
      if (current_strip <= remainder) then
        strip_size = elems_per_strip + 1
      else
        strip_size = elems_per_strip
      end if
      if (count_in_strip == strip_size .and. current_strip < M_strips) then
        current_strip  = current_strip + 1
        count_in_strip = 0
      end if
    end do

    n_parts = M_strips + 1
    if (any(part_id == 0)) stop 'assign_strip_partitions: unassigned elements'
    call partition_summary(part_id, m%n_hex, n_parts)
    deallocate(mat_idx, mat_coord, sort_ord)
  end subroutine assign_strip_partitions

  !==========================================================================
  ! MODO 2: sectores angulares en plano perpendicular a la fibra
  !==========================================================================
  subroutine assign_angular_partitions(m, matrix_tag, fiber_tag, M_sectors, &
                                        fiber_axis, part_id, n_parts)
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: matrix_tag, fiber_tag
    integer,        intent(in)  :: M_sectors   ! número de sectores angulares
    integer,        intent(in)  :: fiber_axis  ! 1=x, 2=y, 3=z
    integer, allocatable, intent(out) :: part_id(:)
    integer,              intent(out) :: n_parts

    integer  :: ie, isec, n_mat
    real(rk) :: cx, cy, cz, ca, cb   ! ca, cb = coordenadas en plano perp.
    real(rk) :: a0, b0               ! centroide del RVE en plano perpendicular
    real(rk) :: theta, dtheta, tmin, tmax
    real(rk), parameter :: pi = 4.0_rk * atan(1.0_rk)

    if (M_sectors < 1) stop 'assign_angular_partitions: M_sectors < 1'

    ! --- centroide del RVE en el plano perpendicular ---
    call rve_centroid_inplane(m, fiber_axis, a0, b0)

    allocate(part_id(m%n_hex))
    part_id = 0
    n_mat   = 0
    dtheta  = 2.0_rk * pi / real(M_sectors, rk)

    do ie = 1, m%n_hex
      if (m%hex_tag(ie) == fiber_tag) then
        part_id(ie) = 1
        cycle
      else if (m%hex_tag(ie) /= matrix_tag) then
        write(*,'(A,I4,A,I6)') '  WARNING: unknown tag ', m%hex_tag(ie), &
                                ' on element ', ie
        cycle
      end if

      ! centroide del elemento en el plano perpendicular a la fibra
      call elem_centroid(m, ie, cx, cy, cz)
      select case (fiber_axis)
        case(1); ca = cy - b0;  cb = cz - a0   ! plano YZ, ángulo=atan2(cz,cy)
        case(2); ca = cx - a0;  cb = cz - b0   ! plano XZ, ángulo=atan2(cz,cx)
        case(3); ca = cx - a0;  cb = cy - b0   ! plano XY, ángulo=atan2(cy,cx)
        case default; stop 'assign_angular_partitions: fiber_axis must be 1,2,3'
      end select

      ! ángulo en [-pi, pi)
      theta = atan2(cb, ca)

      ! mapear a [0, 2*pi)
      if (theta < 0.0_rk) theta = theta + 2.0_rk * pi

      ! sector: k = floor(theta / dtheta) + 1,  clamp a [1, M_sectors]
      isec = int(theta / dtheta) + 1
      if (isec < 1)         isec = 1
      if (isec > M_sectors) isec = M_sectors

      part_id(ie) = isec + 1   ! +1: partición 1 es fibra
      n_mat = n_mat + 1
    end do

    n_parts = M_sectors + 1

    if (any(part_id == 0)) stop 'assign_angular_partitions: unassigned elements'

    write(*,'(A,I7,A,I4,A,I1)') &
      '  [angular] Matrix elems: ', n_mat, '  M_sectors: ', M_sectors, &
      '  fiber_axis: ', fiber_axis
    call partition_summary(part_id, m%n_hex, n_parts)

  end subroutine assign_angular_partitions

  !==========================================================================
  ! MODO 3: sectores angulares + anillos radiales, igual que RVE_Sol.py
  !
  !  part_id = 1 --> fibra
  !  matrix_part = angular_id*n_radial + radial_id + 2
  !  angular_id = 0..n_angular-1, radial_id = 0..n_radial-1
  !==========================================================================
  subroutine assign_radial_angular_partitions(m, matrix_tag, fiber_tag, &
                                              n_angular, n_radial, fiber_axis, &
                                              part_id, n_parts)
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: matrix_tag, fiber_tag
    integer,        intent(in)  :: n_angular, n_radial
    integer,        intent(in)  :: fiber_axis
    integer, allocatable, intent(out) :: part_id(:)
    integer,              intent(out) :: n_parts

    integer  :: ie, n_mat, angular_id, radial_id
    real(rk) :: cx, cy, cz, ca, cb
    real(rk) :: a0, b0, theta, r, r_inner, r_outer, rho, dtheta
    real(rk) :: half_a, half_b
    real(rk), parameter :: pi = 4.0_rk * atan(1.0_rk)
    real(rk), parameter :: eps = 1.0e-14_rk

    if (n_angular < 1) stop 'assign_radial_angular_partitions: n_angular < 1'
    if (n_radial  < 1) stop 'assign_radial_angular_partitions: n_radial < 1'

    call rve_centroid_inplane(m, fiber_axis, a0, b0)
    half_a = 0.5_rk
    half_b = 0.5_rk
    dtheta = 2.0_rk * pi / real(n_angular, rk)

    allocate(part_id(m%n_hex))
    part_id = 0
    n_mat = 0
    r_inner = huge(1.0_rk)

    ! First pass: assign fiber and find the minimum matrix radius.
    do ie = 1, m%n_hex
      if (m%hex_tag(ie) == fiber_tag) then
        part_id(ie) = 1
        cycle
      else if (m%hex_tag(ie) /= matrix_tag) then
        write(*,'(A,I4,A,I6)') '  WARNING: unknown tag ', m%hex_tag(ie), &
                                ' on element ', ie
        cycle
      end if

      call elem_centroid(m, ie, cx, cy, cz)
      call inplane_coords(cx, cy, cz, a0, b0, fiber_axis, ca, cb)
      r = sqrt(ca*ca + cb*cb)
      r_inner = min(r_inner, r)
      n_mat = n_mat + 1
    end do

    if (n_mat == 0) stop 'assign_radial_angular_partitions: no matrix elements'

    ! Second pass: assign matrix partitions.
    do ie = 1, m%n_hex
      if (m%hex_tag(ie) /= matrix_tag) cycle

      call elem_centroid(m, ie, cx, cy, cz)
      call inplane_coords(cx, cy, cz, a0, b0, fiber_axis, ca, cb)

      theta = atan2(cb, ca)
      if (theta < 0.0_rk) theta = theta + 2.0_rk * pi

      r = sqrt(ca*ca + cb*cb)
      r_outer = square_radius(theta, half_a, half_b)
      rho = (r - r_inner) / (r_outer - r_inner + eps)
      rho = max(0.0_rk, min(1.0_rk, rho))

      angular_id = int(theta / dtheta)
      angular_id = max(0, min(n_angular - 1, angular_id))

      radial_id = int(rho * real(n_radial, rk))
      radial_id = max(0, min(n_radial - 1, radial_id))

      part_id(ie) = angular_id*n_radial + radial_id + 2
    end do

    n_parts = 1 + n_angular*n_radial
    if (any(part_id == 0)) stop 'assign_radial_angular_partitions: unassigned elements'

    write(*,'(A,I7,A,I4,A,I4,A,I1)') &
      '  [radial_angular] Matrix elems: ', n_mat, '  n_angular: ', n_angular, &
      '  n_radial: ', n_radial, '  fiber_axis: ', fiber_axis
    call partition_summary(part_id, m%n_hex, n_parts)
  end subroutine assign_radial_angular_partitions

  !==========================================================================
  ! Helpers
  !==========================================================================
  subroutine rve_centroid_inplane(m, fiber_axis, a0, b0)
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: fiber_axis
    real(rk),       intent(out) :: a0, b0
    real(rk) :: xmin, xmax, ymin, ymax, zmin, zmax

    xmin = minval(m%coords(:,1)); xmax = maxval(m%coords(:,1))
    ymin = minval(m%coords(:,2)); ymax = maxval(m%coords(:,2))
    zmin = minval(m%coords(:,3)); zmax = maxval(m%coords(:,3))

    select case (fiber_axis)
      case(1)
        a0 = 0.5_rk*(zmin + zmax)   ! z
        b0 = 0.5_rk*(ymin + ymax)   ! y
      case(2)
        a0 = 0.5_rk*(xmin + xmax)   ! x
        b0 = 0.5_rk*(zmin + zmax)   ! z
      case(3)
        a0 = 0.5_rk*(xmin + xmax)   ! x
        b0 = 0.5_rk*(ymin + ymax)   ! y
      case default
        stop 'rve_centroid_inplane: fiber_axis must be 1,2,3'
    end select

    write(*,'(A,F8.4,A,F8.4,A,I1)') &
      '  RVE in-plane bbox center: (', a0, ', ', b0, ')  fiber_axis=', fiber_axis
  end subroutine rve_centroid_inplane

  subroutine inplane_coords(cx, cy, cz, a0, b0, fiber_axis, ca, cb)
    real(rk), intent(in)  :: cx, cy, cz, a0, b0
    integer,  intent(in)  :: fiber_axis
    real(rk), intent(out) :: ca, cb

    select case (fiber_axis)
      case(1)
        ca = cy - b0
        cb = cz - a0
      case(2)
        ca = cx - a0
        cb = cz - b0
      case(3)
        ca = cx - a0
        cb = cy - b0
      case default
        stop 'inplane_coords: fiber_axis must be 1,2,3'
    end select
  end subroutine inplane_coords

  real(rk) function square_radius(theta, half_a, half_b)
    real(rk), intent(in) :: theta, half_a, half_b
    real(rk) :: c, ss
    real(rk), parameter :: eps = 1.0e-14_rk

    c = abs(cos(theta))
    ss = abs(sin(theta))
    square_radius = min(half_a / max(c, eps), half_b / max(ss, eps))
  end function square_radius

  !---------------------------------------------------------------------------
  ! Centroid of hex8 element = average of 8 nodal coordinates
  !---------------------------------------------------------------------------
  subroutine elem_centroid(m, ie, cx, cy, cz)
    type(MeshData), intent(in)  :: m
    integer,        intent(in)  :: ie
    real(rk),       intent(out) :: cx, cy, cz
    integer :: a, nd
    cx = 0.0_rk; cy = 0.0_rk; cz = 0.0_rk
    do a = 1, 8
      nd = m%hex_conn(ie, a)
      cx = cx + m%coords(nd,1)
      cy = cy + m%coords(nd,2)
      cz = cz + m%coords(nd,3)
    end do
    cx = cx * 0.125_rk
    cy = cy * 0.125_rk
    cz = cz * 0.125_rk
  end subroutine elem_centroid

  !---------------------------------------------------------------------------
  ! Ascending insertion argsort
  !---------------------------------------------------------------------------
  subroutine argsort_ins(arr, n, idx)
    integer,  intent(in)    :: n
    real(rk), intent(in)    :: arr(n)
    integer,  intent(inout) :: idx(n)
    integer  :: i, j, tmp
    real(rk) :: key
    do i = 2, n
      tmp = idx(i); key = arr(tmp); j = i - 1
      do while (j >= 1)
        if (arr(idx(j)) <= key) exit
        idx(j+1) = idx(j); j = j - 1
      end do
      idx(j+1) = tmp
    end do
  end subroutine argsort_ins

  !---------------------------------------------------------------------------
  ! Print element count per partition
  !---------------------------------------------------------------------------
  subroutine partition_summary(part_id, n_hex, n_parts)
    integer, intent(in) :: part_id(n_hex), n_hex, n_parts
    integer, allocatable :: cnt(:)
    integer :: i
    allocate(cnt(n_parts)); cnt = 0
    do i = 1, n_hex; cnt(part_id(i)) = cnt(part_id(i)) + 1; end do
    write(*,'(A)') '  Partition summary (part | n_elements):'
    write(*,'(A,I6,A)') '    Part  1 (fiber)   : ', cnt(1), ' elems'
    do i = 2, n_parts
      write(*,'(A,I4,A,I6,A)') '    Part ', i, ' (sector): ', cnt(i), ' elems'
    end do
    deallocate(cnt)
  end subroutine partition_summary

  !---------------------------------------------------------------------------
  ! Save partition assignment to text file  (elem_id  part_id)
  !---------------------------------------------------------------------------
  subroutine write_partition_file(filename, part_id, n_hex)
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: part_id(n_hex), n_hex
    integer :: u, ios, i
    open(newunit=u, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) stop 'write_partition_file: cannot open'
    write(u,'(A)') '# elem_id   part_id'
    do i = 1, n_hex
      write(u,'(I7,2X,I5)') i, part_id(i)
    end do
    close(u)
    write(*,'(A,A)') '  Partition file written: ', trim(filename)
  end subroutine write_partition_file

end module ebh_partitions_strips

!=============================================================================
! MAIN PROGRAM
!=============================================================================
program ebh_offline_main
  use ebh_kinds
  use ebh_hex8
  use ebh_mesh,              only: MeshData, read_gmsh2
  use ebh_partitions_strips, only: assign_strip_partitions, &
                                    assign_angular_partitions, &
                                    assign_radial_angular_partitions, &
                                    write_partition_file
  use ebh_influence,         only: compute_influence_functions
  use ebh_io,                only: write_E_inf, write_P_inf, write_vtk_hex
  implicit none

  !==========================================================================
  ! USER PARAMETERS
  !==========================================================================

  ! -- Files --
  character(len=256), parameter :: MESH_FILE   = 'modelo_rve.msh'
  character(len=256), parameter :: E_OUT_FILE  = 'E_inf.dat'
  character(len=256), parameter :: P_OUT_FILE  = 'P_inf.dat'
  character(len=256), parameter :: PART_FILE   = 'partitions.dat'  ! set '' to skip

  ! -- Material properties (indexed by Gmsh physical volume tag) --
  !    tag 1 = matrix,  tag 2 = fiber  (change if your Gmsh tags differ)
  integer,  parameter :: N_PHASES   = 2
  integer,  parameter :: MATRIX_TAG = 1
  integer,  parameter :: FIBER_TAG  = 2

  real(rk), parameter :: PHASE_E (N_PHASES) = [   100.0_rk,  &  ! matrix [MPa]
                                                10000.0_rk  ]  ! fiber  [MPa]
  real(rk), parameter :: PHASE_NU(N_PHASES) = [   0.3_rk,  &
                                                   0.3_rk  ]
  integer :: i

  ! -- Partitioning --
  !
  !  PART_TYPE = 'strip'
  !    Strips lineales a lo largo de STRIP_AXIS.
  !    Para fibra en Z: usar STRIP_AXIS=3 (planos XY, perpendiculares a fibra).
  !    M_PARTS = número de strips (paper usa 1, 4, 8, 16, 32, 96)
  !
  !  PART_TYPE = 'angular'
  !    Sectores angulares en plano perpendicular a la fibra.
  !    Para fibra en Z: FIBER_AXIS=3, sectores en plano XY.
  !    Respeta la simetría cilíndrica → E_inf isótropo transversal.
  !    M_PARTS = número de sectores angulares (recomendado: 4, 8, 16)
  !
  character(len=16), parameter :: PART_TYPE  = 'radial_angular'
  integer,           parameter :: N_ANGULAR  = 2
  integer,           parameter :: N_RADIAL   = 1
  integer,           parameter :: M_PARTS    = N_ANGULAR * N_RADIAL
  integer,           parameter :: STRIP_AXIS = 3
  integer,           parameter :: FIBER_AXIS = 3

  !==========================================================================

  type(MeshData) :: m
  real(rk), allocatable :: E_inf(:,:,:), P_inf(:,:,:,:)
  integer,  allocatable :: part_id(:)
  integer  :: n_parts
  integer  :: wall_start, wall_end, wall_rate
  real(rk) :: wall_total
  real(rk), parameter :: pi = 4.0d0*atan(1.0d0)
  real(rk) :: t0, t1, L1(6,6), L2(6,6), LIJKL(6,6)
  real(rk) :: A1,a2

  call system_clock(wall_start, wall_rate)
  write(*,*) pi
  ! A1 = pi*0.29152828269495934**2
  A1 = pi*0.35**2
  A2 = 1 - A1

  write(*,*)
  write(*,'(A)') '================================================='
  write(*,'(A)') '  EBH Offline  –  3D Hex8 FEM'
  write(*,'(A)') '  Fish & Cui (2025): Eigenstate-Based Homogenization'
  write(*,'(A)') '================================================='
  write(*,*)

  ! -------------------------------------------------------------------
  ! 1. Read mesh (Gmsh format 2 ASCII)
  ! -------------------------------------------------------------------
  write(*,'(A,A)') 'Reading mesh: ', trim(MESH_FILE)
  call read_gmsh2(MESH_FILE, m)
  write(*,*)

  ! -------------------------------------------------------------------
  ! 2. Assign partitions:
  !      part 1         = all fiber elements (single partition)
  !      parts 2..M+1   = M_STRIPS matrix strips along STRIP_AXIS
  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------
  ! 2. Assign partitions
  ! -------------------------------------------------------------------
  if (trim(PART_TYPE) == 'radial_angular') then
    write(*,'(A,I4,A,I4,A,I1,A)') 'Assigning radial-angular partitions: ', &
      N_ANGULAR, ' angular x ', N_RADIAL, ' radial (fiber_axis=', &
      FIBER_AXIS, ') + 1 fiber partition ...'
    call assign_radial_angular_partitions(m, MATRIX_TAG, FIBER_TAG, &
                                          N_ANGULAR, N_RADIAL, FIBER_AXIS, &
                                          part_id, n_parts)
  else if (trim(PART_TYPE) == 'angular') then
    write(*,'(A,I4,A,I1,A)') 'Assigning ', M_PARTS, &
      ' angular sectors (fiber_axis=', FIBER_AXIS, ') + 1 fiber partition ...'
    call assign_angular_partitions(m, MATRIX_TAG, FIBER_TAG, M_PARTS, &
                                    FIBER_AXIS, part_id, n_parts)
  else
    write(*,'(A,I4,A)') 'Assigning ', M_PARTS, &
      ' linear strips + 1 fiber partition ...'
    call assign_strip_partitions(m, MATRIX_TAG, FIBER_TAG, M_PARTS, &
                                  STRIP_AXIS, part_id, n_parts)
  end if


  write(*,'(A,I4,A)') '  Total partitions: ', n_parts, ' (= M_PARTS + 1 fiber)'
  write(*,*)

  ! Optionally save partition file for inspection / post-processing in Paraview
  if (len_trim(PART_FILE) > 0) &
    call write_partition_file(PART_FILE, part_id, m%n_hex)
  write(*,*)

  ! -------------------------------------------------------------------
  ! 3. Compute influence functions (offline stage)
  !
  !    Solves  6 + 6*n_parts  linear systems:
  !      Set I  (6 RHS)          --> E^{kl,B}_ij
  !      Set II (6*n_parts RHS)  --> P^{kl,B,A}_ij
  ! -------------------------------------------------------------------
  call cpu_time(t0)
  write(*,'(A)') 'Computing influence functions ...'
  write(*,'(A,I5,A)') '  Total linear solves: ', 6 + 6*n_parts, &
                       '  (6 macro + 6 x n_parts eigenstrain)'
  write(*,*)
  call compute_influence_functions(m, N_PHASES, PHASE_E, PHASE_NU, &
                                   part_id, n_parts, E_inf, P_inf)

  call cpu_time(t1)


  call isotropic_D(phase_E(1), phase_nu(1), L2)
  call isotropic_D(phase_E(2), phase_nu(2), L1)

  LIJKL = MATMUL(L2,E_inf(:,:,2))*A2 + MATMUL(L1,E_inf(:,:,1))*A1 !!! matriz y luego fibra 

    write(*,*)LIJKL(1,:)
  write(*,*)LIJKL(2,:)
  write(*,*)LIJKL(3,:)
  write(*,*)LIJKL(4,:)
  write(*,*)LIJKL(5,:)
  write(*,*)LIJKL(6,:)

  write(*,*)
  write(*,'(A,F8.2,A)') 'Offline stage done.  CPU time: ', t1-t0, ' s'
  write(*,*)

  ! -------------------------------------------------------------------
  ! 4. Write output
  !
  !    E_inf.dat : rows (B-1)*6+ij, cols kl  --> E^{kl,B}_{ij}
  !    P_inf.dat : rows (B-1)*6+ij, cols (A-1)*6+kl  --> P^{kl,B,A}_{ij}
  !
  !    These are the only two arrays needed by the EBH online stage.
  ! -------------------------------------------------------------------
  write(*,'(A)') 'Writing output files ...'
  call write_E_inf(E_OUT_FILE, E_inf, n_parts)
  call write_P_inf(P_OUT_FILE, P_inf, n_parts)

  call system_clock(wall_end)
  wall_total = real(wall_end - wall_start, rk) / real(wall_rate, rk)

  write(*,*)
  write(*,'(A)') '================================================='
  write(*,'(A,F10.3,A)') '  Total wall time: ', wall_total, ' s'
  write(*,'(A,F10.3,A)') '  Influence CPU time: ', t1-t0, ' s'
  write(*,'(A,I4)')     '  n_parts:    ', n_parts
  write(*,'(A,A)')      '  PART_TYPE:  ', trim(PART_TYPE)
  write(*,'(A,I7)')     '  n_nodes:    ', m%n_nodes
  write(*,'(A,I7)')     '  n_hex:      ', m%n_hex
  write(*,'(A)') '  Output: E_inf.dat   P_inf.dat'
  write(*,'(A)') '================================================='
  write(*,*)

end program ebh_offline_main
