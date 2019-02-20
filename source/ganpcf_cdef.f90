! C function declarations

type, bind(C) :: float3
    real(kind = 4) :: x, y, z
end type

interface
    function create_npcf_c(timesRans, numShells, volBox, rMin, rMax) bind(C, name="create_npcf")
        use iso_c_binding
        implicit none
        type(c_ptr) :: create_npcf_c
        integer(c_int), value :: timesRans
        integer(c_int), value :: numShells
        real(c_double), value :: volBox
        real(c_double), value :: rMax
        real(c_double), value :: rMin
    end function
    
    subroutine delete_npcf_c(obj) bind(C, name="delete_npcf")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: obj
    end subroutine
    
    function get_shells_c(obj, shells) bind(C, name="get_shells")
        use iso_c_binding
        implicit none
        integer(c_int) :: get_shells_c
        type(c_ptr), value :: obj
        real(c_double), dimension(:) :: shells
    end function
    
    function get_num_triangles_c(obj) bind(C, name="get_num_triangles")
        use iso_c_binding
        implicit none
        integer(c_int) :: get_num_triangles_c
        type(c_ptr), value :: obj
    end function
    
    function get_triangles_c(obj, tris) bind(C, name="get_triangles")
        use iso_c_binding
        import :: float3
        integer(c_int) :: get_triangles_c
        type(c_ptr), value :: obj
        type(float3), dimension(:) :: tris
    end function
    
    function set_num_particles_c(obj, numParts) bind(C, name="set_num_particles")
        use iso_c_binding
        implicit none
        integer(c_int) :: set_num_particles_c
        type(c_ptr), value :: obj
        integer(c_int), value :: numParts
    end function
    
    function get_2pt_size_c(obj) bind(C, name="get_2pt_size")
        use iso_c_binding
        implicit none
        integer(c_int) :: get_2pt_size_c
        type(c_ptr), value :: obj
    end function
    
    function get_3pt_size_c(obj) bind(C, name="get_3pt_size")
        use iso_c_binding
        implicit none
        integer(c_int) :: get_3pt_size_c
        type(c_ptr), value :: obj
    end function
    
    function calculate_correlations_c(obj, galaxies) bind(C, name="calculate_correlations")
        use iso_c_binding
        import :: float3
        integer(c_int) :: calculate_correlations_c
        type(c_ptr), value :: obj
        type(float3), dimension(:) :: galaxies
    end function
    
    function get_2pt_c(obj, twoPoint) bind(C, name="get_2pt")
        use iso_c_binding
        implicit none
        integer(c_int) :: get_2pt_c
        type(c_ptr), value :: obj
        real(c_double), dimension(:) :: twoPoint
    end function
    
    function get_3pt_c(obj, threePoint) bind(C, name="get_3pt")
        use iso_c_binding
        implicit none
        integer(c_int) :: get_3pt_c
        type(c_ptr), value :: obj
        real(c_double), dimension(:) :: threePoint
    end function
end interface
