module libnpcf
    use iso_c_binding
    
    private
    public :: npcf
    
    include "ganpcf_cdef.f90"
    
    type npcf
        private
        type(c_ptr) :: ptr
        
        contains
#ifdef __GNUC__
        procedure :: delete => delete_npcf_polymorph
#else
        final :: delete_npcf
#endif
        procedure :: getShells => get_shells
        procedure :: getNumTriangles => get_num_triangles
        procedure :: getTriangles => get_triangles
        procedure :: setNumParticles => set_num_particles
        procedure :: get2ptSize => get_2pt_size
        procedure :: get3ptSize => get_3pt_size
        procedure :: calculateCorrelations => calculate_correlations
        procedure :: get2pt => get_2pt
        procedure :: get3pt => get_3pt
    end type
    
    interface npcf
        procedure create_npcf
    end interface
    
    contains
        function create_npcf(timesRans, numShells, volBox, rMin, rMax)
            implicit none
            type(npcf) :: create_npcf
            integer, intent(in) :: timesRans
            integer, intent(in) :: numShells
            double precision, intent(in) :: volBox
            double precision, intent(in) :: rMin
            double precision, intent(in) :: rMax
            create_npcf%ptr = create_npcf_c(timesRans, numShells, volBox, rMin, rMax)
        end function
        
        subroutine delete_npcf(this)
            implicit none
            type(npcf) :: this
            call delete_npcf_c(this%ptr)
        end subroutine
        
        subroutine delete_npcf_polymorph(this)
            implicit none
            class(npcf) :: this
            call delete_npcf_c(this%ptr)
        end subroutine
        
        integer function get_shells(this, shells)
            implicit none
            class(npcf) :: this
            double precision, dimension(:) :: shells
            get_shells = get_shells_c(this%ptr, shells)
        end function
        
        integer function get_num_triangles(this)
            implicit none
            class(npcf) :: this
            get_num_triangles = get_num_triangles_c(this%ptr)
        end function
        
        integer function get_triangles(this, tris)
            implicit none
            class(npcf) :: this
            type(float3), dimension(:) :: tris
            get_triangles = get_triangles_c(this%ptr, tris)
        end function
        
        integer function set_num_particles(this, numParts)
            implicit none
            class(npcf) :: this
            integer, intent(in) :: numParts
            set_num_particles = set_num_particles_c(this%ptr, numParts)
        end function
        
        integer function get_2pt_size(this)
            implicit none
            class(npcf) :: this
            get_2pt_size = get_2pt_size_c(this%ptr)
        end function
        
        integer function get_3pt_size(this)
            implicit none
            class(npcf) :: this
            get_3pt_size = get_3pt_size_c(this%ptr)
        end function
        
        integer function calculate_correlations(this, galaxies)
            implicit none
            class(npcf) :: this
            type(float3), dimension(:) :: galaxies
            calculate_correlations = calculate_correlations_c(this%ptr, galaxies)
        end function
        
        integer function get_2pt(this, twoPoint)
            implicit none
            class(npcf) :: this
            double precision, dimension(:) :: twoPoint
            get_2pt = get_2pt_c(this%ptr, twoPoint)
        end function
        
        integer function get_3pt(this, threePoint)
            implicit none
            class(npcf) :: this
            double precision, dimension(:) :: threePoint
            get_3pt = get_3pt_c(this%ptr, threePoint)
        end function
        
end module
        
