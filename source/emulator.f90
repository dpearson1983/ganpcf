program emulator
    use libnpcf
    
    type, bind(C) :: float3
        real(kind = 4) :: x, y, z
    end type
    
    type(npcf) :: corr
    type(float3), dimension(:), allocatable :: pos, triangles
    integer :: s1, i, timesRans, numShells, numTriangles, dummy
    double precision :: r_max, r_min, V_box
    double precision, dimension (:), allocatable :: threePoint, twoPoint, shells
    
    ! Information for setting up the correlation function object
    timesRans = 10
    numShells = 32
    V_box = 2500.0**3
    r_max = 32.0
    r_min = 0.0
    
    ! Initialize the library object
    corr = npcf(timesRans, numShells, V_box, r_min, r_max)
    
    ! This function determines the number of unique triangles such that r1 <= r2 <= r3. This will remain
    ! the same for all realizations
    numTriangles = corr%getNumTriangles()
    
    ! Set up storage for the correlation functions
    allocate(shells(numShells))
    allocate(twoPoint(numShells))
    allocate(triangles(numTriangles))
    allocate(threePoint(numTriangles))
    
    ! These two functions return the central radii for the spherical shells for the 2-point function
    ! and the central side lengths for the triangles of the 3-point function
    dummy = corr%getShells(shells)
    dummy = corr%getTriangles(triangles)
    
    ! This part of the code emulates what you will need to do in the looping structure of the MCMC chain
        
        ! Instead of reading from a file, your code will generate the list of galaxies
        open(1, file = 'central_output.txt', status = 'old')
        read (1,*) s1
        allocate(pos(s1)) 
        
        do i = 1,s1
            read(1,*)pos(i)
        end do
        close(1)
        
        ! You need to tell the npcf library how many galaxies there are in a particular realization
        dummy = corr%setNumParticles(s1)
        
        ! Next we tell the library to compute the correlations given the list of galaxies
        dummy = corr%calculateCorrelations(pos)
        
        ! Then we can copy the results back to the Fortran side
        dummy = corr%get2pt(twoPoint)
        dummy = corr%get3pt(threePoint)
        
        ! Your code would then compute the likelihood and move to the next parameter vector guess
        
    ! Output the 2 and 3 point functions for testing purposes
    open(2, file = 'twoPoint.dat', action = 'write', status = 'replace')
    do i = 1,numshells
        write (2,*) shells(i), twoPoint(i)
    end do
    close(2)
    
    open(3, file = 'threePoint.dat', action = 'write', status = 'replace')
    do i = 1,numTriangles
        write (3,*) triangles(i), threePoint(i)
    end do
    close(3)
    
    ! Free the memory of dynamically allocated arrays
    deallocate(shells)
    deallocate(triangles)
    deallocate(twoPoint)
    deallocate(threePoint)
    deallocate(pos)
    
#ifdef __GNUC__
    call corr%delete
#endif
    
end program emulator
