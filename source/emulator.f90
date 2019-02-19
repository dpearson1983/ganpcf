program emulator
    
    type, bind(C) :: float3
        real(kind = 4) :: x, y, z
    end type
    
    type(float3), dimension(:), allocatable :: pos
    integer :: s1, i
    
    open(1, file = 'LNKNLogsVelFortran_01.dat', status = 'old')
    read (1,*) s1
    allocate(pos(s1))
    
    do i = 1,s1
        read(1,*)pos(i)
    end do
    close(1)
    
    do i = 1,10
        print *, pos(i)
    end do
    
end program emulator
