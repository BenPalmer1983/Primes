! gfortran -Ofast -march=native -fopenmp -o prime-bitarray-parallel prime-bitarray-parallel.f08

! export OMP_NUM_THREADS=4
! export OMP_WAIT_POLICY=active
! export OMP_PROC_BIND=true

! export OMP_NUM_THREADS=8 && export OMP_WAIT_POLICY=active && export OMP_PROC_BIND=true && ./prime-bitarray-parallel
! export OMP_NUM_THREADS=10 && export OMP_WAIT_POLICY=active && export OMP_PROC_BIND=true && ./prime-bitarray-parallel

! Original code from solution by Thomas Jollans
!    https://github.com/PlummersSoftwareLLC/Primes/blob/drag-race/PrimeFortran/solution_2/prime-bitarray.f08
! Added static array for smaller sieves - wastes memory, but I did get a speed increase
! Added openmp to process sieves in parallel

module PrimeModule
    use iso_fortran_env
    implicit none

    integer(kind=int64), dimension(9), parameter :: validated_sieve_sizes = (/ &
        10, &
        100, &
        1000, &
        10000, &
        100000, &
        1000000, &
        10000000, &
        100000000, &
        1000000000 &
    /)
    integer(kind=int64), dimension(9), parameter :: valid_prime_counts = (/ &
        4, &
        25, &
        168, &
        1229, &
        9592, &
        78498, &
        664579, &
        5761455, &
        50847534 &
    /)

    type PrimeSieve
        private
        integer(kind=int8), dimension(:), allocatable :: raw_bits
        integer(kind=int8) :: raw_bits_stack(1:1000001)
        integer(kind=int64) :: sieve_size
        logical :: force_heap = .false.
        logical :: use_heap = .false.
    contains
        procedure, public :: initialize => primesieve_initialize
        procedure, public :: run_sieve => primesieve_run_sieve
        procedure, public :: validate_results => primesieve_validate_results
        procedure, public :: print_results => primesieve_print_results
        procedure, public :: count_primes => primesieve_count_primes
        procedure, private :: get_bit => primesieve_get_bit
        procedure, private :: clear_bits => primesieve_clear_bits
        procedure, private :: get_bit_stack => primesieve_get_bit_stack
        procedure, private :: clear_bits_stack => primesieve_clear_bits_stack
        final :: primesieve_destructor
    end type
    
contains

    subroutine primesieve_initialize(this, sieve_size)
        class(PrimeSieve), intent(inout) :: this
        integer(kind=int64), intent(in) :: sieve_size

        this%sieve_size = sieve_size
        !this%force_heap = .true.
        
        if(this%sieve_size > 16000000 .or. this%force_heap)then
          this%use_heap = .true.
        end if

        if(this%use_heap)then
          allocate(this%raw_bits(sieve_size / 16 + 1))
          this%raw_bits = -1_int8
        else
          this%raw_bits_stack(1:sieve_size / 16 + 1) = -1_int8
        end if
    end subroutine

    subroutine primesieve_destructor(this)
        type(PrimeSieve), intent(inout) :: this

        if(this%use_heap)then
          deallocate(this%raw_bits)
        end if
    end subroutine

    subroutine primesieve_run_sieve(this)
        class(PrimeSieve), intent(inout) :: this

        if(this%use_heap)then
          call primesieve_run_sieve_heap(this)
        else
          call primesieve_run_sieve_stack(this)
        end if
    end subroutine

    subroutine primesieve_run_sieve_heap(this)
        class(PrimeSieve), intent(inout) :: this
        integer(kind=int64) :: factor, q, num

        factor = 3
        q = int(sqrt(real(this%sieve_size)))

        do while (factor <= q)
            do num = factor, this%sieve_size - 1, 2
                if (this%get_bit(num)) then
                    factor = num
                    exit
                end if
            end do

            call this%clear_bits(factor**2, this%sieve_size - 1, factor * 2)

            factor = factor + 2
        end do
    end subroutine

    subroutine primesieve_run_sieve_stack(this)
        class(PrimeSieve), intent(inout) :: this
        integer(kind=int64) :: factor, q, num

        factor = 3
        q = int(sqrt(real(this%sieve_size)))

        do while (factor <= q)
            do num = factor, this%sieve_size - 1, 2
                if (this%get_bit_stack(num)) then
                    factor = num
                    exit
                end if
            end do

            call this%clear_bits_stack(factor**2, this%sieve_size - 1, factor * 2)

            factor = factor + 2
        end do
    end subroutine

    function primesieve_get_bit(this, num) result(bit)
        class(PrimeSieve), intent(in) :: this
        integer(kind=int64), intent(in) :: num
        logical :: bit

        if (iand(num, 1_int64) == 1_int64) then
            ! odd number
            bit = btest(this%raw_bits(num/16 + 1), iand(num/2, 7_int64))
        else
            ! even number
            bit = .false.
        end if
    end function


    function primesieve_get_bit_stack(this, num) result(bit)
        class(PrimeSieve), intent(in) :: this
        integer(kind=int64), intent(in) :: num
        logical :: bit

        if (iand(num, 1_int64) == 1_int64) then
            ! odd number
            bit = btest(this%raw_bits_stack(num/16 + 1), iand(num/2, 7_int64))
        else
            ! even number
            bit = .false.
        end if
    end function

    subroutine primesieve_clear_bits(this, first, last, step)
        class(PrimeSieve), intent(inout) :: this
        integer(kind=int64), intent(in) :: first, last, step
        integer(kind=int64) :: bitidx, idx, bitstep, last_bitidx
        integer(kind=int8) :: bitmask

        ! There is one bit per two natural numbers as we're not keeping track
        ! of the evens.
        bitidx = first/2
        bitstep = step/2
        last_bitidx = last/2

        ! Rather than calculating the bit mask separately each time, we can just
        ! shift it on every iteration
        bitmask = not(int(lshift(1, iand(bitidx, 7_int64)), kind=int8))

        do while (bitidx <= last_bitidx)
            idx = rshift(bitidx, 3) + 1 ! bitidx / 8 + 1
            
            this%raw_bits(idx) = iand(this%raw_bits(idx), bitmask)

            bitmask = ishftc(bitmask, bitstep)
            bitidx = bitidx + bitstep
        end do
    end subroutine

    subroutine primesieve_clear_bits_stack(this, first, last, step)
        class(PrimeSieve), intent(inout) :: this
        integer(kind=int64), intent(in) :: first, last, step
        integer(kind=int64) :: bitidx, idx, bitstep, last_bitidx
        integer(kind=int8) :: bitmask

        ! There is one bit per two natural numbers as we're not keeping track
        ! of the evens.
        bitidx = first/2
        bitstep = step/2
        last_bitidx = last/2

        ! Rather than calculating the bit mask separately each time, we can just
        ! shift it on every iteration
        bitmask = not(int(lshift(1, iand(bitidx, 7_int64)), kind=int8))

        do while (bitidx <= last_bitidx)
            idx = rshift(bitidx, 3) + 1 ! bitidx / 8 + 1
            
            this%raw_bits_stack(idx) = iand(this%raw_bits_stack(idx), bitmask)

            bitmask = ishftc(bitmask, bitstep)
            bitidx = bitidx + bitstep
        end do
    end subroutine

    function primesieve_validate_results(this) result(is_valid)
        class(PrimeSieve), intent(in) :: this
        logical :: is_valid
        integer(kind=int64) :: count, i

        is_valid = .false.
        count = this%count_primes()  
        
        do i = 1, size(valid_prime_counts)
            if (validated_sieve_sizes(i) == this%sieve_size) then
                is_valid = (count == valid_prime_counts(i))
            end if
        end do
    end function



    subroutine primesieve_print_results(this)
        class(PrimeSieve), intent(in) :: this
        integer(kind=int64) :: i

        if (this%sieve_size < 2) then
            return
        end if

        print *, 2
        
        do i = 3, this%sieve_size -1
            if (this%get_bit(i)) then
                print *, i
            end if
        end do
    end subroutine

    function primesieve_count_primes(this) result(count)
        class(PrimeSieve), intent(in) :: this
        integer(kind=int64) :: count, i

        if (this%sieve_size < 2) then
            count = 0
            return
        end if

        count = 1

        if(this%use_heap)then
          do i = 3, this%sieve_size -1
            if (this%get_bit(i)) then
              count = count + 1
            end if
          end do
        else         
          do i = 3, this%sieve_size -1
            if (this%get_bit_stack(i)) then
              count = count + 1
            end if
          end do
        end if
    end function

end module

program PrimeFortran
    use iso_fortran_env
    use PrimeModule
    USE omp_lib
    implicit none

    integer(kind=int64), parameter :: sieve_size = 1000000
    integer(kind=int64), parameter :: benchmark_secs = 5

    integer(kind=int64) :: start_clock, end_clock, cur_clock, &
                       clock_count_rate, clock_count_max
    integer(kind=int64) :: iters = 0
    integer(kind=int64) :: t_iters = 0
    integer(kind=int64) :: tn, tcount
    real :: time_elapsed

    call system_clock(start_clock, clock_count_rate, clock_count_max)
    end_clock = start_clock + benchmark_secs * clock_count_rate
    cur_clock = start_clock


    do while (cur_clock < end_clock)
      !$OMP PARALLEL SHARED(tcount, iters) PRIVATE(t_iters)
      tcount = OMP_GET_NUM_THREADS()
      !$OMP DO
      do tn = 1, tcount
        call run_prime_sieve()
        t_iters = t_iters + 1
      end do
      !$OMP END DO
      !$OMP CRITICAL
      iters = iters + t_iters
      !$OMP END CRITICAL
      !$OMP END PARALLEL
      call system_clock(cur_clock)
    end do

    time_elapsed = real(cur_clock - start_clock) / real(clock_count_rate)

    CALL results(sieve_size, iters, time_elapsed)

contains

    subroutine run_prime_sieve()
        type(PrimeSieve) :: sieve

        call sieve%initialize(sieve_size)
        call sieve%run_sieve()
    end subroutine


    subroutine results(sieve_size, iters, time_elapsed)
        integer(kind=int64), INTENT(in) :: sieve_size 
        integer(kind=int64), INTENT(in) :: iters
        real, INTENT(in) :: time_elapsed

        logical :: valid
        type(PrimeSieve) :: sieve

        call sieve%initialize(sieve_size)
        call sieve%run_sieve()
        valid = sieve%validate_results()
        if(valid)then
          !write(*,"(A)") "Valid"
        else
          write(*,"(A)") "Not Valid"
        end if

        write (*, "(A,I0,A,F0.8,A,I0,A)") &
            "benpalmer1983-bits-par;", iters, ";", time_elapsed, ";", tcount, &
            ";algorithm=base,faithful=yes,bits=1"
    end subroutine

end program
