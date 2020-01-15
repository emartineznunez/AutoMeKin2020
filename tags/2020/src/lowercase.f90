function string_tolower( string ) result (new) 
    character(len=*)           :: string 

    character(len=len(string)) :: new 

    integer                    :: i 
    integer                    :: k 

    length = len(string) 
    new    = string 
    do i = 1,len(string) 
        k = iachar(string(i:i)) 
        if ( k >= iachar('A') .and. k <= iachar('Z') ) then 
            k = k + iachar('a') - iachar('A') 
            new(i:i) = achar(k) 
        endif 
    enddo 
end function string_tolower
