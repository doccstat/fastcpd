module cfuncs
   interface
      subroutine cnlprt(msg, plen) bind(C, name = 'cnlprt_C')
         use iso_c_binding
         character(kind = c_char) :: msg(*)
         integer(kind = c_int)    :: plen
      end subroutine cnlprt
      subroutine h100s(i1, i2, d1, d2, d3, d4, a1, a2, d5) bind(C, name = 'h100s_C')
         use iso_c_binding
         character(kind = c_char) :: a1, a2
         real(kind = c_double)    :: d1, d2, d3, d4, d5
         integer(kind = c_int)    :: i1, i2
      end subroutine h100s
      subroutine h100l(i1, i2, d1, d2, d3, d4, a1, a2, d5, d6, d7) bind(C, name = 'h100l_C')
         use iso_c_binding
         character(kind = c_char) :: a1, a2
         real(kind = c_double)    :: d1, d2, d3, d4, d5, d6, d7
         integer(kind = c_int)    :: i1, i2
      end subroutine h100l
   end interface
end module cfuncs
