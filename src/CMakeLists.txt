cmake_minimum_required (VERSION 2.6)
project (bmsap)

# load Fortran compiler
enable_language (Fortran)

set (CMAKE_Fortran_FLAGS "-ffixed-line-length-none")

add_executable(${PROJECT_NAME}
   angle.f
   cmshif.f
   eigsrt.f
   flow.f
   gridin.f
   iniflw.f
   jacobi.f
   obsvalue.f
   oscarnxtev.f
   reftrans.f
   sqr.f
   bms_ap.f
   cubic.f
   evalue.f
   input.f
   normvalue.f
   obstring.f
   output.f
   ratio.f
   readwt.f
  )

# install compiled binary
# install (TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
install (TARGETS ${PROJECT_NAME} DESTINATION bin)

