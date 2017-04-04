add_executable(${PROJECT_NAME}
  iniflw.f 
  flow.f
  oscarnxtev.f 
  bms_ap.f 
  obstring.f 
  cmshif.f 
  angle.f
  input.f 
  gridin.f 
  output.f 
  obsvalue.f 
  evalue.f
  normvalue.f
  cubic.f 
  ratio.f 
  jacobi.f 
  eigsrt.f 
  sqr.f 
  reftrans.f
  readwt.f
  )

# install compiled binary
install (TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
