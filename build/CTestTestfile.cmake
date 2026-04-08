# CMake generated Testfile for 
# Source directory: C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper
# Build directory: C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test(wigner_compat_test "C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/build/Debug/wigner_compat_test.exe")
  set_tests_properties(wigner_compat_test PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/CMakeLists.txt;86;add_test;C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test(wigner_compat_test "C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/build/Release/wigner_compat_test.exe")
  set_tests_properties(wigner_compat_test PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/CMakeLists.txt;86;add_test;C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test(wigner_compat_test "C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/build/MinSizeRel/wigner_compat_test.exe")
  set_tests_properties(wigner_compat_test PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/CMakeLists.txt;86;add_test;C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test(wigner_compat_test "C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/build/RelWithDebInfo/wigner_compat_test.exe")
  set_tests_properties(wigner_compat_test PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/CMakeLists.txt;86;add_test;C:/Users/ppl/Documents/Analysis/PollenPhaseTransitionPaper/CMakeLists.txt;0;")
else()
  add_test(wigner_compat_test NOT_AVAILABLE)
endif()
