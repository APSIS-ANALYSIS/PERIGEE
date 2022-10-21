# Load appropriate system files into the configuration file

if( $ENV{MACHINE_NAME} MATCHES "poincare")
  message(STATUS "JL's Mac laptop banach")
  include(${CMAKE_CURRENT_LIST_DIR}/poincare.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "kolmogorov")
  message(STATUS "Mac laptop Kolmogorov")
  include(${CMAKE_CURRENT_LIST_DIR}/kolmogorov.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "ladyzhenskaya")
  message(STATUS "Linux desktop Ladyzhenskaya")
  include(${CMAKE_CURRENT_LIST_DIR}/ladyzhenskaya.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "taiyi")
  message(STATUS "Tai-Yi")
  include(${CMAKE_CURRENT_LIST_DIR}/taiyi.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "qiming")
  message(STATUS "Qi-Ming")
  include(${CMAKE_CURRENT_LIST_DIR}/qiming.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "sus_jiashen")
  message(STATUS "jiashen's Taiyi")
  include(${CMAKE_CURRENT_LIST_DIR}/sus_jiashen.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "jiashen_linux")
  message(STATUS "Jiashen's Linux")
  include(${CMAKE_CURRENT_LIST_DIR}/jiashen_linux.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "Sun")
  message(STATUS "Sun's Linux")
  include(${CMAKE_CURRENT_LIST_DIR}/Sun-thinkpad.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "Yujie_ThinkStationP920")
  message(STATUS "Sun's ThinkStationP920")
  include(${CMAKE_CURRENT_LIST_DIR}/ThinkStationP920_Sun.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "ty-sun")
  message(STATUS "Sun's Taiyi")
  include(${CMAKE_CURRENT_LIST_DIR}/Sun_taiyi.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "tianhe2a")
  message(STATUS "TianHe2A")
  include(${CMAKE_CURRENT_LIST_DIR}/Tianhe2A.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "lqs-P920")
  message(STATUS "Luqs' ThinkStationP920")
  include(${CMAKE_CURRENT_LIST_DIR}/luqs-P920.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "huangjy-P920")
  message(STATUS "Huangjy' ThinkStationP920")
  include(${CMAKE_CURRENT_LIST_DIR}/huangjy-P920.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "luqs-ty")
  message(STATUS "Luqs' Taiyi")
  include(${CMAKE_CURRENT_LIST_DIR}/luqs-ty.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "lqs-tianhe")
  message(STATUS "Luqs' TianHe-2")
  include(${CMAKE_CURRENT_LIST_DIR}/luqs-th2.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "TianHe2A-pp545")
  message(STATUS "TianHe2")
  include(${CMAKE_CURRENT_LIST_DIR}/Tianhe2A-pp545.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "BSCC-T6")
  message(STATUS "BSCC-T6")
  include(${CMAKE_CURRENT_LIST_DIR}/BSCC-T6.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "NC-E")
  message(STATUS "NC-E ningxia")
  include(${CMAKE_CURRENT_LIST_DIR}/NC-E.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "HXMLinux")
  message(STATUS "Huangxm's Linux")
  include(${CMAKE_CURRENT_LIST_DIR}/huangxm-Linux.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "HXMTaiyi")
  message(STATUS "Huangxm's Taiyi")
  include(${CMAKE_CURRENT_LIST_DIR}/huangxm-Taiyi.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "HXMmac")
  message(STATUS "Huangxm's Mac")  
  include(${CMAKE_CURRENT_LIST_DIR}/huangxm-Mac.cmake)
else($ENV{MACHINE_NAME} MATCHES "poincare")
  message(STATUS "The system cannot be identified.")
endif( $ENV{MACHINE_NAME} MATCHES "poincare")

# End of the file
