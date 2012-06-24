cmake_minimum_required(VERSION 2.8)
find_package ( PkgConfig REQUIRED )
project (aeolus)
add_executable(aeolus
addsynth.cc
asection.cc
a_thread.cc
audio.cc
audiowin.cc
button.cc
clalsadrv.cc
division.cc
editwin.cc
enumip.cc
exp2ap.cc
functionwin.cc
imidi.cc
instrwin.cc
itc_ctrl.cc
itc_ip1q.cc
itc_mesg.cc
lfqueue.cc
main.cc
mainwin.cc
mclist.cc
menuwin.cc
meter.cc
midimatrix.cc
midiwin.cc
model.cc
multislider.cc
p_thread.cc
rankwave.cc
reverb.cc
rngen.cc
scale.cc
scales.cc
scroll.cc
slave.cc
slider.cc
styles.cc
textip.cc
textmsg.cc
tiface.cc
xdisplay.cc
xdraw.cc
xhandler.cc
xiface.cc
xresman.cc
xwindow.cc
)

#XIFACE_O =	styles.o mainwin.o midiwin.o audiowin.o instrwin.o editwin.o midimatrix.o \
#		multislider.o functionwin.o xiface.o addsynth.o

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")                            

set(aeolus_VERSION_MAJOR 0)
set(aeolus_VERSION_MINOR 8.4)
set(CMAKE_BUILD_TYPE Release)
set(VERSION "${aeolus_VERSION_MAJOR}.${aeolus_VERSION_MINOR}")
set(PREFIX "${VERSION}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -MMD")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -MP")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=\"\\\"${VERSION}\\\"\"")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLIBDIR=\"${PREFIX}/${LIBDIR}\"")

set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -pthread")

pkg_check_modules ( ALSA alsa>=0.9.1 )
pkg_check_modules ( FREETYPE freetype>=1.0.0 )
find_package (ALSA)
if (ALSA_FOUND)
  include_directories(${ALSA_INCLUDE_DIRS})
  target_link_libraries (aeolus ${ALSA_LIBRARIES})
endif (ALSA_FOUND)


pkg_check_modules ( JACK jack )
set ( JACK_SUPPORT ${JACK_FOUND} )

if (JACK_FOUND)
  include_directories(${JACK_INCLUDE_DIRS})
  target_link_libraries (aeolus ${JACK_LIBRARIES})
endif (JACK_FOUND)

find_package (Freetype)
if (FREETYPE_FOUND)
  include_directories(${FREETYPE_INCLUDE_DIRS})
  target_link_libraries (aeolus ${FREETYPE_LIBRARIES})
endif (FREETYPE_FOUND)

if(UNIX AND NOT APPLE)
    target_link_libraries(aeolus rt)
    target_link_libraries(aeolus Xft)
    target_link_libraries(aeolus X11)
    target_link_libraries(aeolus readline)
endif()
