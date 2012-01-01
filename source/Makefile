#  Copyright (C) 2003-2010 Fons Adriaensen <fons@kokkinizita.net>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


PREFIX = /usr
LIBDIR = lib$(SUFFIX)
VERSION = 0.8.4
DISTDIR = aeolus-$(VERSION)
CPPFLAGS += -O3 -Wall -MMD -MP -DVERSION=\"$(VERSION)\" -DLIBDIR=\"$(PREFIX)/$(LIBDIR)/aeolus\"

all:	aeolus aeolus_x11.so aeolus_txt.so


AEOLUS_O =	main.o audio.o model.o slave.o imidi.o addsynth.o scales.o \
		reverb.o asection.o division.o rankwave.o rngen.o exp2ap.o lfqueue.o
aeolus:	LDLIBS += -lclalsadrv -lclthreads -lpthread -ljack -lasound -ldl -lrt
aeolus: LDFLAGS += -L$(PREFIX)/$(LIBDIR)
aeolus:	$(AEOLUS_O)
	g++ $(LDFLAGS) -o $@ $(AEOLUS_O) $(LDLIBS)

addsynth.o:	CPPFLAGS += -fpic -D_REENTRANT
$(AEOLUS_O):
-include $(AEOLUS_O:%.o=%.d)


XIFACE_O =	styles.o mainwin.o midiwin.o audiowin.o instrwin.o editwin.o midimatrix.o \
		multislider.o functionwin.o xiface.o addsynth.o
aeolus_x11.so:	CPPFLAGS += -fpic -D_REENTRANT -I/usr/X11R6/include -I/usr/include/freetype2
aeolus_x11.so:	LDLIBS += -pthread -lclthreads -lclxclient -lXft -lX11
aeolus_x11.so:	LDFLAGS += -shared -L$(PREFIX)/$(LIBDIR) -L/usr/X11R6/$(LIBDIR)
aeolus_x11.so: $(XIFACE_O) $(LIBCLX)
	g++ $(LDFLAGS) -o $@ $(XIFACE_O) $(LDLIBS)

$(XIFACE_O):
-include $(XIFACE_O:%.o=%.d)


TIFACE_O =	tiface.o
aeolus_txt.so:	CPPFLAGS += -fpic -D_REENTRANT
aeolus_txt.so:	LDLIBS += -pthread -lclthreads -lreadline
aeolus_txt.so:	LDFLAGS += -shared -L$(PREFIX)/$(LIBDIR)
aeolus_txt.so: $(TIFACE_O)
	g++ $(LDFLAGS) -o $(@) $(TIFACE_O) $(LDLIBS)

$(TIFACE_O):
-include $(TIFACE_O:%.o=%.d)


install:	aeolus aeolus_x11.so aeolus_txt.so 
	/usr/bin/install -d $(DESTDIR)$(PREFIX)/$(LIBDIR)/aeolus
	/usr/bin/install -d $(DESTDIR)$(PREFIX)/bin
	/usr/bin/install -m 755 aeolus $(DESTDIR)$(PREFIX)/bin
	/usr/bin/install -m 755 aeolus_x11.so $(DESTDIR)$(PREFIX)/$(LIBDIR)/aeolus
	/usr/bin/install -m 755 aeolus_txt.so $(DESTDIR)$(PREFIX)/$(LIBDIR)/aeolus


clean:
	/bin/rm -f *~ *.o *.d *.a *.so aeolus

