# ===========================================================================
#  Makefile bms_ap
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile		[uqmd.f]
##			TYPE	=	operating system	['uname']
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  


# ----------------------------------------- 

ifeq "$(TYPE)" ""
   TYPE := $(shell uname)
endif


ifeq "$(TYPE)" "AIX" 

FC      	= 	xlf
LD              =       xlf
RM              =       rm
O		=	.o
#FFLAGS          =       -O  -qextname  -qcheck -qextchk -bloadmap:loadmap.out
#LDFLAGS         =       -O  -qextname  -qcheck -qextchk -bloadmap:loadmap.out
FFLAGS     	=       -O2 -qextname 
LDFLAGS    	=       -O2 -qextname
SYSTEMFILES	=	$(SRCAIX)

endif

ifeq "$(TYPE)" "GNU"

FC              =       gfortran
LD		=	gfortran
RM              =       rm
O               =       .o
FFLAGS          =       -O3 -march=nocona -m64  
LDFLAGS         =       -O3 -march=nocona -m64
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "Linux" 

FC              =       gfortran
LD              =       gfortran
RM              =       rm 
O               =       .o
FFLAGS          =       -O
LDFLAGS         =       -O 
SYSTEMFILES     =       $(SRCGNU)

endif


ifeq "$(TYPE)" "HP-UX"

FC              =       f77
LD		=	f77
RM		=	rm
O               =       .o
FFLAGS          =       -w
LDFLAGS         =       -g
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "IRIX64"

FC              =       f77
LD		=	f77
RM		=	rm
O               =       .o
FFLAGS          =       -g -n32 -trapuv -C 
LDFLAGS         =       -g -n32
SYSTEMFILES     =       $(SRCSGI)

endif

ifeq "$(TYPE)" "PURE"

FC              =       f77
LD              =       purify f77
RM              =       rm 
O               =       .o 
FFLAGS          =       -w 
LDFLAGS         =       -g 
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "ALPHA"

FC              =       f77
LD              =       f77
RM              =       rm 
O               =       .o 
FFLAGS          =       -V -g  -C -align dcommons   -check overflow
LDFLAGS         =       -V -g  -C -align dcommons   -check overflow
SYSTEMFILES     =       $(SRCALPHA)

endif

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	bms_ap
endif

SRC		=	iniflw.f flow.f\
			oscarnxtev.f bms_ap.f obstring.f cmshif.f angle.f\
			input.f gridin.f output.f obsvalue.f evalue.f normvalue.f \
			cubic.f ratio.f jacobi.f eigsrt.f sqr.f reftrans.f readwt.f

INC		= 	ucoms.f hicoms.f

# -------------------------------------------------

OBJDIR		=	obj_$(TYPE)
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN).$(TYPE)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.f:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(LD) $(LDFLAGS) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS) loadmap.out

distclean:	
		-rm $(TARGET) loadmap.out
		-rm -r obj_*

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------

./uqmdnxtev.f:	ucoms.f 
./bms_ap.f:	ucoms.f hicoms.f oscarnxtev.f readwt.f
./obsvalue.f: 	ucoms.f	hicoms.f
./evalue.f:	ucoms.f hicoms.f
./frag.f: 	ucoms.f hicoms.f
./normvalue.f:	ucoms.f hicoms.f
./output.f:	ucoms.f hicoms.f
./reftrans.f:	ucoms.f hicoms.f
./obstring.f:	obsvalue.f evalue.f
./readwt.f:     ucoms.f
