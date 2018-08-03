##############################################################################
# 	make DemoFoam
# 	make DemoMCell
#	make -f Makefile DemoFoam
#	make -f Makefile DemoMCell
#	make -f Makefile map-eps
##############################################################################
all:	DemoFoam
##############################################################################
.SUFFIXES:  .eps .dvi .tex .txi .html .bib .blg .bbl .aux .txp
#---------------------------------------------------------------------------
XPS    = gv -geometry +0+10 -bg white -fg black  -magstep 0 
XPS4   = gv -bg white -fg black  -magstep -4 
#---------------------------------------------------------------------------
########### Linux flags:
FFLAGS  =  -fno-automatic -ffixed-line-length-132 -fbackslash -C
########### HP  flags:
### +es for 100 or more columns, +B for backslash in strings, +K  for global SAVE
#FFLAGS  =  -K +es +B -O
########### IBM/AIX  flags:
# fixed=132 for 132 columns,
# qextname for C-style external names with underscore at the end of external name
#FFLAGS  =  -qfixed=132 -qextname -O -C -qextchk
#FFLAGS  =  -qfixed=132 -qextname -g -C -qextchk -qinitauto=FF -qflttrap=overflow:underflow:zerodivide:invalid:enable
########### DEC/ALPHA  flags:
#FFLAGS  =  -extend_source -qfixed=120 -qextname -O  -C -qextchk  
#---------------------------------------------------------------------------
#FFLAGS  =  -qfixed=132 -qextname -g  -C -qextchk -qflttrap=overflow:zerodivide:invalid:enable  -qsigtrap
#LDFLAGS =  -bloadmap:LinkMap
F77 = gfortran
LD = $(F77)
LDFLAGS  =
##############################################################################
.f.o:
	$(F77) $(FFLAGS) -c $<
############# postscript translation (one page)
.txp.dvi:
	(latex $< )
############# option E is for tight bounding on one
.dvi.eps:
	(dvips -E $< -o $@;	$(XPS4) $@ )
###	(dvips -E $< -o )
##############################################################################
SRCLIST= " "
##############################################################################
FoamA.o:	FoamA.h
MCellA.o:	MCellA.h
PseuMar.o:	PseuMar.h
GLK.o:		GLK.h
##############################################################################
COMMAND3 = DemoFoam.exe
OBJECTS3 = DemoFoam.o FoamA.o PseuMar.o GLK.o ranlux.o WtLimit.o
$(COMMAND3): $(OBJECTS3)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS3)
#########
DemoFoam:	$(COMMAND3)
	(time ./$(COMMAND3); diff output-DemoFoam-linux output-DemoFoam | more )
#############################
DemoFoamMap:	DemoFoam
	(latex demo-cell-map.txp; dvips  -E demo-cell-map.dvi -o demo-cell-map.eps;  $(XPS) demo-cell-map.eps&)
##############################################################################
##############################################################################
COMMAND5 = DemoMCell.exe
OBJECTS5 = DemoMCell.o MCellA.o PseuMar.o GLK.o ranlux.o WtLimit.o
$(COMMAND5): $(OBJECTS5)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS5)
#########
DemoMCell:	$(COMMAND5)
	(time ./$(COMMAND5); diff output-DemoMCell-linux output-DemoMCell | more )
#############################
DemoMCellMap:	DemoMCell
	(latex demo-cell-map.txp; dvips  -E demo-cell-map.dvi -o demo-cell-map.eps;  $(XPS) demo-cell-map.eps&)
##############################################################################
line_count:
	(cat *.f *.h | wc -l )
##############################################################################
