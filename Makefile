# make these work for both my remotes and my local

#user=$(shell whoami)
#ifeq ($(user),leander)
#  PATHTOCLASS = /usr/local/class_public-3.0.1
#  PATHTOFFTW = /usr/local/fftw-3.3.9
#else
#  PATHTOCLASS = /home/lthiele/class_public
#  PATHTOFFTW = /usr/local/fftw/gcc/3.3.4
#endif

#PATHTOGSL=/opt/apps/intel18/gsl/2.6
#PATHTOGSL=/opt/apps/intel24/gsl/2.7.1
#PATHTOGSL=/work/07833/tg871330/stampede3/software/gsl-2.6
PATHTOGSL=${TACC_GSL_DIR}
PATHTOCLASS=/work2/07833/tg871330/stampede3/software/class_public-3.2.0
PATHTOFFTW=${TACC_FFTW3_DIR}
#PATHTOFFTW=/moto/opt/fftw3/fftw-3.3.10
#PATHTOFFTW=/opt/apps/intel18/impi18_0/fftw3/3.3.8
#PATHTOFFTW=/opt/apps/intel24/impi21/fftw3/3.3.10
#PATHTOFFTW=/work/07833/tg871330/stampede3/software/fftw-3.3.8

CC = gcc
ARCHIVE = libhmpdf.a
SHARED = libhmpdf.so

CFLAGS = --std=gnu99 -fPIC -Wall -Wextra -Wpedantic -Wno-variadic-macros -Winline -DHAVE_INLINE -DDEBUG -DARICO20
#CFLAGS = --std=gnu99 -fPIC -Wall -Wextra -Wpedantic -Wno-variadic-macros -Winline -DHAVE_INLINE -DDEBUG -DARICO20 -DHMFSWAP
#CFLAGS = --std=gnu99 -fPIC -Wall -Wextra -Wpedantic -Wno-variadic-macros -Winline -DHAVE_INLINE -DDEBUG -DARICO20 -DSAVE_PROF -DSAVE_HMF -DSAVE_SIGMA_NU

OPTFLAGS = -O4 -ggdb3 -ffast-math
OMPFLAGS = -fopenmp

INCLUDE = -I./include 
INCLUDE += -I$(PATHTOCLASS)/include \
           -I$(PATHTOCLASS)/external/HyRec2020 \
           -I$(PATHTOCLASS)/external/RecfastCLASS \
           -I$(PATHTOCLASS)/external/heating \
           -I$(PATHTOFFTW)/api \
	   -I$(PATHTOFFTW)/include \
	   -I$(PATHTOGSL)/include \

LINKER = -L$(PATHTOCLASS)
LINKER += -L$(PATHTOFFTW)/lib -L$(PATHTOGSL)/lib -lclass -lgsl -lgslcblas -lm -lfftw3

SRCDIR = ./src
OBJDIR = ./obj
OUTDIR = ./lib
SODIR  = .

.PHONY: directories

all: directories $(SHARED) $(ARCHIVE)

$(SHARED): $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
	$(CC) -shared -o $(SODIR)/$@ $^ $(LINKER) $(OMPFLAGS)

$(ARCHIVE): $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
	ar -r -o $(OUTDIR)/$@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(OPTFLAGS) $(OMPFLAGS) -o $@ $<

directories: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)
	mkdir -p $(OUTDIR)

.PHONY: python
python:
	echo $(shell pwd) > hmpdf/PATHTOHMPDF.txt
	pip install . --user

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o
	rmdir $(OBJDIR)
	rm $(SODIR)/$(SHARED)
	rm $(OUTDIR)/$(ARCHIVE)
	rmdir $(OUTDIR)
