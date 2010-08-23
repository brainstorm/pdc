F90= ifort -FR
OMPF90= $(F90) -openmp -fpp
CPPFLAG=
FOPT= -O2 -g 
FOPT= -O2
RTC=rtc.o
ESSLSMPLIB=-L/pdc/vol/mkl/9.1.023/lib/em64t -lmkl_lapack -lmkl -lguide -lpthread -Wl,--rpath,/pdc/vol/mkl/9.1.023/lib/em64t
ESSLLIB=$(ESSLSMPLIB)
# For the IBM Nighthawk use the following
#F90= xlf90
#OMPF90= xlf90_r -qsmp=omp
#CPPFLAG= -WF,
#FOPT= -O2 -g  ## -g Breaks the compiler on the twod target
#FOPT= -O2
#ESSLSMPLIB=-lesslsmp
#ESSLLIB=-lessl


default:
	@echo "Try: make Prog1Ex1  or make Prog1Ex2 or similar..."

all: Prog1 Prog2 twod-serial twod Prog4 Prog5

Prog1: Prog1Ex1  Prog1Ex2  Prog1Ex3  Prog1Ex4
	@echo "All targets Prog1 made"

Prog2: Prog2Serial Prog2Ex1 Prog2Ex2 Prog2Ex3 Prog2Ex4 Prog2Ex5
	@echo "All targets Prog2 made"

Prog4: Prog4Serial Prog4OMP
	@echo "All targets Prog4 made"

Prog5: Prog5Serial Prog5OMP
	@echo "All targets Prog5 made"

.f.o:
	$(OMPF90) $(FOPT) -c $?

twod: twod.f omp_smeds.o
	$(OMPF90) $(FOPT) -o $@ $^ $>

twod-serial: twod.f
	$(F90) $(FOPT) -o $@ $^ $>

Prog1Ex1: Prog1.F omp_smeds.o dummy.o
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=1 -o $@ $^ $>

Prog1Ex2: Prog1.F omp_smeds.o dummy.o
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=2 -o $@ $^ $>

Prog1Ex3: Prog1.F omp_smeds.o dummy.o
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=3 -o $@ $^ $>

Prog1Ex4: Prog1.F omp_smeds.o dummy.o
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=4 -o $@ $^ $>

Prog2Serial: Prog2.F omp_smeds.o dummy.o $(RTC)
	$(F90) $(FOPT) -o $@ $^ $> $(ESSLLIB) $(RTC)

Prog2Ex1: Prog2.F omp_smeds.o dummy.o $(RTC)
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=1 -o $@ $^ $> $(ESSLLIB) $(RTC)

Prog2Ex2: Prog2.F omp_smeds.o dummy.o $(RTC)
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=2 -o $@ $^ $> $(ESSLSMPLIB) $(RTC)

Prog2Ex3: Prog2.F omp_smeds.o dummy.o $(RTC)
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=3 -o $@ $^ $> $(ESSLSMPLIB) $(RTC)

Prog2Ex4: Prog2.F omp_smeds.o dummy.o $(RTC)
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=4 -o $@ $^ $> $(ESSLSMPLIB) $(RTC)

Prog2Ex5: Prog2.F omp_smeds.o dummy.o
	$(OMPF90) $(FOPT) $(CPPFLAG)-DEXAMPLE=5 -o $@ $^ $> $(ESSLSMPLIB) $(RTC)

Prog4Serial: Prog4.f omp_smeds.o $(RTC)
	$(F90) $(FOPT) -o $@ $^ $> $(RTC)

Prog4OMP: Prog4.f omp_smeds.o $(RTC)
	$(OMPF90) $(FOPT) -o $@ $^ $> $(RTC)

Prog5Serial: Prog5.f omp_smeds.o $(RTC)
	$(F90) $(FOPT) -o $@ $^ $> $(RTC)

Prog5OMP: Prog5.f omp_smeds.o $(RTC)
	$(OMPF90) $(FOPT) -o $@ $^ $> $(RTC)

clean:
	rm -f Prog1Ex? Prog2Ex? *.o *.mod twod twod-serial \
	Prog4Serial Prog4OMP Prog5Serial Prog5OMP Prog5Serial
