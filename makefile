#makefile for the PAM
MAIN = main
RM = rm -f
DIR = SRC3
OBJ =  $(DIR)/funct.o $(DIR)/$(MAIN).o $(DIR)/dinftr.o $(DIR)/makegrid.o $(DIR)/spheat.o $(DIR)/uhf.o



F90=ifort -O3 -c
link=ifort -O3   

EXENAME = ./xpam



$(EXENAME): $(OBJ)
	$(link) $(OBJ) -o  $(EXENAME)

$(DIR)/funct.o: $(DIR)/funct.f
	cd $(DIR);$(F90) funct.f;cd ..
$(DIR)/$(MAIN).o: $(DIR)/$(MAIN).f
	cd $(DIR);$(F90) $(MAIN).f;cd ..
$(DIR)/makegrid.o: $(DIR)/makegrid.f
	cd $(DIR);$(F90) makegrid.f;cd ..
$(DIR)/dinftr.o: $(DIR)/dinftr.f
	cd $(DIR);$(F90) dinftr.f;cd ..
$(DIR)/spheat.o: $(DIR)/spheat.f
	cd $(DIR);$(F90) spheat.f;cd ..
$(DIR)/uhf.o: $(DIR)/uhf.f
	cd $(DIR);$(F90) uhf.f;cd ..

clean:
	$(RM) $(OBJ) 
realclean:
	$(RM) $(OBJ) ./xpam

