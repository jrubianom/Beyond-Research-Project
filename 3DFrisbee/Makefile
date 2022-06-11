include local_config.mk

#Compiling parameters
#path en spack
#CXX = mpic++
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS) -I$(GENERAL)
RUN = mpirun -np $(PROCCESORS) ./
#RUN = ./
SOURCES = $(wildcard code/$(file).cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)
SOURCESMAIN = $(wildcard temporal/*.cpp)
DEPENDENCIESMAIN = $(SOURCESMAIN:temporal/%.cpp=.objects_temporal/%.o)
PATH_GH_REPOS = /home/live/repos/Beyond-Research-Project

.PHONY: all main mesh graph clean oclean

all: $(file).x

test:
	@echo $(MFEM_INSTALL_DIR)
	@echo $(GENERAL)

$(file).x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$@
	cp $@ executable

main: $(DEPENDENCIESMAIN)
	@echo -e 'Compiling' $@.x '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@.x
	@echo -e 'Done!\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$@.x
	cp $@.x executable
	cp temporal/* .Backup

mesh: current_mesh.geo
	gmsh current_mesh.geo -format msh2 -o mesh.msh -3 > /dev/null

.objects/%.o: code/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

.objects_temporal/%.o: temporal/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

gp:
	cp -rf $(dir) $(PATH_GH_REPOS)
	git -C $(PATH_GH_REPOS) add $(dir)
	git -C  $(PATH_GH_REPOS) commit -m"updated"
	git -C  $(PATH_GH_REPOS) push

#Activae glvis server
actgl: main
	xterm -hold -e ~/repos/GLVIS/glvis-4.1/build/./glvis &
	./main.x

opengl:
	xterm -hold -e ~/repos/GLVIS/glvis-4.1/build/./glvis &

pushex:
	cp code/examples/$(ex) ./code
	emacs code/$(ex) &



prvwpull:
	cp /media/sf_Paraview_samples/temporal_prvw ./temporal_prvw
prvwpush:
	rm -rf /media/sf_Paraview_samples/temporal_prvw
	mkdir /media/sf_Paraview_samples/temporal_prvw
	mv -r temporal_prvw/* /media/sf_Paraview_samples/temporal_prvw

mvall:
	cp -r /media/sf_Paraview_samples/*.msh ../data/

execlean:
	rm -rf executable/*.x

xclean:
	rm -rf *.x

oclean:
	@rm -rf .objects/*.o
