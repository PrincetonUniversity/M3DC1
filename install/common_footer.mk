%.o : %.c
	$(CC) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.cpp
	$(CXX) $< -c $(CFLAGS) $(INCLUDE) -o $@

%.o : %.f90
	$(F90) $< -c $(F90FLAGS) -o $@

%.o : %.F90
	$(F90) $< -c $(F90FLAGS) -o $@
