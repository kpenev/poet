ALGLIB_OBJECTS	= $(addsuffix .o, $(ALGLIB_MODULES))
ALGLIB_DEPS		= $(addsuffix .d, $(ALGLIB_MODULES)) \
				  $(addsuffix .d-e, $(ALGLIB_MODULES))

.PHONY: make_obj

make_obj: $(ALGLIB_OBJECTS)

.cpp.o:
	$(COMPILER) $(CFLAGS) -c -MMD -MP $<
	@sed -i -e '1s,\($*\)\.o[ :]*,\1.o $*.d: ,' $*.d

.PHONY: clean

clean:
	@echo "Deleting ALGLib object files ..."
	@rm -f $(ALGLIB_OBJECTS) $(ALGLIB_DEPS) 
