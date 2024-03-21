TARGETDIR = ./bin
TARGETNAME = a.out

SUBDIRS = 

# Library paths
LIB_PATHS = 

# Library names
LIB_NAMES += pthread

# C++
CXX = g++-13
CXXFLAGS += -std=c++17
INCLUDES = ./include

# Linker
LINK = $(CXX)
LDFLAGS += -std=c++17

#
# Compile options
#
# C
ifeq "$(CXX)" "icpx"
  CXXFLAGS += 
  ifdef OPT
	CXXFLAGS += -Ofast -fopenmp
	LDFLAGS += -fopenmp
  else
	CXXFLAGS += -O0 -g3 -Wuninitialized -DMYDEBUG 
  endif
else
  CXXFLAGS += 
  ifdef OPT
	CXXFLAGS += -Ofast -fopenmp
	LDFLAGS += -fopenmp
  else
	CXXFLAGS += -O0 -g3 -fbounds-check -Wuninitialized -DMYDEBUG # -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += # -fsanitize=address -fno-omit-frame-pointer
  endif
endif

# Directory
SRCDIR = ./src
OBJDIR = ./obj

#
# 
#
LDFLAGS += $(addprefix -L,$(LIB_PATHS)) $(addprefix -l,$(LIB_NAMES))

TARGET = $(TARGETDIR)/$(TARGETNAME)


ALL_SRCS := $(wildcard $(SRCDIR)/*.cpp)

TEST_SRCS = $(wildcard $(SRCDIR)/test_*.cpp)

NON_TEST_SRCS := $(filter-out $(TEST_SRCS), $(ALL_SRCS))

OBJS = $(addprefix $(OBJDIR)/, $(patsubst %.cpp,%.o,$(notdir $(NON_TEST_SRCS))))
DEPS += $(addprefix $(OBJDIR)/, $(patsubst %.cpp,%.d,$(notdir $(NON_TEST_SRCS))))
DEPS += $(addprefix $(OBJDIR)/, $(patsubst %.cpp,%.d,$(notdir $(TEST_SRCS))))


TEST_BINS = $(patsubst $(SRCDIR)/test_%.cpp,$(TARGETDIR)/test_%,$(TEST_SRCS))

.PHONY: $(SUBDIRS) all clean c
all: $(SUBDIRS) $(TARGET) #doxygen

$(SUBDIRS):
	@$(MAKE) -C $@

$(TARGET): $(OBJS)
	@-mkdir -p $(TARGETDIR)
	$(LINK) $(OBJS) $(LDFLAGS) -o $@


# C++
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(addprefix -I,$(INCLUDES))  -o $@ -c $<

# 
clean: c
	@for SUBDIR in $(SUBDIRS); do $(MAKE) clean -C $$SUBDIR ; done

# 
c:
	@-rm -rf $(TARGET) $(OBJDIR) $(TEST_BINS) $(gtestdir)
  

-include $(DEPS)

$(OBJDIR)/%.d: $(SRCDIR)/%.cpp
	@-mkdir -p $(OBJDIR)
	@g++ -std=c++17 $(addprefix -I,$(INCLUDES)) -MM -MP -MT '$(OBJDIR)/$*.o' $< > $(OBJDIR)/$*.d