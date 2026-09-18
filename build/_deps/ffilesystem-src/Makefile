NAME := fs_cli
FNAME := filesystem_cli
LIB := libfilesystem.a

cpp := 1

ifeq ($(origin CC),default)
CC = gcc
endif

ifeq ($(origin CXX),default)
CXX = g++
endif

FC := gfortran

# Clang, including AppleClang:
# make CC=clang CXX=clang++

# NVHPC:
# make CC=nvc CXX=nvc++ FC=nvfortran

# Intel oneAPI
# make CC=icx CXX=icpx FC=ifx

BUILD_DIR := build-make

INC := -Iinclude/

cpp = 1

# optional, but useful
cfeat =
cppfeat =

oflags = -O2 -DNDEBUG

CXXFLAGS := -std=c++17 $(oflags) $(cppfeat) $(INC) -DHAVE_CXX_FILESYSTEM -Dffilesystem_extra
CFLAGS := $(oflags) $(cfeat) $(INC)
FFLAGS := $(oflags)

ARFLAGS := rcs

comdir = src/
fdir = $(comdir)fortran/

COMM_SRCS = \
    $(comdir)absolute.cpp \
	$(comdir)pure.cpp \
	$(comdir)copy.cpp \
	$(comdir)inquire.cpp \
	$(comdir)c.cpp \
	$(comdir)disk.cpp \
	$(comdir)equivalent.cpp \
	$(comdir)env.cpp \
	$(comdir)executable.cpp \
	$(comdir)home.cpp \
	$(comdir)lang.cpp \
	$(comdir)lexical.cpp \
	$(comdir)limits.cpp \
	$(comdir)log.cpp \
	$(comdir)memory.cpp \
	$(comdir)mkdir.cpp \
	$(comdir)move.cpp \
	$(comdir)normalize.cpp \
	$(comdir)os.c \
	$(comdir)parent.cpp \
	$(comdir)partition.cpp \
	$(comdir)permissions.cpp \
	$(comdir)platform.cpp \
	$(comdir)relative.cpp \
	$(comdir)resolve.cpp \
	$(comdir)size.cpp \
	$(comdir)space.cpp \
	$(comdir)symlink.cpp \
	$(comdir)tempdir.cpp \
	$(comdir)time.cpp \
	$(comdir)touch.cpp \
	$(comdir)ulimit.cpp \
	$(comdir)which.cpp \
	$(comdir)windows.cpp \
	$(comdir)extra/case.cpp \
	$(comdir)extra/compiler.cpp \
	$(comdir)extra/component.cpp \
	$(comdir)extra/cygwin.cpp \
	$(comdir)extra/exepath.cpp \
	$(comdir)extra/libpath.cpp \
	$(comdir)extra/locale.cpp \
	$(comdir)extra/owner.cpp \
	$(comdir)extra/random.cpp \
	$(comdir)extra/removable.cpp \
	$(comdir)extra/shell.cpp \
	$(comdir)extra/sysctl.cpp \
	$(comdir)extra/uid.cpp \
	$(comdir)extra/uname.cpp \
	$(comdir)extra/winsock.cpp

OBJS := $(COMM_SRCS:%=$(BUILD_DIR)/%.o)

fbd = $(BUILD_DIR)/$(fdir)

lib = $(BUILD_DIR)/$(LIB)
main = $(BUILD_DIR)/$(NAME)

ifeq ($(OS),Windows_NT)
	SHELL := pwsh.exe
	.SHELLFLAGS := -Command
	LDFLAGS := -lws2_32 -lOle32 -lShell32 -luuid -lUserenv -lSecur32 -lShlwapi
	RM := Remove-Item -Recurse
	MKDIR := New-Item -ItemType Directory -Force -Path
	MKDIR_QUIET := | Out-Null
	COMMENT = ".SHELLFLAGS -Command needed to get Make to use powershell rather than cmd"
else
	RM := rm -rf
	MKDIR := mkdir -p
	MKDIR_QUIET :=
endif

ifeq (icpx,$(findstring icpx,$(CXX)))
  CXXLIBS := -lstdc++
else ifeq (nvc++, $(findstring nvc++,$(CXX)))
  CXXLIBS := -lstdc++
else ifeq (clang++,$(findstring clang++,$(CXX)))
  CXXLIBS := -lc++
else ifeq  (g++,$(findstring g++,$(CXX)))
  CXXLIBS := -lstdc++
endif

ifeq (gfortran,$(findstring gfortran,$(FC)))
  FFLAGS += -J$(BUILD_DIR)
endif

.PHONY: $(main)
.DEFAULT_GOAL := all

$(OBJS): Makefile

all: mbd $(lib) $(main) $(main_f)

mbd: $(fbd)

$(fbd):
	@$(MKDIR) $(fbd) $(MKDIR_QUIET)

$(lib): $(OBJS) $(FOBJS)
	$(AR) $(ARFLAGS) $@ $?

# cosmocc wants OBJS instead of $(lib). The latter is OK for standard compilers.
$(main): app/main.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $< $(LDFLAGS)

$(BUILD_DIR)/%.c.o: %.c
	@$(MKDIR) $(dir $@) $(MKDIR_QUIET)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.cpp.o: %.cpp
	@$(MKDIR) $(dir $@) $(MKDIR_QUIET)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Fortran
ifdef fortran
main_f = $(BUILD_DIR)/$(FNAME)

FSRCS = $(fdir)filesystem.F90
FOBJS := $(FSRCS:%=$(BUILD_DIR)/%.o)

$(main_f): app/fortran/main.f90 $(FOBJS) $(lib)
	$(FC) $(FFLAGS) $(FOBJS) $(lib) -o $@ $< $(LDFLAGS) $(CXXLIBS)

$(BUILD_DIR)/%.f90.o: %.f90
	@$(MKDIR) $(dir $@) $(MKDIR_QUIET)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.F90.o: %.F90
	@$(MKDIR) $(dir $@) $(MKDIR_QUIET)
	$(FC) $(FFLAGS) -c $< -o $@
endif

clean:
	$(RM) $(BUILD_DIR)
