# Makefile for Ecoacoustic Vamp Plugins
# Supports: Linux, macOS, Windows (MinGW), WASM (Emscripten)

PLUGIN_NAME = vamp-ecoacoustics

# Default to host OS detection
ifdef OS
	UNAME_S := $(OS)
else
	UNAME_S := $(shell uname -s)
endif

# Platform overrides
# Usage: make PLATFORM=linux
#        make PLATFORM=macos
#        make PLATFORM=windows
#        make PLATFORM=wasm

ifeq ($(PLATFORM),)
	ifeq ($(OS),Windows_NT)
		PLATFORM = windows
	else ifneq ($(findstring MINGW,$(UNAME_S)),)
		PLATFORM = windows
	else ifneq ($(findstring CYGWIN,$(UNAME_S)),)
		PLATFORM = windows
	else ifneq ($(findstring MSYS,$(UNAME_S)),)
		PLATFORM = windows
	else ifeq ($(UNAME_S),Darwin)
		PLATFORM = macos
	else
		PLATFORM = linux
	endif
endif



CXX ?= g++

# If set to 1, delegate Windows builds to CMake/Visual Studio to produce
# a binary that matches MSVC's runtime/exports. Default to MSVC for
# Windows builds to produce portable Windows artifacts.
MSVC_BUILD ?= 1

# Platform-sensitive flags
PIC_FLAG :=
ifeq ($(PLATFORM),linux)
PIC_FLAG = -fPIC
endif
ifeq ($(PLATFORM),macos)
PIC_FLAG = -fPIC
endif

# Threads control: set THREADS=0 to compile with -DDISABLE_THREADS
THREAD_FLAG :=
ifeq ($(THREADS),0)
THREAD_FLAG = -DDISABLE_THREADS
endif

CXXFLAGS = -std=c++11 -O3 -Wall $(PIC_FLAG) $(THREAD_FLAG) -I. -Isrc -Iext

ifeq ($(PLATFORM),windows)
	CXXFLAGS += -DPOCKETFFT_NO_MULTITHREADING=1
endif

SOURCES = src/ACIPlugin.cpp src/ACIaccPlugin.cpp src/ADIaccPlugin.cpp src/AEIaccPlugin.cpp src/BIaccPlugin.cpp src/THPlugin.cpp src/SHPlugin.cpp src/HPlugin.cpp src/NDSIPlugin.cpp src/FFT.cpp src/plugins.cpp src/PluginAdapter.cpp src/RealTime.cpp
ifeq ($(PLATFORM),windows)
ifeq ($(MSVC_BUILD),1)
OBJECTS = $(SOURCES:.cpp=.obj)
else
OBJECTS = $(SOURCES:.cpp=.o)
endif
else
OBJECTS = $(SOURCES:.cpp=.o)
endif

ifeq ($(PLATFORM),macos)
	TARGET = $(PLUGIN_NAME).dylib
	LDFLAGS = -dynamiclib -install_name $(TARGET)
else ifeq ($(PLATFORM),windows)
	TARGET = $(PLUGIN_NAME).dll
	# Windows builds default to MSVC via CMake. The MSVC build path is
	# invoked when `MSVC_BUILD=1`. Legacy MinGW-specific static runtime
	# flags have been removed to avoid producing MinGW-dependent artifacts.
	LDFLAGS = -shared -Wl,--retain-symbols-file=vamp-plugin.list
else ifeq ($(PLATFORM),wasm)
	CXX = emcc
	TARGET = $(PLUGIN_NAME).wasm
	# Emscripten specific flags for wasm plugin as side module
	CXXFLAGS = -std=c++11 -O3 -fPIC $(THREAD_FLAG) -I. -Isrc -Iext -DPOCKETFFT_NO_MULTITHREADING=1
	LDFLAGS = -s SIDE_MODULE=1 -O3 -s EXPORT_ALL=1 -s LINKABLE=1 --no-entry
else
	# Linux
	TARGET = $(PLUGIN_NAME).so

	# Use version script only if present (some linkers/toolchains may not expect it)
	VERSION_SCRIPT := vamp-plugin.map
	VERSION_SCRIPT_EXISTS :=
	ifneq ($(wildcard $(VERSION_SCRIPT)),)
	VERSION_SCRIPT_EXISTS = 1
	endif

	ifeq ($(VERSION_SCRIPT_EXISTS),1)
	LDFLAGS = -shared -Wl,-soname,$(TARGET) -Wl,--version-script=$(VERSION_SCRIPT)
	else
	LDFLAGS = -shared -Wl,-soname,$(TARGET)
	endif
endif

all: $(TARGET)

ifeq ($(PLATFORM),windows)
ifeq ($(MSVC_BUILD),1)
all:
	@echo "Building with MSVC via CMake (MSVC_BUILD=1)..."
	cmake -S . -B build -G "Visual Studio 17 2022" || cmake -S . -B build
	cmake --build build --config Release --target $(TARGET) -- /m:1
	@echo "Copying MSVC-built DLL into repo root..."
	if exist build\Release\$(TARGET) copy /Y build\Release\$(TARGET) $(TARGET) >nul
	@echo "MSVC build complete."

endif
endif

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
ifeq ($(PLATFORM),windows)
	-del /Q $(subst /,\,$(OBJECTS)) $(subst /,\,$(TARGET))
else
	rm -f $(OBJECTS) $(TARGET) $(PLUGIN_NAME).wasm $(PLUGIN_NAME).js
endif
ifeq ($(PLATFORM),wasm)
	rm -f $(PLUGIN_NAME).wasm $(PLUGIN_NAME).js $(PLUGIN_NAME).wast
endif

install: $(TARGET)
	@echo "Installing for $(PLATFORM)..."
ifeq ($(PLATFORM),windows)
	@echo "Copy $(TARGET) to your Vamp Plugins directory (e.g., C:/Program Files/Vamp Plugins/)"
else ifeq ($(PLATFORM),macos)
	@echo "Copy $(TARGET) to ~/Library/Audio/Plug-Ins/Vamp/"
else ifeq ($(PLATFORM),linux)
	@echo "Copy $(TARGET) to ~/vamp/"
else ifeq ($(PLATFORM),wasm)
	@echo "Use $(TARGET) in a WebAssembly-enabled Vamp host or web application"
endif

.PHONY: all clean install
