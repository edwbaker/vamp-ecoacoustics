# Makefile for Ecoacoustic Vamp Plugins
# Supports: Linux, macOS, Windows (MinGW), WASM (Emscripten)

PLUGIN_NAME = vamp-aci

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

CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -fPIC -I. -Isrc -Iext

SOURCES = src/ACIPlugin.cpp src/AmplitudeIndexPlugin.cpp src/TemporalEntropyPlugin.cpp src/ACImtPlugin.cpp src/FFT.cpp src/plugins.cpp src/PluginAdapter.cpp src/RealTime.cpp
OBJECTS = $(SOURCES:.cpp=.o)

ifeq ($(PLATFORM),macos)
	TARGET = $(PLUGIN_NAME).dylib
	LDFLAGS = -dynamiclib -install_name $(TARGET)
else ifeq ($(PLATFORM),windows)
	TARGET = $(PLUGIN_NAME).dll
	LDFLAGS = -shared -Wl,--retain-symbols-file=vamp-plugin.list -static-libgcc -static-libstdc++ -Wl,-Bstatic -lwinpthread -Wl,-Bdynamic
else ifeq ($(PLATFORM),wasm)
	CXX = emcc
	TARGET = $(PLUGIN_NAME).wasm
	# Emscripten specific flags
	CXXFLAGS = -std=c++11 -O3 -fPIC -I. -Isrc -Iext -s ALLOW_MEMORY_GROWTH=1
	LDFLAGS = -s SIDE_MODULE=1 -O3 -s EXPORTED_FUNCTIONS="['_vampGetPluginDescriptor']"
else
	# Linux
	TARGET = $(PLUGIN_NAME).so
	LDFLAGS = -shared -Wl,-soname,$(TARGET) -Wl,--version-script=vamp-plugin.map
endif

all: $(TARGET)

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

install: $(TARGET)
	@echo "Installing for $(PLATFORM)..."
ifeq ($(PLATFORM),windows)
	@echo "Copy $(TARGET) to your Vamp Plugins directory (e.g., C:/Program Files/Vamp Plugins/)"
else ifeq ($(PLATFORM),macos)
	@echo "Copy $(TARGET) to ~/Library/Audio/Plug-Ins/Vamp/"
else ifeq ($(PLATFORM),linux)
	@echo "Copy $(TARGET) to ~/vamp/"
endif

.PHONY: all clean install
