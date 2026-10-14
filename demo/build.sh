#!/bin/sh
# Build the browser demo into demo/:
#   vamp-ecoacoustics.wasm  the plugins, built as by `make PLATFORM=wasm`
#   vamp-host.js/.wasm      the Vamp host from vamp-wasm (https://github.com/edwbaker/vamp-wasm)
#
# Needs Emscripten (emcc on PATH, or EMSDK set to an emsdk install) and a vamp-wasm
# checkout, found at VAMP_WASM_DIR (default: ../vamp-wasm, next to this repository).
#
#   sh demo/build.sh
#   python -m http.server 8000 -d demo      # then open http://localhost:8000
set -e

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DEMO="$ROOT/demo"
VAMP_WASM_DIR=${VAMP_WASM_DIR:-$ROOT/../vamp-wasm}

if ! command -v emcc >/dev/null 2>&1 && [ -n "$EMSDK" ]; then
    PATH="$EMSDK/upstream/emscripten:$PATH"
    export EM_CONFIG="${EM_CONFIG:-$EMSDK/.emscripten}"
fi
command -v emcc >/dev/null 2>&1 || { echo "emcc not found: install Emscripten or set EMSDK" >&2; exit 1; }
[ -f "$VAMP_WASM_DIR/src/WasmHost.cpp" ] || { echo "vamp-wasm not found at $VAMP_WASM_DIR (set VAMP_WASM_DIR)" >&2; exit 1; }

echo "Building plugins: demo/vamp-ecoacoustics.wasm"
cd "$ROOT"
emcc -std=c++11 -O3 -fPIC -I. -Isrc -Iext -DPOCKETFFT_NO_MULTITHREADING=1 \
    -s SIDE_MODULE=1 -s EXPORT_ALL=1 --no-entry \
    -o "$DEMO/vamp-ecoacoustics.wasm" src/*.cpp

echo "Building host: demo/vamp-host.js, demo/vamp-host.wasm"
cd "$VAMP_WASM_DIR"
SDK=inst/vamp/src
# em++, not emcc: the host must link libc++/libc++abi for itself and the plugins.
# Only the host SDK: Emscripten would resolve the plugin's own copy of the plugin SDK
# (vamp-sdk/*.cpp) to the host's, and the two copies differ.
em++ -O3 -g0 -std=c++17 -Iinst/vamp \
    -s MAIN_MODULE=1 -s MODULARIZE=1 -s EXPORT_NAME=VampHostModule \
    -s ALLOW_MEMORY_GROWTH=1 -s ENVIRONMENT=web,worker,node \
    -s EXPORTED_RUNTIME_METHODS=FS -lembind \
    -o "$DEMO/vamp-host.js" \
    $SDK/vamp-hostsdk/PluginBufferingAdapter.cpp $SDK/vamp-hostsdk/PluginChannelAdapter.cpp \
    $SDK/vamp-hostsdk/PluginHostAdapter.cpp $SDK/vamp-hostsdk/PluginInputDomainAdapter.cpp \
    $SDK/vamp-hostsdk/PluginLoader.cpp $SDK/vamp-hostsdk/PluginSummarisingAdapter.cpp \
    $SDK/vamp-hostsdk/PluginWrapper.cpp $SDK/vamp-hostsdk/RealTime.cpp \
    $SDK/vamp-hostsdk/Files.cpp src/WasmHost.cpp

echo "Done. Serve demo/ over HTTP, e.g.: python -m http.server 8000 -d demo"
