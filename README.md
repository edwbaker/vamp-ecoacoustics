# Ecoacoustic Vamp Plugins

[![Build](https://github.com/edwbaker/vamp-ecoacoustics/actions/workflows/build.yml/badge.svg)](https://github.com/edwbaker/vamp-ecoacoustics/actions/workflows/build.yml)
[![Release](https://img.shields.io/github/v/release/edwbaker/vamp-ecoacoustics)](https://github.com/edwbaker/vamp-ecoacoustics/releases)
[![License: GPL v2+](https://img.shields.io/badge/license-GPL%20v2%2B-blue)](#license)

High-performance [Vamp plugins](https://vamp-plugins.org/) for ecoacoustic analysis, implementing standard acoustic indices used in bioacoustics and soundscape ecology.

## Features

- **Fast**: C++ implementations significantly faster than R  equivalents
- **Accurate**: Validated against `soundecology` and `seewave` R packages
- **Cross-platform**: Builds for Windows, macOS, Linux, and WebAssembly
- **Vamp ecosystem**: Works with Audacity, and any Vamp host

## Implemented Indices

| Plugin ID | Index | Reference Implementation |
|-----------|-------|-------------------------|
| `aci` | Acoustic Complexity Index | `seewave::ACI` |
| `aci-acc` | ACI (accumulated) | `soundecology::acoustic_complexity` |
| `adi-acc` | Acoustic Diversity Index | `soundecology::acoustic_diversity` |
| `aei-acc` | Acoustic Evenness Index | `soundecology::acoustic_evenness` |
| `bi-acc` | Bioacoustic Index | `soundecology::bioacoustic_index` |
| `ndsi` | NDSI | `soundecology::ndsi` |
| `sh` | Spectral Entropy | `seewave::sh` |
| `th` | Temporal Entropy | `seewave::th` |
| `h` | Total Entropy | `seewave::H` |

## Installation

### Pre-built Binaries

Download from [GitHub Releases](https://github.com/edwbaker/vamp-ecoacoustics/releases):

- **Windows**: `vamp-ecoacoustics.dll`
- **macOS**: `vamp-ecoacoustics.dylib`
- **Linux**: `vamp-ecoacoustics.so`

Copy to your Vamp plugins folder:
- **Windows**: `C:\Program Files\Vamp Plugins\` or `%USERPROFILE%\Vamp Plugins\`
- **macOS**: `/Library/Audio/Plug-Ins/Vamp/` or `~/Library/Audio/Plug-Ins/Vamp/`
- **Linux**: `/usr/lib/vamp/` or `/usr/local/lib/vamp/` or `~/.vamp/`

### Building from Source

#### Windows (MSVC)
```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

#### macOS / Linux
```bash
make PLATFORM=macos   # or PLATFORM=linux
```

#### WebAssembly
```bash
make PLATFORM=wasm
```

This builds an Emscripten side module, which a WebAssembly Vamp host loads with `dlopen()`. See [In a web browser](#in-a-web-browser) for a demo page that does this.

## Usage

### With R (ReVAMP package)

Makes use of the [ReVAMP](https://revamp.ebaker.me.uk) package to run Vamp plugins from R.

```r
library(ReVAMP)
library(tuneR)

# Load audio
wav <- readWave("recording.wav")

# Compute indices
aci <- runPlugin(wav, "vamp-ecoacoustics:aci-acc")
adi <- runPlugin(wav, "vamp-ecoacoustics:adi-acc")
aei <- runPlugin(wav, "vamp-ecoacoustics:aei-acc")
bi  <- runPlugin(wav, "vamp-ecoacoustics:bi-acc")
ndsi <- runPlugin(wav, "vamp-ecoacoustics:ndsi")

# With custom parameters
aei_custom <- runPlugin(wav, "vamp-ecoacoustics:aei-acc",
                        params = list(maxFreq = 8000, dbThreshold = -40))
```

### In a web browser

[`demo/`](demo/) is a web page that runs all nine indices in the browser. It loads the WebAssembly build of the plugins into the [vamp-wasm](https://github.com/edwbaker/vamp-wasm) host, which runs in a Web Worker. Recordings are analysed locally, WAV files keep their own sample rate, and blocks are passed to the plugins the way ReVAMP's `runPlugin()` passes them. Two synthetic soundscapes are built in for trying it out.

Building it needs [Emscripten](https://emscripten.org) and a vamp-wasm checkout next to this repository (or set `VAMP_WASM_DIR`):

```bash
sh demo/build.sh        # or: make demo
python -m http.server 8000 -d demo
```

Then open http://localhost:8000. The page has to be served over HTTP: browsers won't run it from `file://`.

## Plugin Parameters

### ACI (aci, aci-acc)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `minFreq` | Minimum frequency (Hz) | 0 |
| `maxFreq` | Maximum frequency (Hz) | Nyquist |
| `nbWindows` | Number of time windows to divide the file into (aci-acc). Set to `0` to auto-select windows so each is ≈ 5 s | 12 |

### ADI (adi-acc)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `maxFreq` | Maximum frequency (Hz) | 10000 |
| `dbThreshold` | Threshold in dB (relative to max) | -50 |
| `freqStep` | Frequency band width (Hz) | 1000 |

### AEI (aei-acc)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `maxFreq` | Maximum frequency (Hz) | 10000 |
| `dbThreshold` | Threshold in dB (relative to max) | -50 |
| `freqStep` | Frequency band width (Hz) | 1000 |

### BI (bi-acc)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `minFreq` | Minimum frequency (Hz) | 2000 |
| `maxFreq` | Maximum frequency (Hz) | 8000 |

### NDSI
| Parameter | Description | Default |
|-----------|-------------|---------|
| `anthroMin` | Anthrophony min frequency (Hz) | 1000 |
| `anthroMax` | Anthrophony max frequency (Hz) | 2000 |
| `bioMin` | Biophony min frequency (Hz) | 2000 |
| `bioMax` | Biophony max frequency (Hz) | 11000 |

## Validation

All plugins are validated against their R reference implementations. Example results:

| Index | Test Files | Max Difference |
|-------|------------|----------------|
| ACI | 91 | < 0.001% |
| ADI | 93 | < 0.01% |
| AEI | 91 | < 0.02% |
| BI | 91 | < 0.01% |
| NDSI | 91 | < 0.01% |

## Testing

```bash
make test                                  # Linux / macOS / MinGW
cmake -S . -B build -A x64                 # Windows (MSVC)
cmake --build build --config Release --target test_indices
ctest --test-dir build -C Release --output-on-failure
```

`tests/test_indices.cpp` runs all nine plugins the way a Vamp host does (preferred block size, zero-padded final block, `getRemainingFeatures()`) on deterministic synthetic signals, so there are no audio fixtures. It checks:

- **API**: identifiers, channel limits, output counts, parameter defaults.
- **Known answers**: a tone in the biophony band gives NDSI ≈ +1, flat noise gives ADI ≈ ln(10) and AEI ≈ 0, a steady tone has TH ≈ 1, ACI is unchanged by overall gain, stereo equals two mono runs, `reset()` is repeatable.
- **Baselines** (`tests/golden.inc`): recorded values that catch unintended numerical changes, including differences between platforms. These come from this code, not from R. After a deliberate change, regenerate them with `tests/test_indices --write-golden tests/golden.inc`.

## Performance

Typical speedups vs R implementations (varies by file):

| Index | vs soundecology | vs seewave |
|-------|-----------------|------------|
| ACI | ~8x | ~89x |
| ADI | ~3-6x | - |
| AEI | ~2-4x | - |
| BI | ~7x | - |
| NDSI | ~30x | - |
| H | - | ~4x |
| SH | - | ~2x |
| TH | - | ~2x |

## References

- Pieretti N, Farina A, Morri FD (2011). A new methodology to infer the singing activity of an avian community: the Acoustic Complexity Index (ACI). *Ecological Indicators*, 11, 868-873.
- Villanueva-Rivera LJ, Pijanowski BC, Doucette J, Pekin B (2011). A primer of acoustic analysis for landscape ecologists. *Landscape Ecology*, 26, 1233-1246.
- Sueur J, Aubin T, Simonis C (2008). Seewave, a free modular tool for sound analysis and synthesis. *Bioacoustics*, 18, 213-226.

## Citing

If you use these plugins, please cite them using the metadata in [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button). Each GitHub release is archived on Zenodo, which gives it a citable DOI.

## License

GPL (>= 2)

## Author

[Ed Baker](https://ebaker.me.uk)

## Contributing

Issues, feature requests and pull requests welcome at [GitHub](https://github.com/edwbaker/vamp-ecoacoustics).
