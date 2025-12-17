# Ecoacoustic Vamp Plugins

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
                        params = list(maxFreq = 8, dbThreshold = -40))
```

## Plugin Parameters

### ACI (aci, aci-acc)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `minFreq` | Minimum frequency (kHz) | 0 |
| `maxFreq` | Maximum frequency (kHz) | Nyquist |

### ADI (adi-acc)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `maxFreq` | Maximum frequency (kHz) | 10 |
| `dbThreshold` | Threshold in dB (relative to max) | -50 |
| `freqStep` | Frequency band width (Hz) | 1000 |

### AEI (aei-acc)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `maxFreq` | Maximum frequency (kHz) | 10 |
| `dbThreshold` | Threshold in dB (relative to max) | -50 |
| `freqStep` | Frequency band width (Hz) | 1000 |

### BI (bi-acc)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `minFreq` | Minimum frequency (kHz) | 2 |
| `maxFreq` | Maximum frequency (kHz) | 8 |

### NDSI
| Parameter | Description | Default |
|-----------|-------------|---------|
| `anthroMin` | Anthrophony min frequency (kHz) | 1 |
| `anthroMax` | Anthrophony max frequency (kHz) | 2 |
| `bioMin` | Biophony min frequency (kHz) | 2 |
| `bioMax` | Biophony max frequency (kHz) | 11 |

## Validation

All plugins are validated against their R reference implementations. Example results:

| Index | Test Files | Max Difference |
|-------|------------|----------------|
| ACI | 91 | < 0.001% |
| ADI | 93 | < 0.01% |
| AEI | 91 | < 0.02% |
| BI | 91 | < 0.01% |
| NDSI | 91 | < 0.01% |

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

## License

GPL (>= 2)

## Author

[Ed Baker](https://ebaker.me.uk)

## Contributing

Issues, feature requests and pull requests welcome at [GitHub](https://github.com/edwbaker/vamp-ecoacoustics).
