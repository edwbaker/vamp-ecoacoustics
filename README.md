# Ecoacoustic Vamp Plugins

## Overview

This repository contains a collection of Vamp plugins for ecoacoustic analysis, implementing various indices commonly used in bioacoustics. These plugins are designed to mirror the behavior of standard R (`seewave`, `tuneR`) and Python (`maad`) implementations, providing high-performance C++ alternatives.

## Implemented Indices

1.  **ACI (Acoustic Complexity Index)**: Measures the variability of intensities. Mirrors `seewave::ACI`.
2.  **ADI (Acoustic Diversity Index)**: Measures the diversity of the spectrum. Mirrors `seewave::ADI`.
3.  **AEI (Acoustic Evenness Index)**: Measures the evenness of the spectrum. Mirrors `seewave::AEI`.
4.  **BI (Bioacoustic Index)**: Calculates the area under the curve in specific frequency bands. Mirrors `seewave::bio`.
5.  **SH (Spectral Entropy)**: Shannon entropy of the spectral distribution. Mirrors `seewave::sh`.
6.  **TH (Temporal Entropy)**: Shannon entropy of the amplitude envelope. Mirrors `seewave::th`.
7.  **H (Total Entropy)**: Product of Spectral and Temporal Entropy ($H = SH \times TH$). Mirrors `seewave::H`.
8.  **NDSI (Normalized Difference Soundscape Index)**: Ratio of biophony to anthropophony. Mirrors `seewave::NDSI`.

## Installation

## Implementation Details

The plugin exactly follows the seewave::ACI() algorithm:

1. Computes a Short-Term Fourier Transform (STFT) of the audio
2. Optionally filters to a specific frequency range
3. Divides the spectrogram into temporal windows
4. For each frequency bin within each temporal window:
   - Calculates the differences between consecutive time frames
   - Normalizes by the sum of intensities in that frequency bin
   - Sums the absolute values
5. Returns the total ACI across all frequency bins and temporal windows

## Parameters

### minfreq (Minimum Frequency)
- **Description**: Lower frequency limit in kHz
- **Default**: 0 (use full frequency range)
- **Range**: 0 to Nyquist frequency
- **Corresponds to**: `flim[1]` in seewave::ACI()

### maxfreq (Maximum Frequency)
- **Description**: Upper frequency limit in kHz
- **Default**: 0 (use full frequency range)
- **Range**: 0 to Nyquist frequency
- **Corresponds to**: `flim[2]` in seewave::ACI()

### nbwindows (Number of Windows)
- **Description**: Number of temporal windows to divide the recording
- **Default**: 1
- **Range**: 1 to 100
- **Corresponds to**: `nbwindows` in seewave::ACI()

## Usage in R with ReVAMP

```r
library(ReVAMP)
library(tuneR)

# Load audio
wave <- readWave("audio.wav")

# Basic usage - compute ACI with default parameters
result <- runPlugin("vamp-example-plugins:aci", wave, 
                   blockSize = 512, stepSize = 512)
aci_value <- result$aci[1]

# With multiple temporal windows
result <- runPlugin("vamp-example-plugins:aci", wave,
                   params = list(nbwindows = 4),
                   blockSize = 512, stepSize = 512)

# With frequency limits (2-6 kHz)
result <- runPlugin("vamp-example-plugins:aci", wave,
                   params = list(minfreq = 2, maxfreq = 6),
                   blockSize = 512, stepSize = 512)

# With overlap (50%)
result <- runPlugin("vamp-example-plugins:aci", wave,
                   blockSize = 512, stepSize = 256)
```

## Comparison with seewave::ACI()

The plugin is designed to produce equivalent results to seewave::ACI():

```r
library(seewave)
library(ReVAMP)

# seewave version
aci_seewave <- ACI(wave, wl = 512, ovlp = 0, nbwindows = 1)

# Vamp plugin version
aci_vamp <- runPlugin("vamp-example-plugins:aci", wave,
                     blockSize = 512, stepSize = 512)$aci[1]
```

### Parameter Mapping

| seewave parameter | Vamp parameter | Notes |
|------------------|----------------|-------|
| `wl` | `blockSize` | Window length for STFT |
| `ovlp` | calculated from `stepSize` | ovlp=50 means stepSize = blockSize/2 |
| `flim[1]` | `minfreq` | Lower frequency in kHz |
| `flim[2]` | `maxfreq` | Upper frequency in kHz |
| `nbwindows` | `nbwindows` | Number of temporal windows |
| `wn` | N/A | Vamp SDK handles windowing automatically |

## Technical Notes

- **Input Domain**: Frequency domain (plugin receives magnitude spectra from Vamp SDK)
- **Output**: Single value representing total ACI
- **Processing**: All frames are accumulated, final ACI computed in `getRemainingFeatures()`
- **Memory**: Stores entire spectrogram in memory for temporal window analysis

## References

Pieretti N, Farina A, Morri FD (2011) A new methodology to infer the singing activity of an avian community: the Acoustic Complexity Index (ACI). *Ecological Indicators*, 11, 868-873.

Farina A, Pieretti N, Piccioli L (2011) The soundscape methodology for long-term bird monitoring: a Mediterranean Europe case-study. *Ecological Informatics*, 6, 354-363.

Sueur J, Aubin T, Simonis C (2008). Seewave, a free modular tool for sound analysis and synthesis. *Bioacoustics*, 18: 213-226.

## License

MIT License

## Author

ReVAMP package (after Pieretti et al. 2011 and seewave::ACI implementation)
