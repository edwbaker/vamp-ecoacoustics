/*
    Ecoacoustic Vamp Plugins - browser demo: spectrogram and frequency-band overlay

    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

const FFT_SIZE = 1024;
const DB_RANGE = 80;
const MAX_COLUMNS = 1600;

// matplotlib's inferno colour map, interpolated from 11 anchors into a 256-entry table
const INFERNO = ['#000004', '#160b39', '#420a68', '#6a176e', '#932667', '#bc3754',
                 '#dd513a', '#f37819', '#fca50a', '#f6d746', '#fcffa4'];
const LUT = (() => {
    const rgb = INFERNO.map(hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16)));
    const lut = new Uint8ClampedArray(256 * 3);
    for (let i = 0; i < 256; i++) {
        const x = i / 255 * (rgb.length - 1);
        const k = Math.min(rgb.length - 2, Math.floor(x));
        for (let c = 0; c < 3; c++) lut[i * 3 + c] = rgb[k][c] + (rgb[k + 1][c] - rgb[k][c]) * (x - k);
    }
    return lut;
})();

const HANN = Float64Array.from({ length: FFT_SIZE }, (_, i) => 0.5 - 0.5 * Math.cos(2 * Math.PI * i / (FFT_SIZE - 1)));

// In-place iterative radix-2 FFT
function fft(re, im) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
        let bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            [re[i], re[j]] = [re[j], re[i]];
            [im[i], im[j]] = [im[j], im[i]];
        }
    }
    for (let len = 2; len <= n; len <<= 1) {
        const half = len >> 1;
        const wr = Math.cos(-2 * Math.PI / len), wi = Math.sin(-2 * Math.PI / len);
        for (let i = 0; i < n; i += len) {
            let cr = 1, ci = 0;
            for (let k = 0; k < half; k++) {
                const a = i + k, b = a + half;
                const tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
                re[b] = re[a] - tr; im[b] = im[a] - ti;
                re[a] += tr; im[a] += ti;
                [cr, ci] = [cr * wr - ci * wi, cr * wi + ci * wr];
            }
        }
    }
}

function sizeCanvas(canvas, maxWidth = Infinity) {
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.max(1, Math.min(maxWidth, Math.round(canvas.clientWidth * dpr)));
    canvas.height = Math.max(1, Math.round(canvas.clientHeight * dpr));
    return canvas.getContext('2d');
}

// Linear-frequency spectrogram of `samples`, 0 Hz at the bottom and Nyquist at the top
export function drawSpectrogram(canvas, samples) {
    const ctx = sizeCanvas(canvas, MAX_COLUMNS);
    const { width, height } = canvas;
    const bins = FFT_SIZE / 2;
    const hop = Math.max(0, samples.length - FFT_SIZE) / width;
    // Average a few frames per column when columns are far apart, so short calls aren't missed
    const frames = Math.max(1, Math.min(4, Math.floor(hop / FFT_SIZE)));
    const re = new Float64Array(FFT_SIZE), im = new Float64Array(FFT_SIZE);
    const power = new Float32Array(width * bins);
    let max = -Infinity;

    for (let x = 0; x < width; x++) {
        const column = power.subarray(x * bins, (x + 1) * bins);
        for (let f = 0; f < frames; f++) {
            const start = Math.floor(x * hop + f * hop / frames);
            for (let i = 0; i < FFT_SIZE; i++) {
                re[i] = (samples[start + i] || 0) * HANN[i];
                im[i] = 0;
            }
            fft(re, im);
            for (let k = 0; k < bins; k++) column[k] += re[k] * re[k] + im[k] * im[k];
        }
        for (let k = 0; k < bins; k++) {
            column[k] = 10 * Math.log10(column[k] / frames + 1e-20);
            if (column[k] > max) max = column[k];
        }
    }

    const image = ctx.createImageData(width, height);
    for (let y = 0; y < height; y++) {
        const k = Math.min(bins - 1, Math.floor((1 - (y + 0.5) / height) * bins));
        for (let x = 0; x < width; x++) {
            const level = Math.max(0, Math.min(1, (power[x * bins + k] - max + DB_RANGE) / DB_RANGE));
            const c = Math.round(level * 255) * 3;
            const o = (y * width + x) * 4;
            image.data[o] = LUT[c];
            image.data[o + 1] = LUT[c + 1];
            image.data[o + 2] = LUT[c + 2];
            image.data[o + 3] = 255;
        }
    }
    ctx.putImageData(image, 0, 0);
}

const TONES = { anthro: '#fbbf24', bio: '#4ade80' };

function niceStep(span, maxTicks, steps) {
    return steps.find(s => span / s <= maxTicks) || steps[steps.length - 1];
}

// Axes, plus the bands an index uses ([{lo, hi, label, tone, dividers}] in Hz); the rest is dimmed
export function drawOverlay(canvas, { duration, nyquist, bands = [], caption = '' }) {
    const ctx = sizeCanvas(canvas);
    const dpr = window.devicePixelRatio || 1;
    ctx.scale(dpr, dpr);
    const w = canvas.width / dpr, h = canvas.height / dpr;
    const y = f => h * (1 - Math.max(0, Math.min(f, nyquist)) / nyquist);

    ctx.font = '11px system-ui, -apple-system, "Segoe UI", sans-serif';
    ctx.textBaseline = 'middle';
    ctx.shadowColor = 'rgba(0, 0, 0, 0.85)';

    if (bands.length) {
        ctx.fillStyle = 'rgba(0, 0, 0, 0.6)';
        ctx.fillRect(0, 0, w, h);
        for (const band of bands) ctx.clearRect(0, y(band.hi), w, y(band.lo) - y(band.hi));
        for (const band of bands) {
            const colour = TONES[band.tone] || '#ffffff';
            ctx.strokeStyle = colour;
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            for (const f of [band.lo, band.hi]) {
                ctx.moveTo(0, Math.round(y(f)) + 0.5);
                ctx.lineTo(w, Math.round(y(f)) + 0.5);
            }
            ctx.stroke();
            if (band.dividers) {
                ctx.setLineDash([3, 4]);
                ctx.globalAlpha = 0.55;
                ctx.beginPath();
                for (const f of band.dividers) {
                    ctx.moveTo(0, Math.round(y(f)) + 0.5);
                    ctx.lineTo(w, Math.round(y(f)) + 0.5);
                }
                ctx.stroke();
                ctx.setLineDash([]);
                ctx.globalAlpha = 1;
            }
            if (band.label) {
                ctx.fillStyle = colour;
                ctx.shadowBlur = 4;
                ctx.textAlign = 'right';
                ctx.fillText(band.label, w - 8, (y(band.lo) + y(band.hi)) / 2);
                ctx.shadowBlur = 0;
            }
        }
    }

    ctx.fillStyle = '#ffffff';
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.7)';
    ctx.lineWidth = 1;
    ctx.shadowBlur = 3;

    // Frequency ticks (kHz) up the left edge
    ctx.textAlign = 'left';
    const kHz = niceStep(nyquist / 1000, 6, [0.5, 1, 2, 5, 10, 20, 50, 100, 200]);
    for (let f = kHz; f * 1000 < nyquist * 0.97; f += kHz) {
        const ty = Math.round(y(f * 1000)) + 0.5;
        ctx.beginPath();
        ctx.moveTo(0, ty);
        ctx.lineTo(5, ty);
        ctx.stroke();
        ctx.fillText(`${f} kHz`, 8, ty);
    }

    // Time ticks (s) along the bottom
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    const secs = niceStep(duration, Math.max(2, Math.floor(w / 90)), [0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 30, 60, 120, 300, 600, 1800, 3600]);
    for (let t = secs; t < duration * 0.98; t += secs) {
        const tx = Math.round(t / duration * w) + 0.5;
        ctx.beginPath();
        ctx.moveTo(tx, h);
        ctx.lineTo(tx, h - 5);
        ctx.stroke();
        ctx.fillText(t >= 60 ? `${Math.floor(t / 60)}:${String(Math.round(t % 60)).padStart(2, '0')}` : `${+t.toFixed(1)} s`, tx, h - 7);
    }

    if (caption) {
        ctx.textAlign = 'left';
        ctx.textBaseline = 'top';
        ctx.font = '600 12px system-ui, -apple-system, "Segoe UI", sans-serif';
        ctx.shadowBlur = 4;
        ctx.fillText(caption, 60, 8);
    }
    ctx.shadowBlur = 0;
}
