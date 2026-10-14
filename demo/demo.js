/*
    Ecoacoustic Vamp Plugins - browser demo

    Loads a recording (or synthesises one), hands it to the Vamp host running in
    worker.js, and shows the value of each index.

    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

import { decodeAudio, encodeWav, synthesize, EXAMPLES } from './audio.js';
import { drawSpectrogram, drawOverlay } from './spectrogram.js';

const $ = selector => document.querySelector(selector);

// Display order: complexity, diversity and evenness, bioacoustic, soundscape, entropy
const ORDER = ['aci-acc', 'aci', 'adi-acc', 'aei-acc', 'bi-acc', 'ndsi', 'h', 'sh', 'th'];

const ABOUT = {
    'aci-acc': 'How much intensity changes from one frame to the next within each frequency bin, relative to the total, summed across the recording. Modulated sound such as bird song raises it; steady noise does not.',
    'aci': 'The same measure computed as seewave does, over the whole file by default.',
    'adi-acc': 'Shannon diversity of the frequency bands louder than the dB threshold. Higher when sound fills many bands.',
    'aei-acc': 'Gini coefficient of the same bands. Higher when sound is concentrated in a few bands, lower when it is spread evenly.',
    'bi-acc': 'Area under the mean spectrum between the frequency limits, above its quietest point. Designed to track bird activity.',
    'ndsi': 'Biophony against anthrophony: +1 when all the energy is in the biophony band, −1 when it is all in the anthrophony band.',
    'h': 'Spectral entropy × temporal entropy. Near 1 for noise, near 0 for a pure tone.',
    'sh': 'Evenness of the mean spectrum. Near 1 when energy is spread across all frequencies.',
    'th': 'Evenness of the amplitude envelope over time. Near 1 for steady sound, lower for sparse or pulsed sound.',
};

// The R functions each plugin is validated against
const REFERENCE = {
    'aci-acc': 'soundecology::acoustic_complexity',
    'aci': 'seewave::ACI',
    'adi-acc': 'soundecology::acoustic_diversity',
    'aei-acc': 'soundecology::acoustic_evenness',
    'bi-acc': 'soundecology::bioacoustic_index',
    'ndsi': 'soundecology::ndsi',
    'h': 'seewave::H',
    'sh': 'seewave::sh',
    'th': 'seewave::th',
};

// The frequency ranges (Hz) an index uses, from its current parameter values
const wholeRange = (hz, nyquist) => [{ lo: hz('minFreq'), hi: hz('maxFreq') || nyquist }];
const steppedRange = hz => {
    const max = hz('maxFreq'), step = hz('freqStep'), dividers = [];
    for (let f = step; f < max; f += step) dividers.push(f);
    return [{ lo: 0, hi: max, dividers, label: `${dividers.length + 1} bands of ${formatKHz(step)} kHz` }];
};
const BANDS = {
    'aci-acc': wholeRange,
    'aci': wholeRange,
    'adi-acc': steppedRange,
    'aei-acc': steppedRange,
    'bi-acc': hz => [{ lo: hz('minFreq'), hi: hz('maxFreq') }],
    'ndsi': hz => [
        { lo: hz('anthroMin'), hi: hz('anthroMax'), label: 'anthrophony', tone: 'anthro' },
        { lo: hz('bioMin'), hi: hz('bioMax'), label: 'biophony', tone: 'bio' },
    ],
};

const state = {
    plugins: [],     // plugin descriptions from the worker, in display order
    overrides: {},   // plugin key -> {parameter: value} the user has changed
    audio: null,     // {name, sampleRate, channels, duration, url, download, note}
    gen: 0,          // bumped for each new recording
    runs: {},        // plugin key -> id of its latest run
    results: {},     // plugin key -> result or {error}
    highlight: null, // plugin key whose bands are drawn on the spectrogram
    ready: false,    // the worker has loaded the host and plugins
};
let nextRunId = 1;

// ---- Status -----------------------------------------------------------------

function status(text, kind = 'info') {
    const el = $('#status');
    el.textContent = text;
    el.dataset.kind = kind;
}

// ---- Worker -----------------------------------------------------------------

let worker;
if (location.protocol === 'file:') {
    status('Browsers won’t run WebAssembly workers from file:// pages. Serve this folder over HTTP, e.g. python -m http.server 8000 -d demo', 'error');
} else {
    worker = new Worker('worker.js');
    worker.onerror = event => {
        event.preventDefault();
        status(`The Vamp host failed to start (${event.message || 'see the console'}). Has demo/build.sh been run?`, 'error');
    };
    worker.onmessage = ({ data }) => handlers[data.type](data);
}

const handlers = {
    ready({ plugins }) {
        state.ready = true;
        status(`Ready: ${plugins.length} indices loaded from vamp-ecoacoustics.wasm.`, 'ready');
        for (const el of document.querySelectorAll('[data-needs-host]')) el.disabled = false;
        worker.postMessage({ type: 'describe', gen: state.gen, sampleRate: state.audio?.sampleRate || 44100 });
    },

    described({ gen, plugins }) {
        if (gen !== state.gen) return;
        renderCards(plugins);
        if (state.audio) runAll();
    },

    progress({ id, fraction }) {
        const key = Object.keys(state.runs).find(k => state.runs[k] === id);
        if (key) card(key)?.style.setProperty('--progress', fraction);
    },

    result(result) {
        if (state.runs[result.key] !== result.id) return; // superseded
        state.results[result.key] = result;
        showResult(result.key);
        updateSummary();
    },

    error({ id, key, message }) {
        if (key && state.runs[key] === id) {
            state.results[key] = { error: message };
            showResult(key);
            updateSummary();
        } else {
            status(message, 'error');
        }
    },
};

// ---- Loading audio ----------------------------------------------------------

async function loadFile(file) {
    if (!file) return;
    status(`Reading ${file.name}…`);
    try {
        const audio = await decodeAudio(await file.arrayBuffer(), file.name);
        setAudio({ name: file.name, url: URL.createObjectURL(file), ...audio });
    } catch (err) {
        status(`Could not read ${file.name}: ${err.message}`, 'error');
    }
}

function loadExample(kind) {
    const example = EXAMPLES[kind], sampleRate = 44100;
    const channels = [synthesize(kind, sampleRate)];
    const url = URL.createObjectURL(encodeWav(channels, sampleRate));
    setAudio({ name: example.file, sampleRate, channels, url, download: example.file, note: 'Synthetic example' });
}

function setAudio(audio) {
    if (state.audio) URL.revokeObjectURL(state.audio.url);
    audio.duration = audio.channels[0].length / audio.sampleRate;
    state.audio = audio;
    state.gen++;
    state.runs = {};
    state.results = {};

    $('#recording').hidden = false;
    $('#player').src = audio.url;
    const save = $('#save');
    save.hidden = !audio.download;
    if (audio.download) {
        save.href = audio.url;
        save.download = audio.download;
    }
    const channels = audio.channels.length === 1 ? 'mono' : audio.channels.length === 2 ? 'stereo' : `${audio.channels.length} channels`;
    $('#file-name').textContent = audio.name;
    $('#file-details').textContent = [formatDuration(audio.duration), `${+(audio.sampleRate / 1000).toFixed(3)} kHz`, channels, audio.note]
        .filter(Boolean).join(' · ');
    drawAll();

    status(`Computing indices for ${audio.name}…`);
    worker.postMessage({ type: 'audio', gen: state.gen, sampleRate: audio.sampleRate, channels: audio.channels });
    // Parameter ranges and some block sizes depend on the sample rate, so describe again before running
    worker.postMessage({ type: 'describe', gen: state.gen, sampleRate: audio.sampleRate });
}

// ---- Running ----------------------------------------------------------------

function runAll() {
    for (const plugin of state.plugins) run(plugin.key);
    updateSummary();
}

function run(key) {
    if (!state.audio) return;
    const id = nextRunId++;
    state.runs[key] = id;
    delete state.results[key];
    const el = card(key);
    el.classList.add('running');
    el.style.setProperty('--progress', 0);
    worker.postMessage({ type: 'run', gen: state.gen, id, key, params: state.overrides[key] || {} });
}

// ---- Cards ------------------------------------------------------------------

const card = key => document.querySelector(`.card[data-key="${CSS.escape(key)}"]`);
const pluginFor = key => state.plugins.find(p => p.key === key);

function renderCards(plugins) {
    const rank = p => (ORDER.indexOf(p.identifier) + 1 || ORDER.length + 1);
    state.plugins = [...plugins].sort((a, b) => rank(a) - rank(b));
    $('#results').replaceChildren(...state.plugins.map(buildCard));
    drawBands();
}

function buildCard(plugin) {
    const el = $('#card-template').content.firstElementChild.cloneNode(true);
    el.dataset.key = plugin.key;
    el.querySelector('h3').textContent = plugin.name;
    el.querySelector('.id').textContent = plugin.identifier;
    el.querySelector('.about').textContent = ABOUT[plugin.identifier] || plugin.description;
    el.querySelector('.reference').textContent = REFERENCE[plugin.identifier] ? `Matches ${REFERENCE[plugin.identifier]}()` : '';
    if (plugin.parameters.length) buildParams(el, plugin);
    else el.querySelector('.params').remove();

    el.addEventListener('mouseenter', () => highlight(plugin.key));
    el.addEventListener('mouseleave', () => highlight(null));
    el.addEventListener('focusin', () => highlight(plugin.key));
    return el;
}

function buildParams(el, plugin) {
    const fields = el.querySelector('.fields');
    for (const param of plugin.parameters) {
        const label = document.createElement('label');
        label.title = param.description;
        const name = document.createElement('span');
        name.textContent = param.name;
        const input = document.createElement('input');
        input.type = 'number';
        input.min = param.minValue;
        input.max = param.maxValue;
        input.step = param.isQuantized ? param.quantizeStep : 'any';
        input.value = round(valueOf(plugin, param.identifier));
        input.dataset.param = param.identifier;
        input.disabled = !worker;
        const unit = document.createElement('small');
        unit.textContent = param.unit;
        label.append(name, input, unit);
        fields.append(label);
    }
    fields.addEventListener('change', ({ target }) => {
        const param = plugin.parameters.find(p => p.identifier === target.dataset.param);
        let value = Number(target.value);
        if (target.value === '' || !Number.isFinite(value)) value = param.defaultValue;
        value = Math.min(param.maxValue, Math.max(param.minValue, value));
        target.value = round(value);
        setOverride(plugin, param.identifier, value === param.defaultValue ? undefined : value);
    });
    el.querySelector('.reset').addEventListener('click', () => {
        delete state.overrides[plugin.key];
        for (const input of fields.querySelectorAll('input')) {
            input.value = round(valueOf(plugin, input.dataset.param));
        }
        setOverride(plugin);
    });
    updateParamSummary(el, plugin);
}

function setOverride(plugin, identifier, value) {
    if (identifier) {
        const overrides = state.overrides[plugin.key] ||= {};
        if (value === undefined) delete overrides[identifier];
        else overrides[identifier] = value;
        if (!Object.keys(overrides).length) delete state.overrides[plugin.key];
    }
    updateParamSummary(card(plugin.key), plugin);
    if (state.highlight === plugin.key) drawBands();
    run(plugin.key);
}

function updateParamSummary(el, plugin) {
    const changed = Object.keys(state.overrides[plugin.key] || {}).length;
    const summary = el.querySelector('.params summary');
    if (!summary) return;
    summary.textContent = `Parameters (${plugin.parameters.length})${changed ? ` · ${changed} changed` : ''}`;
    el.querySelector('.reset').hidden = !changed;
}

function valueOf(plugin, identifier) {
    const override = state.overrides[plugin.key]?.[identifier];
    return override ?? plugin.parameters.find(p => p.identifier === identifier)?.defaultValue;
}

// One row per output, pairing "x" with "x_right" and ignoring right-channel outputs for mono audio
function outputRows(plugin, features) {
    const stereo = state.audio.channels.length > 1;
    const rows = [];
    plugin.outputs.forEach((output, i) => {
        const values = (features[i] || []).filter(f => f.values.length).map(f => f.values[0]);
        if (!values.length) return;
        const isRight = output.identifier.endsWith('_right');
        if (isRight && !stereo) return;
        const left = isRight && rows.find(r => r.id === output.identifier.slice(0, -'_right'.length));
        if (left) left.values.push(...values);
        else rows.push({ id: output.identifier, name: output.name, values });
    });
    return rows;
}

function showResult(key) {
    const el = card(key), plugin = pluginFor(key), result = state.results[key];
    if (!el || !plugin) return;
    el.classList.remove('running');
    const value = el.querySelector('.value'), outputs = el.querySelector('.outputs'), meta = el.querySelector('.meta');
    outputs.replaceChildren();

    if (result.error) {
        value.replaceChildren(Object.assign(document.createElement('span'), { className: 'error', textContent: result.error }));
        meta.textContent = '';
        return;
    }

    const [main, ...rest] = outputRows(plugin, result.features);
    value.replaceChildren(...(main ? channelValues(main.values) : [Object.assign(document.createElement('span'), { className: 'error', textContent: 'No output' })]));
    for (const row of rest) {
        const dt = document.createElement('dt');
        dt.textContent = row.name.replace(/ \(Right\)$/, '');
        const dd = document.createElement('dd');
        dd.append(...channelValues(row.values));
        outputs.append(dt, dd);
    }
    meta.textContent = `${result.blockSize}-sample blocks · ${formatMs(result.elapsedMs)}`;
}

function channelValues(values) {
    return values.map((v, c) => {
        const span = document.createElement('span');
        if (values.length > 1) {
            const tag = document.createElement('small');
            tag.textContent = values.length === 2 ? 'LR'[c] : `${c + 1}`;
            span.append(tag, ' ');
        }
        span.append(formatValue(v));
        return span;
    });
}

// ---- Summary and export -----------------------------------------------------

function updateSummary() {
    const bar = $('#summary');
    if (!state.audio || !state.plugins.length) {
        bar.hidden = true;
        return;
    }
    bar.hidden = false;
    const done = state.plugins.filter(p => state.results[p.key]);
    $('#download').disabled = done.length < state.plugins.length;
    if (done.length < state.plugins.length) {
        $('#summary-text').textContent = `Computing… ${done.length} of ${state.plugins.length} indices`;
        return;
    }
    const failed = done.filter(p => state.results[p.key].error).length;
    const ms = done.reduce((total, p) => total + (state.results[p.key].elapsedMs || 0), 0);
    const speed = state.audio.duration * 1000 / ms;
    $('#summary-text').textContent = `${state.plugins.length - failed} indices in ${formatMs(ms)}` +
        (Number.isFinite(speed) ? `, ${Math.round(speed).toLocaleString()}× faster than real time` : '') +
        (failed ? ` (${failed} failed)` : '');
    status(`Done: ${state.audio.name}`, 'ready');
}

function csv() {
    const quote = v => /[",\n]/.test(v) ? `"${String(v).replace(/"/g, '""')}"` : String(v);
    const lines = [['file', 'index', 'output', 'channel', 'value']];
    for (const plugin of state.plugins) {
        const result = state.results[plugin.key];
        if (!result?.features) continue;
        for (const row of outputRows(plugin, result.features)) {
            row.values.forEach((v, c) => lines.push([state.audio.name, plugin.identifier, row.id, c + 1, v]));
        }
    }
    return lines.map(line => line.map(quote).join(',')).join('\n') + '\n';
}

// ---- Spectrogram ------------------------------------------------------------

function drawAll() {
    if (!state.audio) return;
    drawSpectrogram($('#spectrogram'), state.audio.channels[0]);
    drawBands();
}

function highlight(key) {
    if (state.highlight === key) return;
    state.highlight = key;
    drawBands();
}

function drawBands() {
    if (!state.audio) return;
    const nyquist = state.audio.sampleRate / 2;
    const plugin = state.highlight && pluginFor(state.highlight);
    let bands = [], caption = '';
    if (plugin) {
        const hz = id => {
            const param = plugin.parameters.find(p => p.identifier === id);
            return valueOf(plugin, id) * (param?.unit === 'kHz' ? 1000 : 1);
        };
        bands = (BANDS[plugin.identifier]?.(hz, nyquist) || []).map(b => ({ ...b, hi: Math.min(b.hi, nyquist) }));
        const range = b => `${formatKHz(b.lo)}–${formatKHz(b.hi)} kHz`;
        caption = `${plugin.identifier.toUpperCase()}: ` + (bands.length
            ? bands.map(b => b.dividers ? `${b.label}, ${range(b)}` : `${b.label ? b.label + ' ' : ''}${range(b)}`).join(', ')
            : 'whole spectrum');
    }
    drawOverlay($('#overlay'), { duration: state.audio.duration, nyquist, bands, caption });
}

// ---- Formatting -------------------------------------------------------------

function formatValue(v) {
    if (!Number.isFinite(v)) return String(v);
    const a = Math.abs(v);
    return a >= 1000 ? v.toFixed(1) : a >= 1 ? v.toFixed(3) : v.toFixed(4);
}
const round = v => +v.toPrecision(6);
const formatKHz = hz => `${+(hz / 1000).toFixed(2)}`;
const formatMs = ms => ms < 1000 ? `${Math.max(1, Math.round(ms))} ms` : `${(ms / 1000).toFixed(2)} s`;
function formatDuration(s) {
    return s < 60 ? `${s.toFixed(1)} s` : `${Math.floor(s / 60)} min ${Math.round(s % 60)} s`;
}

// ---- Wiring -----------------------------------------------------------------

$('#file').addEventListener('change', e => loadFile(e.target.files[0]));
const drop = $('#drop');
drop.addEventListener('dragover', e => {
    e.preventDefault();
    drop.classList.add('over');
});
drop.addEventListener('dragleave', () => drop.classList.remove('over'));
drop.addEventListener('drop', e => {
    e.preventDefault();
    drop.classList.remove('over');
    if (state.ready && e.dataTransfer.files.length) loadFile(e.dataTransfer.files[0]);
});
for (const button of document.querySelectorAll('[data-example]')) {
    button.addEventListener('click', () => loadExample(button.dataset.example));
}
$('#download').addEventListener('click', () => {
    const url = URL.createObjectURL(new Blob([csv()], { type: 'text/csv' }));
    const name = `${state.audio.name.replace(/\.[^.]+$/, '')}-indices.csv`;
    Object.assign(document.createElement('a'), { href: url, download: name }).click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
});

let resizeTimer;
new ResizeObserver(() => {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(drawAll, 150);
}).observe($('.spectrogram'));
