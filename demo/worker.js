/*
    Ecoacoustic Vamp Plugins - browser demo

    Runs the vamp-wasm host (vamp-host.js/.wasm) in a Web Worker and loads
    vamp-ecoacoustics.wasm into it as an Emscripten side module. Working in a
    worker keeps the page responsive, and lets dlopen() compile the plugin
    synchronously (browsers only allow that for small modules on the main thread).

    Messages in:  {type: 'describe', gen, sampleRate}
                  {type: 'audio', gen, sampleRate, channels: Float32Array[]}
                  {type: 'run', gen, id, key, params: {identifier: value}}
    Messages out: {type: 'ready', plugins: [key]}
                  {type: 'described', gen, sampleRate, plugins: [info]}
                  {type: 'progress', id, fraction}
                  {type: 'result', id, key, features, blockSize, stepSize, elapsedMs}
                  {type: 'error', id, message}

    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

importScripts('vamp-host.js');

const LIBRARY = 'vamp-ecoacoustics';
const pluginKeys = [];
let Module = null;
let audio = null;

const ready = (async () => {
    Module = await VampHostModule();

    const response = await fetch(LIBRARY + '.wasm');
    if (!response.ok) throw new Error(`Could not fetch ${LIBRARY}.wasm (HTTP ${response.status})`);
    Module.FS.mkdir('/vamp');
    // PluginLoader only scans for the platform's library suffix, which is .so under Emscripten
    Module.FS.writeFile(`/vamp/${LIBRARY}.so`, new Uint8Array(await response.arrayBuffer()));

    const host = new Module.WasmVampHost(44100);
    host.setEnv('VAMP_PATH', '/vamp');
    const keys = host.listPlugins();
    for (let i = 0; i < keys.size(); i++) pluginKeys.push(keys.get(i));
    keys.delete();
    host.delete();
    if (!pluginKeys.length) throw new Error(`No plugins found in ${LIBRARY}.wasm`);
})();

ready.then(
    () => postMessage({ type: 'ready', plugins: pluginKeys }),
    err => postMessage({ type: 'error', message: err.message || String(err) })
);

function withPlugin(key, sampleRate, fn) {
    const host = new Module.WasmVampHost(sampleRate);
    try {
        if (!host.loadPlugin(key)) throw new Error(`Could not load ${key}`);
        return fn(host);
    } finally {
        host.delete();
    }
}

// Descriptors depend on the sample rate (e.g. frequency limits, ADI's block size)
function describe(sampleRate) {
    return pluginKeys.map(key => withPlugin(key, sampleRate, host => ({ key, ...host.getPluginInfo() })));
}

// Feed the audio block by block the way ReVAMP's runPlugin() does, so results
// match R: preferred block/step sizes, final partial block zero-padded.
function run(id, key, params) {
    const { sampleRate, channels } = audio;
    return withPlugin(key, sampleRate, host => {
        const info = host.getPluginInfo();
        let blockSize = info.preferredBlockSize || 1024;
        const stepSize = info.preferredStepSize || (info.inputDomain === 'frequency' ? blockSize / 2 : blockSize);
        if (stepSize > blockSize) blockSize = info.inputDomain === 'frequency' ? stepSize * 2 : stepSize;

        for (const [identifier, value] of Object.entries(params)) host.setParameter(identifier, value);
        if (!host.initialise(channels.length, stepSize, blockSize)) {
            throw new Error(`${key} rejected ${channels.length} channel(s) with block size ${blockSize}`);
        }

        const features = {};
        const collect = set => {
            if (!set) return;
            for (const output of Object.keys(set)) (features[output] ||= []).push(...set[output]);
        };

        const total = channels[0].length;
        const padded = channels.map(() => new Float32Array(blockSize));
        const started = performance.now();
        let lastProgress = started;
        let start = 0;
        do {
            const end = start + blockSize;
            const block = channels.map((data, c) => {
                if (end <= total) return data.subarray(start, end);
                padded[c].fill(0);
                padded[c].set(data.subarray(start, total));
                return padded[c];
            });
            collect(host.process(block, start / sampleRate));
            start += stepSize;

            const now = performance.now();
            if (now - lastProgress > 100) {
                postMessage({ type: 'progress', id, fraction: Math.min(1, start / total) });
                lastProgress = now;
            }
        } while (start - stepSize + blockSize < total);
        collect(host.getRemainingFeatures());

        return { key, features, blockSize, stepSize, elapsedMs: performance.now() - started };
    });
}

let audioGen = 0;

onmessage = async ({ data }) => {
    try {
        await ready;
        switch (data.type) {
        case 'describe':
            postMessage({ type: 'described', gen: data.gen, sampleRate: data.sampleRate, plugins: describe(data.sampleRate) });
            break;
        case 'audio':
            audio = { sampleRate: data.sampleRate, channels: data.channels };
            audioGen = data.gen;
            break;
        case 'run':
            if (data.gen !== audioGen) return; // queued before a newer file was loaded
            postMessage({ type: 'result', id: data.id, ...run(data.id, data.key, data.params) });
            break;
        }
    } catch (err) {
        postMessage({ type: 'error', id: data.id, key: data.key, message: err.message || String(err) });
    }
};
