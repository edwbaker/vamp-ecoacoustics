// Runs every plugin of a release's WebAssembly library through vamp-wasm in Node.js on the 91 pilot recordings of the
// Acoustic Index User's Guide, at its defaults and as ReVAMP's runPlugin does (samples as 32-bit floats, the plugin's
// preferred block and step size, the final block zero-padded, the first value of the output's first feature), and
// compares each value with the R function's value in the 0.3.0 benchmark (reference.csv), with the same tolerance as
// check_release.R: within 1e-4 of the R value, relative to it, plus 1e-6 absolute.
//   node check_wasm.js <vamp-wasm folder> <vamp-ecoacoustics.wasm> <recordings folder> <results csv> [label]
// The vamp-wasm folder holds a vamp-wasm release zip's files (vamp-host.js, vamp-host.wasm, vamp-multi-host.js).
// Exits 1 if a plugin is not found, a recording is not the benchmark's (md5) or a value does not pass.
const fs = require('fs'), path = require('path'), crypto = require('crypto');
const [, , hostDir, pluginWasm, wavDir, outCsv, labelArg] = process.argv;
const label = labelArg || `${process.platform} Node.js ${process.version}`;
const REL = 1e-4, ABS = 1e-6;
const VampHostModule = require(path.resolve(hostDir, 'vamp-host.js'));
const VampMultiHost = require(path.resolve(hostDir, 'vamp-multi-host.js'));

function readCsv(file) {
  const [head, ...lines] = fs.readFileSync(file, 'utf8').trim().split(/\r?\n/);
  const keys = head.split(',');
  return lines.map(l => Object.fromEntries(l.split(',').map((v, i) => [keys[i], v])));
}

function readWav(file) {
  const b = fs.readFileSync(file);
  if (b.toString('latin1', 0, 4) !== 'RIFF' || b.toString('latin1', 8, 12) !== 'WAVE') throw new Error('not WAV ' + file);
  let pos = 12, fmt = null, data = null;
  while (pos + 8 <= b.length) {
    const id = b.toString('latin1', pos, pos + 4); let size = b.readUInt32LE(pos + 4); const body = pos + 8;
    if (id === 'fmt ') {
      fmt = { format: b.readUInt16LE(body), channels: b.readUInt16LE(body + 2), rate: b.readUInt32LE(body + 4), bits: b.readUInt16LE(body + 14) };
      if (fmt.format === 0xFFFE && size >= 26) fmt.format = b.readUInt16LE(body + 24);
    } else if (id === 'data') { size = Math.min(size, b.length - body); data = [body, size]; break; }
    pos = body + size + (size & 1);
  }
  if (!fmt || !data) throw new Error('no fmt/data ' + file);
  const n = Math.floor(data[1] / (fmt.bits / 8)), x = new Float32Array(n), o = data[0];
  for (let i = 0; i < n; i++) {
    if (fmt.format === 1) {
      if (fmt.bits === 16) x[i] = b.readInt16LE(o + 2 * i) / 32768;
      else if (fmt.bits === 8) x[i] = (b[o + i] - 128) / 128;
      else if (fmt.bits === 24) x[i] = (b.readIntLE(o + 3 * i, 3)) / 8388608;
      else if (fmt.bits === 32) x[i] = b.readInt32LE(o + 4 * i) / 2147483648;
      else throw new Error('bits ' + fmt.bits);
    } else if (fmt.format === 3 && fmt.bits === 32) x[i] = b.readFloatLE(o + 4 * i);
    else throw new Error('format ' + fmt.format);
  }
  const frames = Math.floor(n / fmt.channels), chans = [];
  for (let c = 0; c < fmt.channels; c++) { const ch = new Float32Array(frames); for (let i = 0; i < frames; i++) ch[i] = x[i * fmt.channels + c]; chans.push(ch); }
  return { rate: fmt.rate, chans, frames };
}

function analyse(multi, w, plugin, outId) {
  const host = multi.createHost(w.rate);
  try {
    if (!host.loadPlugin('vamp-ecoacoustics:' + plugin)) throw new Error('load ' + plugin);
    const info = host.getPluginInfo();
    const block = info.preferredBlockSize, step = info.preferredStepSize || block;
    if (block !== step) throw new Error('block != step for ' + plugin);
    if (!host.initialise(w.chans.length, step, block)) throw new Error('initialise ' + plugin);
    const oi = info.outputs.findIndex(o => o.identifier === outId);
    let value = null;
    const take = fs_ => { const f = fs_ && fs_[oi]; if (value === null && f && f.length) { const v = Object.values(f[0].values); if (v.length) value = v[0]; } };
    let read = 0;
    do {
      const blk = w.chans.map(c => c.subarray(read, Math.min(read + block, w.frames)));
      take(host.process(blk, read / w.rate));
      read += block;
    } while (read < w.frames);
    take(host.getRemainingFeatures());
    return value;
  } finally { host.delete(); }
}

(async () => {
  const ref = readCsv(path.join(__dirname, 'reference.csv'));
  const libs = [{ name: 'vamp-ecoacoustics.so', data: new Uint8Array(fs.readFileSync(pluginWasm)) }];
  const multi = await VampMultiHost.create(VampHostModule, libs);
  const listed = multi.listPlugins();
  const missing = [...new Set(ref.map(r => 'vamp-ecoacoustics:' + r.plugin))].filter(k => !listed.includes(k));
  if (missing.length) throw new Error('not found by vamp-wasm: ' + missing.join(', '));
  console.log(`${label}\nNode.js ${process.version}, vamp-wasm from ${hostDir}`);

  const wavs = fs.readdirSync(wavDir, { recursive: true, withFileTypes: true })
    .filter(d => d.isFile() && /\.wav$/i.test(d.name)).map(d => path.join(d.parentPath, d.name));
  const md5 = new Map(wavs.map(f => [f, crypto.createHash('md5').update(fs.readFileSync(f)).digest('hex')]));
  const pathOf = (f, sum) => {
    const hit = wavs.find(p => path.basename(p) === f && md5.get(p) === sum);
    if (!hit) throw new Error(`${f} (md5 ${sum}) is not in ${wavDir}`);
    return hit;
  };

  const rows = [], byFile = new Map();
  for (const r of ref) { if (!byFile.has(r.file)) byFile.set(r.file, []); byFile.get(r.file).push(r); }
  for (const [f, rs] of byFile) {
    const w = readWav(pathOf(f, rs[0].md5));
    for (const r of rs) {
      let v = null;
      try { v = analyse(multi, w, r.plugin, r.output); } catch (e) { console.error(f, r.plugin, e.message); }
      const rv = Number(r.r_value), bv = Number(r.benchmark_value);
      const ok = v !== null && Number.isFinite(v) && Math.abs(v - rv) <= REL * Math.abs(rv) + ABS;
      rows.push({ host: label, file: f, plugin: r.plugin, output: r.output, value: v, r_value: rv, benchmark_value: bv,
                  diff_r: v === null ? null : v - rv, diff_benchmark: v === null ? null : v - bv, ok });
    }
  }
  const cols = Object.keys(rows[0]);
  const q = s => (typeof s === 'string' ? `"${s}"` : s === null ? 'NA' : String(s));
  fs.writeFileSync(outCsv, [cols.map(c => `"${c}"`).join(','), ...rows.map(r => cols.map(c => q(r[c])).join(','))].join('\n') + '\n');

  const plugins = [...new Set(ref.map(r => r.plugin))];
  const max = a => a.reduce((m, x) => (x > m ? x : m), 0);
  const s = plugins.map(p => {
    const x = rows.filter(r => r.plugin === p);
    return { plugin: p, recordings: x.length, failed: x.filter(r => !r.ok).length,
             max_abs_diff_r: max(x.map(r => Math.abs(r.diff_r ?? Infinity))),
             max_rel_diff_r: max(x.filter(r => r.r_value !== 0).map(r => Math.abs(r.diff_r ?? Infinity) / Math.abs(r.r_value))),
             max_abs_diff_benchmark: max(x.map(r => Math.abs(r.diff_benchmark ?? Infinity))) };
  });
  console.table(s.map(x => ({ ...x, max_abs_diff_r: x.max_abs_diff_r.toPrecision(3), max_rel_diff_r: x.max_rel_diff_r.toPrecision(3),
                              max_abs_diff_benchmark: x.max_abs_diff_benchmark.toPrecision(3) })));
  if (process.env.GITHUB_STEP_SUMMARY) {
    fs.appendFileSync(process.env.GITHUB_STEP_SUMMARY, [`### ${label}`, '', `Node.js ${process.version}`, '',
      '| plugin | recordings | failed | largest difference from R | largest relative difference from R | largest difference from the benchmark build |',
      '|---|---|---|---|---|---|',
      ...s.map(x => `| ${x.plugin} | ${x.recordings} | ${x.failed} | ${x.max_abs_diff_r.toPrecision(3)} | ${x.max_rel_diff_r.toPrecision(3)} | ${x.max_abs_diff_benchmark.toPrecision(3)} |`),
      ''].join('\n') + '\n');
  }
  const failed = rows.filter(r => !r.ok).length;
  if (failed) { console.log(`${failed} of ${rows.length} values outside the tolerance`); process.exit(1); }
  console.log(`all ${rows.length} values within the tolerance`);
})().catch(e => { console.error(e.stack || e); process.exit(1); });
