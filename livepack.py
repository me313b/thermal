"""Live Pack v2: an interactive instrument, not a decoration.

Hover any region for its live numbers; click to open a station inspector;
flow arrows, a temperature legend, and a pulsing ring on the weakest link.
Self-diagnosing (errors paint themselves on the canvas), retina-safe,
badges pre-filled server-side. All state injected from the solver.
"""
import json


def live_pack_html(state: dict, stations: list) -> str:
    s = json.dumps(state)
    st_js = json.dumps(stations)
    S = state
    b_q = f"<b>{S['q_kw']:.2f} kW</b> at {S['c_rms']:.2f}C rms"
    b_t = (f"can <b>{S['T_b']:.1f}°C</b> · core <b>{S['T_core']:.1f}°C</b>"
           f"<br>limit {S['T_limit']:.0f}°C")
    b_w = (f"water {S['T_w_in']:.0f}→{S['T_w_in']+S['dT_water']:.1f}°C · "
           f"{S['flow_lpm']:.0f} L/min")
    b_u = (f"{'guided' if S['mode']=='serpentine' else 'oil'} "
           f"<b>{S['u_mm_s']:.1f} mm/s</b> · spread {S['spread']:.1f} °C")
    tpl = """
<div id="lp-root" style="font-family:Inter,-apple-system,'Segoe UI',sans-serif">
<style>
 #lp-root{position:relative;background:#0B1220;border-radius:16px;
   overflow:hidden;box-shadow:0 1px 3px rgba(16,24,40,.12)}
 #lp-cv{width:100%;height:430px;display:block;cursor:crosshair}
 #lp-bar{position:absolute;left:12px;top:10px;display:flex;gap:6px;
   align-items:center;background:rgba(15,23,42,.6);backdrop-filter:blur(8px);
   border:1px solid rgba(255,255,255,.14);border-radius:12px;
   padding:5px 9px;z-index:4}
 .lpb{border:0;background:transparent;color:#CBD5E1;font-weight:600;
   font-size:11.5px;padding:4px 9px;border-radius:8px;cursor:pointer;
   transition:all .15s}
 .lpb:hover{background:rgba(255,255,255,.10);color:#fff}
 .lpb.on{background:rgba(255,255,255,.16);color:#fff}
 #lp-spd{width:76px;accent-color:#8B9CF9}
 #lp-status{position:absolute;right:14px;top:14px;font-size:11px;
   color:#94A3B8;z-index:4}
 .lp-badge{position:absolute;background:rgba(15,23,42,.6);
   backdrop-filter:blur(8px);border:1px solid rgba(255,255,255,.14);
   border-radius:10px;padding:5px 11px;font-size:11.5px;color:#CBD5E1;
   pointer-events:none;z-index:3}
 .lp-badge b{font-size:13px;color:#fff}
 #lp-tip{position:absolute;display:none;background:rgba(15,23,42,.92);
   border:1px solid rgba(255,255,255,.2);border-radius:8px;color:#E2E8F0;
   font-size:11px;padding:5px 9px;pointer-events:none;z-index:6;
   max-width:230px}
 #lp-insp{position:absolute;right:12px;top:52px;width:290px;max-height:330px;
   overflow-y:auto;background:rgba(15,23,42,.92);backdrop-filter:blur(10px);
   border:1px solid rgba(255,255,255,.2);border-radius:14px;color:#E2E8F0;
   font-size:12px;padding:12px 14px;display:none;z-index:5}
 #lp-insp h4{margin:0 0 6px 0;font-size:13px;color:#fff}
 #lp-insp .x{position:absolute;right:9px;top:7px;cursor:pointer;
   color:#94A3B8;font-size:15px}
 #lp-insp .x:hover{color:#fff}
 #lp-insp p{margin:5px 0;line-height:1.45;color:#CBD5E1}
 #lp-insp b{color:#fff}
 #lp-help{position:absolute;left:12px;bottom:52px;max-width:330px;
   background:rgba(15,23,42,.92);border:1px solid rgba(255,255,255,.2);
   border-radius:12px;color:#CBD5E1;font-size:11.5px;padding:10px 13px;
   display:none;z-index:5;line-height:1.5}
</style>
<canvas id="lp-cv"></canvas>
<div id="lp-bar">
  <button class="lpb on" id="lp-play">&#10074;&#10074;</button>
  <input type="range" id="lp-spd" min="0.2" max="4" step="0.1" value="1">
  <button class="lpb on" id="lp-oil">Oil</button>
  <button class="lpb on" id="lp-wat">Water</button>
  <button class="lpb on" id="lp-heat">Glow</button>
  <button class="lpb" id="lp-arr">Arrows</button>
  <button class="lpb" id="lp-exag">Exagg.</button>
  <button class="lpb" id="lp-hlp">?</button>
</div>
<div id="lp-status">live · hover and click the hardware</div>
<div class="lp-badge" id="bq" style="left:12px;bottom:12px">__BQ__</div>
<div class="lp-badge" id="bt" style="right:12px;bottom:12px;text-align:right">__BT__</div>
<div class="lp-badge" id="bw" style="right:12px;bottom:52px;text-align:right">__BW__</div>
<div class="lp-badge" id="bu" style="left:12px;bottom:52px">__BU__</div>
<div id="lp-tip"></div>
<div id="lp-insp"><span class="x" id="lp-x">&times;</span>
<div id="lp-insp-body"></div></div>
<div id="lp-help"><b>How to use this.</b> Everything is live from the
solver. <br>&bull; <b>Hover</b> a cell, tube, the oil, the headspace or the
chiller for its numbers.<br>&bull; <b>Click</b> to pin a full inspector for
that station.<br>&bull; The <b>red pulsing ring</b> marks the weakest link
in the chain right now.<br>&bull; <b>Arrows</b> shows the oil flow field;
<b>Exagg.</b> stretches the colour scale to the pack's own range so small
differences show.<br>Colours are real temperatures - see the scale at the
bottom.</div>
</div>
<script>
(function(){
"use strict";
const status = document.getElementById('lp-status');
function fail(msg){ status.textContent = 'error: ' + msg;
  status.style.color = '#FCA5A5';
  try { const c = document.getElementById('lp-cv').getContext('2d');
    c.setTransform(1,0,0,1,0,0); c.fillStyle='#FCA5A5';
    c.font='13px Inter,sans-serif';
    c.fillText('Live Pack error: '+msg, 14, 26);
  } catch(e){} }
window.addEventListener('error', e => fail(e.message));
try {
const S = __STATE__;
const ST = __STATIONS__;
const cv = document.getElementById('lp-cv');
const cx = cv.getContext('2d');
const tip = document.getElementById('lp-tip');
const insp = document.getElementById('lp-insp');
const inspB = document.getElementById('lp-insp-body');
let W = 900, H = 430;
function fit(){
  W = cv.clientWidth || 900; H = cv.clientHeight || 430;
  const d = window.devicePixelRatio || 1;
  cv.width = Math.round(W * d); cv.height = Math.round(H * d);
  cx.setTransform(d, 0, 0, d, 0, 0);
}
fit();
if (window.ResizeObserver) new ResizeObserver(fit).observe(cv);

let play=true, spd=1, showOil=true, showWat=true, showHeat=true,
    showArr=false, exag=false, sel=null, hover=null, mx=-1, my=-1, tsec=0;
function tog(id, fn){
  const b = document.getElementById(id);
  b.addEventListener('click', function(){
    b.classList.toggle('on'); fn(b.classList.contains('on'), b);
  });
}
tog('lp-play', (v,b)=>{play=v; b.innerHTML=v?'&#10074;&#10074;':'&#9654;';
  status.textContent=v?'live · hover and click the hardware':'paused';});
tog('lp-oil', v=>showOil=v); tog('lp-wat', v=>showWat=v);
tog('lp-heat', v=>showHeat=v); tog('lp-arr', v=>showArr=v);
tog('lp-exag', v=>exag=v);
tog('lp-hlp', v=>{document.getElementById('lp-help').style.display=
  v?'block':'none';});
document.getElementById('lp-spd').addEventListener('input',
  e=>spd=parseFloat(e.target.value));
document.getElementById('lp-x').addEventListener('click',
  ()=>{sel=null; insp.style.display='none';});

function tRange(){
  const lo = exag ? Math.min(S.T_w_in, S.T_amb) : 15;
  const hi = exag ? Math.max(S.T_core, S.T_b + 1) : 60;
  return [lo, hi];
}
function tcol(T, a){
  const [lo, hi] = tRange();
  let f = Math.min(Math.max((T - lo) / (hi - lo), 0), 1);
  const r = Math.round(80 + f * (241 - 80)),
        g = Math.round(150 + f * (82 - 150)),
        b = Math.round(240 + f * (82 - 240));
  return 'rgba(' + r + ',' + g + ',' + b + ',' + a + ')';
}
function rrect(x, y, w, h, r, fill, stroke){
  cx.beginPath();
  cx.moveTo(x + r, y);
  cx.arcTo(x + w, y, x + w, y + h, r);
  cx.arcTo(x + w, y + h, x, y + h, r);
  cx.arcTo(x, y + h, x, y, r);
  cx.arcTo(x, y, x + w, y, r);
  cx.closePath();
  if (fill) cx.fill(); if (stroke) cx.stroke();
}

// ---- geometry (recomputed each frame; regions kept for hit tests) ----
let R = {};
function geom(){
  const mL = 24, mR = 92, top = 58, bot = 34;
  const bx = mL, by = top, bw = W - mL - mR, bh = H - top - bot;
  const oilTop = by + bh * (1 - S.fill_frac);
  const cTop = by + bh * (1 - S.cell_top_frac);
  const cBot = by + bh * (1 - S.cell_bot_frac);
  const tubeY = S.interstitial ? (cTop + cBot) / 2
                               : Math.max(oilTop + 22, by + 26);
  R = {bx, by, bw, bh, oilTop, cTop, cBot, tubeY,
       n: S.n_rows_draw, nt: S.n_tubes_draw,
       chx: W - mR + 14, chy: by + 8, chw: mR - 26, chh: 74};
  R.pitch = bw / (R.n + 0.6);
  R.cw = R.pitch * S.d_over_p;
  R.tp = bw / (R.nt + 1);
}

// ---- flow field: k counter-rotating buoyant rolls in the oil ----
function field(px, py){                 // px,py normalised in oil region
  const u = S.u_mm_s * 4.0 * spd;
  if (S.mode === 'serpentine'){
    const lane = Math.floor(py * R.n);
    return {vx: (lane % 2 ? -1 : 1) * u * 2.4,
            vy: Math.sin(px * 14 + lane * 2) * 2.5};
  }
  const k = S.mode === 'stirred' ? 1 : 3;
  const A = u * 2.2;
  return {vx:  A * Math.sin(k * Math.PI * px) * Math.cos(Math.PI * py),
          vy: -A * Math.cos(k * Math.PI * px) * Math.sin(Math.PI * py)};
}

const NP = 170, PW = 80;
const oilP = [], watP = [];
for (let i = 0; i < NP; i++)
  oilP.push({x: Math.random(), y: Math.random()});
for (let i = 0; i < PW; i++)
  watP.push({x: Math.random(), lane: i % S.n_tubes_draw});

function hitTest(x, y){
  if (x >= R.chx && x <= R.chx + R.chw && y >= R.chy &&
      y <= R.chy + R.chh) return 'chiller';
  for (let j = 0; j < R.nt; j++){
    const tx = R.bx + R.tp * (j + 1);
    if (Math.hypot(x - tx, y - R.tubeY) < 15) return 'tubes';
  }
  if (y >= R.cTop && y <= R.cBot){
    for (let i = 0; i < R.n; i++){
      const cxl = R.bx + R.pitch * (0.4 + i) + (R.pitch - R.cw) / 2;
      if (x >= cxl && x <= cxl + R.cw) return 'cells';
    }
    if (x > R.bx && x < R.bx + R.bw) return 'oil';
  }
  if (x > R.bx && x < R.bx + R.bw && y > R.oilTop &&
      y < R.by + R.bh) return 'oil';
  if (x > R.bx && x < R.bx + R.bw && y > R.by && y < R.oilTop)
    return 'head';
  return null;
}
cv.addEventListener('mousemove', e=>{
  const r = cv.getBoundingClientRect();
  mx = e.clientX - r.left; my = e.clientY - r.top;
  hover = hitTest(mx, my);
  if (hover){
    const st = ST.find(s => s.id === hover);
    tip.innerHTML = '<b>' + st.title + '</b><br>' + st.hint
                    + '<br><span style="color:#94A3B8">click to inspect'
                    + '</span>';
    tip.style.display = 'block';
    tip.style.left = Math.min(mx + 14, W - 240) + 'px';
    tip.style.top = Math.min(my + 12, H - 70) + 'px';
    cv.style.cursor = 'pointer';
  } else { tip.style.display = 'none'; cv.style.cursor = 'crosshair'; }
});
cv.addEventListener('mouseleave', ()=>{ tip.style.display='none';
  hover=null; });
cv.addEventListener('click', ()=>{
  if (!hover) { sel=null; insp.style.display='none'; return; }
  sel = hover;
  const st = ST.find(s => s.id === sel);
  inspB.innerHTML = '<h4>' + st.title + '</h4>' + st.body;
  insp.style.display = 'block';
});

let last = performance.now();
function frame(now){
  try {
    const dt = Math.min((now - last) / 1000, 0.05); last = now;
    tsec += dt;
    geom();
    const {bx, by, bw, bh, oilTop, cTop, cBot, tubeY, n, nt, pitch, cw,
           tp} = R;
    cx.fillStyle = '#0B1220'; cx.fillRect(0, 0, W, H);
    // enclosure + headspace + oil
    cx.strokeStyle = 'rgba(255,255,255,.5)'; cx.lineWidth = 2;
    cx.strokeRect(bx - 7, by - 7, bw + 14, bh + 14);
    cx.fillStyle = '#0E1830'; cx.fillRect(bx, by, bw, oilTop - by);
    cx.fillStyle = 'rgba(245,158,11,0.10)';
    cx.fillRect(bx, oilTop, bw, by + bh - oilTop);
    cx.strokeStyle = 'rgba(245,158,11,.6)'; cx.lineWidth = 1;
    cx.beginPath(); cx.moveTo(bx, oilTop);
    cx.lineTo(bx + bw, oilTop); cx.stroke();
    cx.fillStyle = '#64748B'; cx.font = '10px Inter,sans-serif';
    cx.fillText('headspace (' + S.headspace_mm.toFixed(0) + ' mm N\\u2082)',
                bx + 8, by + 13);
    // cells
    for (let i = 0; i < n; i++){
      const x = bx + pitch * (0.4 + i) + (pitch - cw) / 2;
      if (showHeat){
        const gcx = x + cw / 2, gcy = (cTop + cBot) / 2;
        const gl = cx.createRadialGradient(gcx, gcy, 2, gcx, gcy, cw * 1.7);
        const a = Math.min(0.12 + S.q_kw * 0.05, 0.4);
        gl.addColorStop(0, 'rgba(241,82,82,' + a + ')');
        gl.addColorStop(1, 'rgba(241,82,82,0)');
        cx.fillStyle = gl;
        cx.fillRect(x - cw, cTop - cw, cw * 3, (cBot - cTop) + 2 * cw);
      }
      cx.fillStyle = tcol(S.T_b, 0.96);
      cx.strokeStyle = 'rgba(255,255,255,.30)'; cx.lineWidth = 1;
      rrect(x, cTop, cw, cBot - cTop, 5, true, true);
      cx.fillStyle = tcol(S.T_core, 0.95);
      rrect(x + cw * 0.30, cTop + 4, cw * 0.40, (cBot - cTop) - 8, 4,
            true, false);
      if (S.mode === 'serpentine' && i < n - 1){
        cx.fillStyle = 'rgba(139,156,249,.9)';
        cx.fillRect(x + cw + (pitch - cw) / 2 - 1.2, cTop, 2.4,
                    cBot - cTop);
      }
    }
    // flow arrows
    if (showArr){
      cx.strokeStyle = 'rgba(148,163,184,.7)'; cx.lineWidth = 1.2;
      for (let gi = 1; gi < 13; gi++) for (let gj = 1; gj < 6; gj++){
        const px = gi / 13, py = gj / 6;
        const ax = bx + px * bw,
              ay = oilTop + py * (by + bh - oilTop);
        const v = field(px, py);
        const L = Math.min(Math.hypot(v.vx, v.vy) * 0.35, 13);
        if (L < 2) continue;
        const th = Math.atan2(v.vy, v.vx);
        cx.beginPath(); cx.moveTo(ax, ay);
        cx.lineTo(ax + L * Math.cos(th), ay + L * Math.sin(th));
        cx.stroke();
        cx.beginPath();
        cx.arc(ax + L * Math.cos(th), ay + L * Math.sin(th), 1.6, 0,
               6.283);
        cx.fillStyle = 'rgba(148,163,184,.8)'; cx.fill();
      }
    }
    // tubes + water beads + inlet/outlet labels
    for (let j = 0; j < nt; j++){
      const tx = bx + tp * (j + 1);
      cx.fillStyle = 'rgba(165,180,204,.4)';
      cx.beginPath(); cx.arc(tx, tubeY, 11, 0, 6.283); cx.fill();
      cx.fillStyle = '#D97706';
      cx.beginPath(); cx.arc(tx, tubeY, 5.5, 0, 6.283); cx.fill();
    }
    cx.fillStyle = '#7DD3FC'; cx.font = '10px Inter,sans-serif';
    cx.fillText('water in ' + S.T_w_in.toFixed(0) + '°C', bx + 4,
                tubeY - 16);
    cx.fillStyle = tcol(S.T_w_in + S.dT_water, 1);
    cx.fillText('out ' + (S.T_w_in + S.dT_water).toFixed(1) + '°C',
                bx + bw - 74, tubeY - 16);
    if (showWat){
      for (const p of watP){
        if (play){ p.x += dt * S.flow_norm * 0.25 * spd;
                   if (p.x > 1) p.x -= 1; }
        const tx = bx + tp * (p.lane + 1);
        const ang = p.x * 6.283;
        cx.fillStyle = tcol(S.T_w_in + S.dT_water * p.x, 0.95);
        cx.beginPath();
        cx.arc(tx + Math.cos(ang) * 3.4, tubeY + Math.sin(ang) * 3.4,
               1.8, 0, 6.283);
        cx.fill();
      }
    }
    // oil particles in rolls
    if (showOil){
      const topFrac = 0.0;
      for (const p of oilP){
        if (play){
          const v = field(p.x, p.y);
          p.x += v.vx * dt / bw * 3.2;
          p.y += v.vy * dt / (by + bh - oilTop) * 3.2;
          if (p.x < 0.005) p.x = 0.005; if (p.x > 0.995) p.x = 0.995;
          if (p.y < 0.005) p.y = 0.005; if (p.y > 0.995) p.y = 0.995;
        }
        const px = bx + p.x * bw,
              py = oilTop + p.y * (by + bh - oilTop);
        cx.fillStyle = 'rgba(251,191,36,.7)';
        cx.beginPath(); cx.arc(px, py, 1.7, 0, 6.283); cx.fill();
      }
    }
    // chiller glyph
    cx.fillStyle = 'rgba(110,119,240,.18)';
    cx.strokeStyle = 'rgba(139,156,249,.8)'; cx.lineWidth = 1.4;
    rrect(R.chx, R.chy, R.chw, R.chh, 10, true, true);
    cx.fillStyle = '#C7D2FE'; cx.font = '600 10px Inter,sans-serif';
    cx.fillText('CHILLER', R.chx + 8, R.chy + 16);
    cx.fillStyle = '#94A3B8'; cx.font = '9.5px Inter,sans-serif';
    cx.fillText(S.chil_duty_kw.toFixed(2) + ' kW', R.chx + 8, R.chy + 32);
    cx.fillText('COP ' + S.chil_cop.toFixed(1), R.chx + 8, R.chy + 45);
    cx.fillText(S.chil_el_w.toFixed(0) + ' W el', R.chx + 8, R.chy + 58);
    const fan = tsec * 3 * spd;
    cx.strokeStyle = '#8B9CF9';
    for (let a = 0; a < 3; a++){
      cx.beginPath();
      cx.arc(R.chx + R.chw - 16, R.chy + R.chh - 18, 8,
             fan + a * 2.09, fan + a * 2.09 + 1.2);
      cx.stroke();
    }
    // weakest-link pulsing ring
    const pulse = 0.5 + 0.5 * Math.sin(tsec * 4);
    cx.strokeStyle = 'rgba(248,113,113,' + (0.35 + 0.5 * pulse) + ')';
    cx.lineWidth = 2.5;
    if (S.weak_region === 'tubes'){
      const tx = bx + tp * (Math.floor(nt / 2) + 1);
      cx.beginPath();
      cx.arc(tx, tubeY, 17 + 3 * pulse, 0, 6.283); cx.stroke();
      cx.fillStyle = 'rgba(248,113,113,.95)';
      cx.font = '600 11px Inter,sans-serif';
      cx.fillText('weakest: ' + S.weak, tx + 24, tubeY + 4);
    } else {
      const i = Math.floor(n / 2);
      const x = bx + pitch * (0.4 + i) + (pitch - cw) / 2;
      rrect(x - 3 - 2 * pulse, cTop - 3 - 2 * pulse,
            cw + 6 + 4 * pulse, (cBot - cTop) + 6 + 4 * pulse, 7,
            false, true);
      cx.fillStyle = 'rgba(248,113,113,.95)';
      cx.font = '600 11px Inter,sans-serif';
      cx.fillText('weakest: ' + S.weak, x + cw + 10, cTop + 14);
    }
    // hover highlight
    if (hover === 'cells' || sel === 'cells'){
      cx.strokeStyle = 'rgba(255,255,255,.55)'; cx.lineWidth = 1.5;
      cx.strokeRect(bx + pitch * 0.35, cTop - 3, bw - pitch * 0.7,
                    cBot - cTop + 6);
    }
    if (hover === 'chiller' || sel === 'chiller'){
      cx.strokeStyle = 'rgba(255,255,255,.6)'; cx.lineWidth = 1.6;
      rrect(R.chx - 2, R.chy - 2, R.chw + 4, R.chh + 4, 11, false, true);
    }
    // temperature legend
    const lx = bx + bw / 2 - 70, ly = by + bh + 14;
    const [lo, hi] = tRange();
    for (let i = 0; i < 140; i++){
      cx.fillStyle = tcol(lo + (hi - lo) * i / 140, 1);
      cx.fillRect(lx + i, ly, 1, 8);
    }
    cx.fillStyle = '#94A3B8'; cx.font = '9.5px Inter,sans-serif';
    cx.fillText(lo.toFixed(0) + '°C', lx - 26, ly + 8);
    cx.fillText(hi.toFixed(0) + '°C', lx + 146, ly + 8);
  } catch (e){ fail(e.message); return; }
  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);
} catch (e){ fail(e.message); }
})();
</script>
"""
    return (tpl.replace("__STATE__", s).replace("__STATIONS__", st_js)
               .replace("__BQ__", b_q).replace("__BT__", b_t)
               .replace("__BW__", b_w).replace("__BU__", b_u))
