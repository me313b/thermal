"""pack3d.py - an interactive 3D view of the serpentine plate-channel
architecture for the cockpit. Cells stand in a grid; a metal cooling
plate sits in each gap between cell columns; water tubes are bonded along
each plate and run front-to-back (into the page); coolant flows
through-plane, entering the front face and leaving the back, warming as
it goes. Drag to rotate; the flow animates. Hand-rolled orthographic 3D
(rotation matrix + painter's algorithm), so it is fully self-contained.
"""
import json


def pack3d_html(p):
    payload = json.dumps(p)
    return """
<div id="p3root" style="font-family:Inter,system-ui,sans-serif;color:#cbd5e1">
  <canvas id="p3" width="760" height="440"
          style="width:100%;max-width:760px;border-radius:12px;
                 background:radial-gradient(1200px 500px at 50% -10%,
                 #16233f 0%,#0b1220 60%);cursor:grab;touch-action:none">
  </canvas>
  <div style="display:flex;gap:14px;align-items:center;flex-wrap:wrap;
              margin-top:8px;font-size:12.5px">
    <button id="p3spin" style="background:#1e293b;color:#e2e8f0;
      border:1px solid #334155;border-radius:8px;padding:5px 12px;
      cursor:pointer">⏸ pause spin</button>
    <label>flow <input id="p3flow" type="range" min="0" max="100"
      value="55" style="vertical-align:middle;width:130px"></label>
    <label>explode <input id="p3exp" type="range" min="0" max="100"
      value="0" style="vertical-align:middle;width:110px"></label>
    <span id="p3lab" style="color:#94a3b8"></span>
  </div>
  <div style="font-size:11.5px;color:#64748b;margin-top:4px">
    Drag to rotate. Cells (amber) sit in a grid; a metal plate (grey)
    fills each gap; water tubes (blue) run bonded along the plates,
    <b>into the page</b>; coolant flows front→back and warms as it
    collects heat. This is the through-plane path a flat 2D view cannot
    show.
  </div>
</div>
<script>
(function(){
 const P=__PAYLOAD__;
 const cv=document.getElementById('p3'), ctx=cv.getContext('2d');
 const W=cv.width, H=cv.height;
 const lab=document.getElementById('p3lab');
 let yaw=-0.62, pitch=0.42, spin=true, t0=performance.now();
 const flow=document.getElementById('p3flow'),
       expS=document.getElementById('p3exp');
 // --- geometry in mm, centred ---
 const D=P.D, Hc=P.Ht, pitch_mm=P.pitch, plt=P.plate_t, tod=P.tube_od;
 const nx=P.nx, ny=P.ny, ntz=P.tubes_z;   // tubes stacked in z per plate
 const depth=(ny-1)*pitch_mm + D;          // total through-plane depth
 const spanX=(nx-1)*pitch_mm;
 const cx0=-spanX/2, cy0=-depth/2, cz0=-Hc/2;
 // scale to fit
 const R=D/2;
 let scale=Math.min(W/(spanX+depth+2*D), H/(Hc+depth))*0.62;
 const CX=W*0.5, CY=H*0.56;

 function rot(x,y,z,ex){
   // explode: push plates+tubes apart in x a little
   const cz=Math.cos(yaw), sz=Math.sin(yaw);
   let X=x*cz - y*sz, Y=x*sz + y*cz, Z=z;
   const cp=Math.cos(pitch), sp=Math.sin(pitch);
   let Y2=Y*cp - Z*sp, Z2=Y*sp + Z*cp;
   return {sx:CX + X*scale, sy:CY - Z2*scale, d:Y2};
 }
 function tcol(f){ // 0 cool -> 1 hot
   const r=Math.round(56+(239-56)*f), g=Math.round(189+(68-189)*f),
         b=Math.round(248+(68-248)*f); return [r,g,b];
 }
 function rgba(c,a){return 'rgba('+c[0]+','+c[1]+','+c[2]+','+a+')';}

 function build(ex){
   const prim=[];
   const Tc=P.Tcell, Tin=P.Tw_in, Tout=P.Tw_out, lim=P.limit;
   const fcell=Math.max(0,Math.min(1,(Tc-P.Tw_in)/Math.max(lim-P.Tw_in,1)));
   // cells (cylinders): grid nx x ny
   for(let i=0;i<nx;i++) for(let j=0;j<ny;j++){
     const x=cx0+i*pitch_mm, y=cy0+j*(pitch_mm), 
           top=rot(x,y,cz0+Hc), bot=rot(x,y,cz0);
     prim.push({t:'cell',a:top,b:bot,r:R*scale,
                d:(top.d+bot.d)/2, f:fcell, x:x});
   }
   // plates in the gaps between columns (nx-1 of them), y-z sheets
   for(let i=0;i<nx-1;i++){
     const xg=cx0+(i+0.5)*pitch_mm + (ex? (i-(nx-2)/2)*ex*8:0);
     const c=[rot(xg,cy0,cz0),rot(xg,cy0+depth,cz0),
              rot(xg,cy0+depth,cz0+Hc),rot(xg,cy0,cz0+Hc)];
     prim.push({t:'plate',c:c,d:(c[0].d+c[2].d)/2});
     // tubes bonded along this plate, running in y (into page)
     for(let k=0;k<ntz;k++){
       const zt=cz0+Hc*(k+1)/(ntz+1);
       const fr=rot(xg,cy0,zt), bk=rot(xg,cy0+depth,zt);
       prim.push({t:'tube',a:fr,b:bk,r:(tod/2)*scale,
                  d:(fr.d+bk.d)/2, xg:xg, zt:zt});
     }
   }
   return prim.sort((A,B)=>B.d-A.d);   // far first
 }

 // moving coolant markers per tube (parametric along y)
 const parts=[];
 function seedParts(){
   parts.length=0;
   for(let i=0;i<nx-1;i++) for(let k=0;k<P.tubes_z;k++)
     for(let m=0;m<5;m++) parts.push({i:i,k:k,s:Math.random()});
 }
 seedParts();

 function draw(now){
   const ex=+expS.value/100, fv=+flow.value/100;
   scale=Math.min(W/(spanX+depth+2*D), H/(Hc+depth))*0.62*(1-0.12*ex);
   if(spin) yaw += 0.0035;
   ctx.clearRect(0,0,W,H);
   const prim=build(ex);
   for(const o of prim){
     if(o.t==='cell'){
       // body as round-capped thick line
       ctx.lineCap='round'; ctx.lineWidth=2*o.r;
       const g=ctx.createLinearGradient(o.a.sx,o.a.sy,o.b.sx,o.b.sy);
       const c=tcol(o.f);
       g.addColorStop(0,rgba(c,.95)); g.addColorStop(1,'rgba(120,53,15,.95)');
       ctx.strokeStyle=g;
       ctx.beginPath(); ctx.moveTo(o.a.sx,o.a.sy);
       ctx.lineTo(o.b.sx,o.b.sy); ctx.stroke();
       // top cap highlight
       ctx.fillStyle=rgba(c,.6);
       ctx.beginPath(); ctx.ellipse(o.a.sx,o.a.sy,o.r,o.r*Math.max(.25,
         Math.abs(Math.sin(pitch))),0,0,6.2832); ctx.fill();
       ctx.strokeStyle='rgba(251,191,36,.5)'; ctx.lineWidth=1;
       ctx.stroke();
     } else if(o.t==='plate'){
       ctx.beginPath(); ctx.moveTo(o.c[0].sx,o.c[0].sy);
       for(let q=1;q<4;q++) ctx.lineTo(o.c[q].sx,o.c[q].sy);
       ctx.closePath();
       ctx.fillStyle='rgba(148,163,184,.42)';
       ctx.strokeStyle='rgba(203,213,225,.55)'; ctx.lineWidth=1;
       ctx.fill(); ctx.stroke();
     } else { // tube
       ctx.lineCap='round'; ctx.lineWidth=2*o.r;
       ctx.strokeStyle='rgba(14,165,233,.95)';
       ctx.beginPath(); ctx.moveTo(o.a.sx,o.a.sy);
       ctx.lineTo(o.b.sx,o.b.sy); ctx.stroke();
       ctx.lineWidth=Math.max(1,o.r*0.5);
       ctx.strokeStyle='rgba(186,230,253,.5)';
       ctx.beginPath(); ctx.moveTo(o.a.sx,o.a.sy);
       ctx.lineTo(o.b.sx,o.b.sy); ctx.stroke();
     }
   }
   // coolant particles (drawn after tubes, front→back, warming)
   for(const pt of parts){
     pt.s += (0.002+0.02*fv);
     if(pt.s>1){ pt.s-=1; }
     const xg=cx0+(pt.i+0.5)*pitch_mm+(ex?(pt.i-(nx-2)/2)*ex*8:0);
     const zt=cz0+Hc*(pt.k+1)/(P.tubes_z+1);
     const y=cy0+pt.s*depth;
     const pr=rot(xg,y,zt);
     const c=tcol(0.15+0.6*pt.s);
     ctx.fillStyle=rgba(c,.95);
     ctx.beginPath(); ctx.arc(pr.sx,pr.sy,Math.max(1.6,tod*scale*0.28),
       0,6.2832); ctx.fill();
   }
   // inlet / outlet labels
   const inl=rot(cx0+ (nx-2>=0? (0.5)*pitch_mm : 0), cy0, cz0+Hc*0.5);
   const out=rot(cx0+ (nx-2>=0? (0.5)*pitch_mm : 0), cy0+depth, cz0+Hc*0.5);
   ctx.font='12px Inter,sans-serif'; ctx.textAlign='center';
   ctx.fillStyle='#7dd3fc';
   ctx.fillText('coolant in '+P.Tw_in.toFixed(0)+'°C', inl.sx, inl.sy-12);
   ctx.fillStyle='#fca5a5';
   ctx.fillText('out '+P.Tw_out.toFixed(1)+'°C', out.sx, out.sy-12);
   lab.innerHTML='cells '+P.Tcell.toFixed(1)+'°C &nbsp;·&nbsp; plate → tube '+
     '→ water &nbsp;·&nbsp; showing '+nx+'×'+ny+' of '+P.ns_np;
   requestAnimationFrame(draw);
 }

 // drag to rotate
 let drag=false, px=0, py=0;
 cv.addEventListener('pointerdown',e=>{drag=true;spin=false;
   document.getElementById('p3spin').textContent='▶ resume spin';
   px=e.clientX;py=e.clientY;cv.style.cursor='grabbing';});
 window.addEventListener('pointermove',e=>{if(!drag)return;
   yaw+=(e.clientX-px)*0.008; pitch+=(e.clientY-py)*0.006;
   pitch=Math.max(-0.2,Math.min(1.2,pitch)); px=e.clientX;py=e.clientY;});
 window.addEventListener('pointerup',()=>{drag=false;cv.style.cursor='grab';});
 document.getElementById('p3spin').addEventListener('click',function(){
   spin=!spin; this.textContent=spin?'⏸ pause spin':'▶ resume spin';});
 requestAnimationFrame(draw);
})();
</script>
""".replace("__PAYLOAD__", payload)
