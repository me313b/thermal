"""Cockpit v2: the whole design space on the rails.

New in v2: build_geometry, enclosure_calc, build_masses and busbar_props
are ported to the browser too, so pack architecture (format, Ns, Np,
pitch, arrangement, tube plane, tube spec) reshapes the box, tube length,
areas, oil volume and masses honestly - the v9.1 fixed-geometry
simplification is gone. Instruments: a live temperature ladder, a
flight-data-recorder strip of the governing temperature, SNAP reference
deltas on every readout, MAX-C autothrottle (in-browser bisection),
turbulence trim and dry-cooler set autopilots, and scenario presets.
The surrogate audit chip stays: solved at the design point on load and
compared against the injected full-Python result.
"""
import json


def cockpit_html(payload: dict) -> str:
    P = json.dumps(payload)
    tpl = r"""
<div id="cp-root" style="font-family:Inter,-apple-system,'Segoe UI',sans-serif">
<style>
 #cp-root{position:relative;background:#0B1220;border-radius:18px;
  overflow:hidden;color:#CBD5E1;box-shadow:0 1px 3px rgba(16,24,40,.15)}
 #cp-grid{display:grid;grid-template-columns:192px 1fr 192px;gap:0}
 .rail{padding:8px 8px 12px;background:rgba(15,23,42,.55);
  max-height:560px;overflow-y:auto}
 .rail::-webkit-scrollbar{width:6px}
 .rail::-webkit-scrollbar-thumb{background:rgba(255,255,255,.15);
  border-radius:3px}
 details{margin:0 0 7px;border:1px solid rgba(255,255,255,.09);
  border-radius:10px;background:rgba(255,255,255,.03)}
 summary{cursor:pointer;font-size:10px;letter-spacing:.12em;
  color:#8B9CF9;padding:6px 9px;user-select:none;list-style:none}
 summary::before{content:'\25B8  ';font-size:8px}
 details[open] summary::before{content:'\25BE  '}
 .ctl{margin:2px 7px 9px;padding:0}
 .ctl label{display:flex;justify-content:space-between;font-size:10px;
  color:#94A3B8;margin-bottom:2px}
 .ctl label b{color:#fff;font-size:10.5px;text-align:right}
 .ctl input[type=range]{width:100%;accent-color:#8B9CF9;height:14px}
 .seg{display:flex;gap:3px;flex-wrap:wrap}
 .seg button{flex:1;border:1px solid rgba(255,255,255,.12);
  background:transparent;color:#94A3B8;font-size:9px;padding:3px 2px;
  border-radius:6px;cursor:pointer;min-width:40px}
 .seg button.on{background:rgba(139,156,249,.25);color:#fff;
  border-color:#8B9CF9}
 select{width:100%;background:#101A2E;color:#E2E8F0;
  border:1px solid rgba(255,255,255,.12);border-radius:7px;padding:3px;
  font-size:10.5px}
 .sw{display:flex;align-items:center;justify-content:space-between;
  font-size:10.5px;color:#CBD5E1;margin:2px 7px 9px}
 #cp-cv{width:100%;height:430px;display:block;cursor:crosshair}
 #cp-inst{width:100%;height:96px;display:block;
  border-top:1px solid rgba(255,255,255,.08)}
 #cp-pfd{display:flex;gap:7px;flex-wrap:wrap;align-items:stretch;
  padding:9px 12px;background:rgba(15,23,42,.7);
  border-top:1px solid rgba(255,255,255,.1)}
 .pfd{background:rgba(255,255,255,.05);border:1px solid
  rgba(255,255,255,.10);border-radius:11px;padding:6px 11px;min-width:96px}
 .pfd span{display:block;font-size:9px;letter-spacing:.08em;
  color:#8B9CF9}
 .pfd b{font-size:15px;color:#fff}
 .pfd small{color:#94A3B8;font-size:9.5px;display:block}
 .pfd .d{font-size:10px}
 .d.up{color:#F87171}.d.dn{color:#6EE7B7}
 .bar{height:5px;border-radius:4px;background:rgba(255,255,255,.12);
  margin-top:4px;overflow:hidden}
 .bar i{display:block;height:100%;background:#34D399}
 .ann{display:flex;gap:5px;margin-left:auto;align-items:center;
  flex-wrap:wrap}
 .lamp{font-size:9px;font-weight:700;letter-spacing:.05em;
  padding:4px 8px;border-radius:7px;border:1px solid
  rgba(255,255,255,.12);color:#475569;background:rgba(255,255,255,.03)}
 .lamp.red{color:#fff;background:#B91C1C;border-color:#F87171;
  animation:bl .8s infinite alternate}
 .lamp.amb{color:#0B1220;background:#F59E0B;border-color:#FCD34D}
 .lamp.grn{color:#0B1220;background:#34D399;border-color:#6EE7B7}
 @keyframes bl{from{opacity:1}to{opacity:.55}}
 #cp-top{display:flex;align-items:center;gap:8px;padding:8px 12px;
  border-bottom:1px solid rgba(255,255,255,.1);flex-wrap:wrap}
 #cp-top h3{margin:0;font-size:13px;color:#fff;letter-spacing:.05em}
 #cp-chk{font-size:10px;color:#94A3B8}
 #cp-chk b{color:#34D399}
 .tbtn{background:rgba(139,156,249,.16);color:#E2E8F0;
  border:1px solid rgba(139,156,249,.6);border-radius:8px;
  padding:4px 10px;font-size:10.5px;cursor:pointer}
 .tbtn:hover{background:rgba(139,156,249,.32)}
 #cp-copy{margin-left:auto}
 .preset{background:rgba(255,255,255,.05);border:1px solid
  rgba(255,255,255,.14);color:#94A3B8;font-size:9.5px;
  border-radius:7px;padding:3px 8px;cursor:pointer}
 .preset:hover{color:#fff;border-color:#8B9CF9}
 #cp-tip{position:absolute;display:none;background:rgba(15,23,42,.94);
  border:1px solid rgba(255,255,255,.2);border-radius:8px;color:#E2E8F0;
  font-size:11px;padding:5px 9px;pointer-events:none;z-index:6;
  max-width:240px}
 #cp-insp{position:absolute;width:300px;max-height:300px;overflow-y:auto;
  background:rgba(15,23,42,.94);border:1px solid rgba(255,255,255,.2);
  border-radius:14px;color:#E2E8F0;font-size:12px;padding:12px 14px;
  display:none;z-index:5;right:8px;top:8px}
 #cp-insp h4{margin:0 0 6px;font-size:13px;color:#fff}
 #cp-insp .x{position:absolute;right:9px;top:7px;cursor:pointer;
  color:#94A3B8}
 #cp-insp p{margin:5px 0;line-height:1.45}
 #cp-insp b{color:#fff}
 #cp-status{font-size:10px;color:#94A3B8}
</style>
<div id="cp-top"><h3>PACK COCKPIT</h3>
 <span id="cp-chk">surrogate check: ...</span>
 <button class="tbtn" id="cp-snap">SNAP</button>
 <button class="tbtn" id="cp-maxc">A/T&nbsp;MAX-C</button>
 <button class="tbtn" id="cp-turb">TRIM&nbsp;TURB</button>
 <button class="tbtn" id="cp-dry">DRY&nbsp;SET</button>
 <span style="width:6px"></span>
 <button class="preset" data-p="base">Baseline</button>
 <button class="preset" data-p="serp">Serpentine</button>
 <button class="preset" data-p="eco">Economy 30&deg;</button>
 <button class="preset" data-p="c4">4C attempt</button>
 <button class="preset" data-p="ext">Ext pump</button>
 <button class="preset" data-p="reset"
  style="border-color:#F59E0B;color:#FCD34D">&#8635; Reset</button>
 <span id="cp-status">live</span>
 <button class="tbtn" id="cp-copy">&#10697; Copy for Design</button></div>
<div id="cp-grid">
 <div class="rail" id="railL"></div>
 <div style="position:relative"><canvas id="cp-cv"></canvas>
   <canvas id="cp-inst"></canvas>
   <div id="cp-tip"></div>
   <div id="cp-insp"><span class="x" id="cp-x">&times;</span>
     <div id="cp-insp-b"></div></div></div>
 <div class="rail" id="railR"></div>
</div>
<div id="cp-pfd"></div>
</div>
<script>
(function(){
"use strict";
const P = __PAYLOAD__;
const status=document.getElementById('cp-status');
function fail(m){status.textContent='error: '+m;status.style.color='#FCA5A5';
 try{const c=document.getElementById('cp-cv').getContext('2d');
  c.setTransform(1,0,0,1,0,0);c.fillStyle='#FCA5A5';
  c.font='13px Inter';c.fillText('Cockpit error: '+m,14,26);}catch(e){}}
window.addEventListener('error',e=>fail(e.message));
try{
const D=P.design,B=P.base,K=P.consts,FMT=P.formats,
 KT=K.KT,RT=K.RT,WATERS=P.waters;
const FL={};P.fluids.forEach(f=>FL[f.name]=f);
// ------------- controls (the whole design space) -------------
const C={c:D.C1,tamb:D.T_amb,fluid:D.coolant,
 tshape:D.tshape||'round',tw:(D.tw||0.012)*1000,th:(D.th||0.008)*1000,
 pipeD:(D.pipe_id||0.019)*1000,pipeL:D.pipe_len||2.5,
 fmt:D.fmt,ns:D.Ns,np:D.Np,cap:D.cap_Ah,rdc:D.r_dc,kdcir:D.k_dcir,
 pitch:D.pitch*1000,arr:D.arrangement,tplane:D.tube_plane,
 flow:D.flow_lpm,twin:D.T_water_in,nt:D.n_tubes,loop:D.loop_fluid,
 tod:D.tube_od*1000,twall:D.tube_wall*1000,tmat:D.tube_mat,
 fins:D.fins_on,circ:D.circ0,u:D.u0,plt:D.plate_t*1000,
 plc:D.plate_contact,tlim:D.T_limit,limc:D.limit_core,hext:D.h_ext};
// ------------- physics port (v9.1) + geometry/mass ports (v9.2) -------
const g0=9.81;
function nuOfT(f,T){const Tk=Math.max(T,-30)+273.15;
 return f.nu25*Math.exp(f.B*(1/Tk-1/298.15));}
function props(f,T){const nu=nuOfT(f,T),al=f.k/(f.rho*f.cp);
 return{k:f.k,rho:f.rho,cp:f.cp,nu,al,Pr:nu/al,beta:f.beta};}
function Ra_(p,dT,L){dT=Math.max(Math.abs(dT),0.05);
 return g0*p.beta*dT*L*L*L/(p.nu*p.al);}
function nuVert(Ra,Pr){const f=Math.pow(1+Math.pow(0.492/Pr,9/16),8/27);
 return Math.pow(0.825+0.387*Math.pow(Ra,1/6)/f,2);}
function nuHor(Ra,Pr){const f=Math.pow(1+Math.pow(0.559/Pr,9/16),8/27);
 return Math.pow(0.60+0.387*Math.pow(Ra,1/6)/f,2);}
function nuCB(Re,Pr){if(Re<1e-6)return 0;
 const a=0.62*Math.sqrt(Re)*Math.cbrt(Pr),
 b=Math.pow(1+Math.pow(0.4/Pr,2/3),0.25),
 c=Math.pow(1+Math.pow(Re/282000,5/8),4/5);return 0.3+a/b*c;}
const blend=(a,b)=>Math.cbrt(a*a*a+b*b*b);
function gapf(mm){return mm>=6?1:Math.max(0.35,
 Math.pow(Math.max(mm,0.3)/6,0.6));}
// exact build_geometry port
const _RA=[0.125,0.25,0.3333,0.5,1.0],
 _RF=[82.34,72.93,68.36,62.19,56.91],
 _RN=[5.60,4.44,3.96,3.39,2.98];
function _lerp(x,xs,ys){if(x<=xs[0])return ys[0];
 for(let i=1;i<xs.length;i++)if(x<=xs[i]){
  const w=(x-xs[i-1])/(xs[i]-xs[i-1]);
  return ys[i-1]+w*(ys[i]-ys[i-1]);}
 return ys[ys.length-1];}
function tubeSec(){const t=C.twall/1000;
 if(C.tshape==='square'){const ao=C.tw/1000,ai=Math.max(ao-2*t,1e-4);
  return{shape:'square',Ain:ai*ai,Pin:4*ai,Pout:4*ao,Dh:ai,
   metal:ao*ao-ai*ai,AoutCS:ao*ao,contact:ao,fRe:56.91,lamNu:2.98};}
 if(C.tshape==='rect'){const wo=C.tw/1000,ho=C.th/1000,
  wi=Math.max(wo-2*t,1e-4),hi=Math.max(ho-2*t,1e-4),
  A=wi*hi,Pi=2*(wi+hi),asp=Math.min(wi,hi)/Math.max(wi,hi);
  return{shape:'rect',Ain:A,Pin:Pi,Pout:2*(wo+ho),Dh:4*A/Pi,
   metal:wo*ho-wi*hi,AoutCS:wo*ho,contact:ho,
   fRe:_lerp(asp,_RA,_RF),lamNu:_lerp(asp,_RA,_RN)};}
 const od=C.tod/1000,ai=Math.max(od-2*t,1e-3);
 return{shape:'round',Ain:Math.PI*ai*ai/4,Pin:Math.PI*ai,
  Pout:Math.PI*od,Dh:ai,metal:Math.PI/4*(od*od-ai*ai),
  AoutCS:Math.PI*od*od/4,contact:0,fRe:64.0,lamNu:3.66};}
function geometry(){
 const fm=FMT[C.fmt],Dm=fm.d,H=fm.h,p=C.pitch/1000;
 const N=Math.round(C.ns)*Math.round(C.np);
 const gap=(p-Dm)*1000;
 const ncol=Math.ceil(Math.sqrt(N)),nrow=Math.ceil(N/ncol);
 const rp=p*(C.arr==='Hexagonal'?Math.sqrt(3)/2:1);
 const e=D.edge_margin;
 const Lx=ncol*p+2*e,Ly=(nrow-1)*rp+p+2*e;
 const Lz=D.bottom_gap+H+D.tube_zone+D.gas_gap;
 const fillh=Lz-D.gas_gap;
 const ts=tubeSec(),od=C.tod/1000,di=ts.Dh;
 const Lt=Math.max(Lx-2*D.manifold_margin,0.05)*D.passes;
 const nt=Math.round(C.nt);
 const fe=D.end_fraction;
 const Ac=N*(Math.PI*Dm*H+fe*2*Math.PI*Dm*Dm/4);
 const Vbox=Lx*Ly*fillh,Vc=N*Math.PI*Dm*Dm/4*H,
  Vt=nt*Lt*ts.AoutCS;
 const Aext=2*(Lx*Ly+Lx*Lz+Ly*Lz);
 const fpc=Math.max(p*rp-Math.PI*Dm*Dm/4,1e-6),
  blk=1-D.holder_block;
 const Aflow=N*fpc*blk,
  Dh=4*fpc/(Math.PI*Dm)*Math.sqrt(Math.max(blk,0.05));
 const Hloop=C.tplane==='Top of pack'?H/2+D.tube_zone/2:
  C.tplane==='Interstitial (between rows)'?0.008:H/2;
 return{N,gap,ncol,nrow,rp,Lx,Ly,Lz,fillh,od,di,ts,Lt,nt,Ac,Vbox,Vc,Vt,
  Aext,Aflow,Dh,Hloop,Dm,H,pm:p,platelen:Math.max(Lx-2*
  D.manifold_margin,0.1)};}
function finPack(g,h){const r1=g.od/2,r2=r1+D.fin_h,r2c=r2+D.fin_t/2;
 const n=1/D.fin_p,Ab=Math.PI*g.od*Math.max(0,1-D.fin_t/D.fin_p);
 const Af=n*2*Math.PI*(r2c*r2c-r1*r1);
 const m=Math.sqrt(2*Math.max(h,1)/(D.k_fin*D.fin_t));
 const Lc=D.fin_h+D.fin_t/2,phi=1+0.35*Math.log(r2c/r1),x=m*Lc*phi;
 const eta=x>1e-9?Math.tanh(x)/x:1;
 const vol=n*(Math.PI*(r2*r2-r1*r1)*D.fin_t);
 return{Aeff:Ab+eta*Af,vol};}
function plateFin(g,h,on,t,contact){if(!on)return{A:0,eta:0,m:0};
 const srv=Math.min(1,g.nt/Math.max(g.nrow-1,1));contact*=srv;
 const k=205,L=g.H,m=Math.sqrt(2*Math.max(h,5)/(k*t)),
 eta=Math.tanh(m*L)/Math.max(m*L,1e-9),npl=Math.max(g.nrow-1,1);
 return{A:npl*2*L*g.platelen*eta*contact,eta,
  m:npl*L*g.platelen*t*2700};}
function busR(g){const I=D.C1*C.cap*Math.round(C.np);
 const A=Math.max(I/D.bus_J,10);
 const L=Math.round(C.ns)*g.pm*1.15;
 return{R:1.7e-8*L/(A*1e-6),m:A*1e-6*L*8960*1.25,I0:I};}
function enclosure(g){
 const p=Math.max(D.p_des_bar*1e5,900*g0*g.fillh);
 const b=Math.min(g.Lx,g.Ly);
 const t=Math.max(b*Math.sqrt(0.31*p/(D.sigma_MPa*1e6))*D.stiff,0.0015);
 return{t,m:g.Aext*t*2700*1.18};}
function massesCalc(g,finVol,plm){
 const f=FL[C.fluid],fm=FMT[C.fmt];
 const Voil=Math.max(g.Vbox-g.Vc-g.Vt-finVol,1e-4);
 const moil=Voil*f.rho,mcell=g.N*fm.m;
 const mt=RT[C.tmat]*g.nt*g.Lt*g.ts.metal;
 const mf=finVol*(D.fin_mat==='Aluminium'?2700:8940);
 const enc=enclosure(g);
 const mst=D.struct_mass>0?D.struct_mass:enc.m;
 const mh=g.N*D.m_holder_g/1000,bb=busR(g);
 const mp=mcell+moil+mt+mf+mst+mh+bb.m+plm;
 const E=Math.round(C.ns)*Math.round(C.np)*D.v_nom*C.cap/1000;
 return{Voil,moil,mcell,mt,mf,mst,mp,E,whkg:E*1000/mp};}
function thermo(f,g,Q,Toil){const p=props(f,Toil);
 const A=g.Aflow,Dh=g.Dh,Ll=2.2*g.fillh,H=g.Hloop,Kl=K.K_loop;
 const resid=u=>p.beta*g0*H*Q/(u*A*p.cp)
  -(32*p.rho*p.nu*Ll*u/(Dh*Dh)+Kl*0.5*p.rho*u*u);
 let lo=1e-6,hi=0.08;
 for(let i=0;i<50;i++){const m=0.5*(lo+hi);
  if(resid(m)>0)lo=m;else hi=m;}
 return{u:lo,dT:Math.min(Q/(p.rho*lo*A*p.cp),60)};}
function rOfT(T){return C.rdc*Math.exp(-C.kdcir*(T-25));}
function hWater(md,di,L,lamNu,Pin){const W_=WATERS[C.loop];
 lamNu=lamNu||3.66;Pin=Pin||Math.PI*di;
 const{mu,k,cp}=W_;const Pr=mu*cp/k;
 const Re=md>0?4*md/(mu*Pin):0;
 const lam=R=>{const gz=(di/L)*R*Pr;
  return lamNu+0.0668*gz/(1+0.04*Math.pow(gz,2/3));};
 const tur=R=>{const f=Math.pow(0.790*Math.log(R)-1.64,-2);
  return (f/8)*(R-1000)*Pr/(1+12.7*Math.sqrt(f/8)*
   (Math.pow(Pr,2/3)-1));};
 let Nu,reg;
 if(Re<=0){Nu=3.66;reg='no flow';}
 else if(Re<2300){Nu=lam(Re);reg='laminar';}
 else if(Re<3000){const w=(Re-2300)/700;
  Nu=(1-w)*lam(2300)+w*tur(3000);reg='transitional';}
 else{Nu=tur(Re);reg='turbulent';}
 return{h:Nu*k/di,Re,regime:reg};}
function waterPump(g){const W_=WATERS[C.loop];
 const md=C.flow/60*W_.rho/1000,mdt=md/Math.max(g.nt,1),
 Ai=g.ts.Ain,v=mdt/(W_.rho*Ai),
 Re=W_.rho*v*g.di/W_.mu,
 f=Re<2300?g.ts.fRe/Math.max(Re,1):0.316*Math.pow(Re,-0.25),
 dp=(f*g.Lt/g.di+6)*0.5*W_.rho*v*v;
 return dp*(md/W_.rho)/0.35;}
function stirP(f,g,u){if(u<=1e-6)return 0;const p=props(f,35);
 const dp=32*p.rho*p.nu*(2.2*g.fillh)*u/(g.Dh*g.Dh)
  +K.K_loop*0.5*p.rho*u*u;
 return dp*(u*g.Aflow)/0.30;}
function serpP(f,g,u,plt){if(u<=1e-6)return 0;const p=props(f,35);
 const s=Math.max((g.pm-g.Dm-plt)/2,5e-4),
 nch=Math.max(g.nrow-1,1)*2,
 dp=12*p.rho*p.nu*g.platelen*u/(s*s)+3*0.5*p.rho*u*u;
 return dp*(u*s*g.H*nch)/0.35;}
function chiller(Qw,Tin,Tamb){const Tc=Tin-5+273.15,
 Th=Tamb+10+273.15,lift=Math.max(Th-Tc,3),
 COP=Math.max(0.45*Tc/lift,0.4);
 return{COP,Pel:Qw/COP,lift};}
function solve(cOverride){
 const cc=cOverride===undefined?C.c:cOverride;
 const f=FL[C.fluid],g=geometry(),W_=WATERS[C.loop];
 const md=C.flow/60*W_.rho/1000,mdt=md/Math.max(g.nt,1);
 const wat=hWater(mdt,g.di,g.Lt,g.ts.lamNu,g.ts.Pin);
 const Rin=1/Math.max(wat.h*g.ts.Pin*g.Lt*g.nt,1e-9);
 const Rw=g.ts.shape==='round'
  ?Math.log(g.od/g.di)/(2*Math.PI*KT[C.tmat]*g.Lt*g.nt)
  :(C.twall/1000)/(KT[C.tmat]*0.5*(g.ts.Pin+g.ts.Pout)*g.Lt*g.nt);
 const Ratm=1/Math.max(C.hext*g.Aext,1e-9);
 const ext=C.circ==='extpump',
 serp=C.circ==='serpentine'||ext,stir=C.circ==='stirred';
 const uc=(serp||stir)?C.u:0;
 const bb=busR(g);
 let Til=C.twin+8,Tb=Til+6,Twl=C.twin+2,Q=1000,uts=0,dTl=0;
 let Rb=0,Rot=0,Ao=0,pl={A:0,eta:0,m:0},hc=0,ht=0,fv=0;
 const I=cc*C.cap;
 for(let it=0;it<60;it++){
  Q=g.N*I*I*rOfT(Tb)*1e-3;
  Q+=(I*Math.round(C.np))*(I*Math.round(C.np))*bb.R;
  const ts=thermo(f,g,Math.max(Q,1),0.5*(Tb+Til));
  uts=ts.u;dTl=ts.dT;
  const ue=Math.max(uc,uts);
  const pc=props(f,0.5*(Tb+Til));
  const Rac=Ra_(pc,Tb-Til,g.H);
  let hn=nuVert(Rac,pc.Pr)*gapf(g.gap)*pc.k/g.H,hf=0;
  if(ue>1e-6)hf=nuCB(ue*g.Dm/pc.nu,pc.Pr)*pc.k/g.Dm;
  hc=blend(hn,hf)*K.cal;
  const pt=props(f,0.5*(Til+Twl));
  const Rat=Ra_(pt,Til-Twl,g.od);
  let hn2=nuHor(Rat,pt.Pr)*pt.k/g.od,hf2=0;
  if(ue>1e-6)hf2=nuCB(ue*g.od/pt.nu,pt.Pr)*pt.k/g.od;
  ht=blend(hn2,hf2)*K.cal;
  let Aeff;
  if(C.fins&&g.ts.shape==='round'){
   const fp=finPack(g,ht);Aeff=fp.Aeff*g.Lt*g.nt;
   fv=fp.vol*g.Lt*g.nt;}
  else{Aeff=g.ts.Pout*g.Lt*g.nt;fv=0;}
  pl=plateFin(g,ht,serp,C.plt/1000,C.plc);
  Ao=Aeff+pl.A;
  Rb=1/Math.max(hc*g.Ac,1e-9);
  Rot=1/Math.max(ht*Ao,1e-9);
  const Rc=Rot+Rw+Rin;
  let dTr=Q/Math.max(md*W_.cp,1e-9);
  const Ts=C.twin+0.5*Math.min(dTr,60);
  const TilN=(Q+Ts/Rc+C.tamb/Ratm)/(1/Rc+1/Ratm);
  const Qw=(TilN-Ts)/Rc;
  const TwlN=Ts+Qw*(Rin+Rw);
  const TbN=TilN+Q*Rb;
  Til+=0.6*(TilN-Til);Tb+=0.6*(TbN-Tb);Twl+=0.6*(TwlN-Twl);
 }
 const dTw=Q/Math.max(md*W_.cp,1e-9);
 const Qw=(Til-(C.twin+0.5*dTw))/(Rot+Rw+Rin);
 const Qatm=(Til-C.tamb)*C.hext*g.Aext;
 if(serp||C.tplane==='Interstitial (between rows)')dTl*=0.35;
 const spread=dTw+dTl;
 const rCore=1/(4*Math.PI*D.k_rad*g.H);
 const Tcore=Tb+(Q/g.N)*rCore;
 const Pp=waterPump(g);
 const xl=ext?extP(f,g,C.u,C.plt/1000):null;
 const Pc=ext?xl.P:serp?serpP(f,g,C.u,C.plt/1000)
  :stir?stirP(f,g,C.u):0;
 const chl=chiller(Math.max(Qw,1)+Pp,C.twin,C.tamb);
 const ms=massesCalc(g,fv,pl.m);
 const films=[['can-oil film',Rb],['oil-tube film',Rot],
  ['water film',Rin]];
 films.sort((a,b)=>b[1]-a[1]);
 return{g,Tb,Tcore,Til,Q,Qw,Qatm,dTw,dTl,spread,uts,xl,
  ue:Math.max(uc,uts),Rb,Rot,Rin,Rw,hc,ht,hw:wat.h,Re:wat.Re,
  regime:wat.regime,Ao,pl,Pp,Pc,chl,ms,
  weak:films[0][0],
  weakRegion:films[0][0]==='can-oil film'?'cells':'tubes'};}
function Tgov(r){return C.limc?r.Tcore:r.Tb;}
function maxC(){let lo=0.2,hi=6;
 for(let i=0;i<12;i++){const m=0.5*(lo+hi);
  if(Tgov(solve(m))>C.tlim)hi=m;else lo=m;}
 return lo;}
// ------------- audit -------------
let res=solve();
const dchk=res.Tb-B.T_b;
document.getElementById('cp-chk').innerHTML=
 'surrogate check at design point: <b>'+(dchk>=0?'+':'')+
 dchk.toFixed(2)+' °C</b> vs full solver'+
 (Math.abs(dchk)>0.4?' <span style="color:#F59E0B">(review)</span>':'');
// ------------- rails -------------
const railL=document.getElementById('railL'),
 railR=document.getElementById('railR');
const INP={},LBL={},FMTS={};
function group(rail,name,open){const d=document.createElement('details');
 if(open)d.open=true;
 d.innerHTML='<summary>'+name+'</summary>';
 rail.appendChild(d);return d;}
function slider(gr,key,label,min,max,step,fmt,onch){
 const d=document.createElement('div');d.className='ctl';
 d.innerHTML='<label>'+label+'<b></b></label>'+
  '<input type="range" min="'+min+'" max="'+max+'" step="'+step+
  '" value="'+C[key]+'">';
 gr.appendChild(d);
 const inp=d.querySelector('input'),vv=d.querySelector('b');
 INP[key]=inp;LBL[key]=vv;FMTS[key]=fmt;
 vv.textContent=fmt(C[key]);
 inp.addEventListener('input',()=>{C[key]=parseFloat(inp.value);
  vv.textContent=fmt(C[key]);(onch||recalc)();});
 return d;}
function segctl(gr,key,label,opts,onch){
 const d=document.createElement('div');d.className='ctl';
 d.innerHTML='<label>'+label+'</label><div class="seg"></div>';
 const seg=d.querySelector('.seg');
 opts.forEach(o=>{const b=document.createElement('button');
  b.textContent=o[1];if(C[key]===o[0])b.classList.add('on');
  b.addEventListener('click',()=>{C[key]=o[0];
   seg.querySelectorAll('button').forEach(x=>
    x.classList.remove('on'));
   b.classList.add('on');(onch||recalc)();});
  seg.appendChild(b);});
 gr.appendChild(d);return d;}
function selctl(gr,key,label,opts,onch){
 const d=document.createElement('div');d.className='ctl';
 let o='';opts.forEach(x=>{o+='<option'+(x===C[key]?' selected':'')+
  '>'+x+'</option>';});
 d.innerHTML='<label>'+label+'</label><select>'+o+'</select>';
 d.querySelector('select').addEventListener('change',e=>{
  C[key]=e.target.value;(onch||recalc)();});
 gr.appendChild(d);return d;}
function switchctl(gr,key,label,onch){
 const d=document.createElement('div');d.className='sw';
 d.innerHTML='<span>'+label+'</span><input type="checkbox" '+
  (C[key]?'checked':'')+'>';
 d.querySelector('input').addEventListener('change',e=>{
  C[key]=e.target.checked;(onch||recalc)();});
 gr.appendChild(d);return d;}
function setCtl(key,val){C[key]=val;
 if(INP[key]){INP[key].value=val;LBL[key].textContent=FMTS[key](val);}}
// LEFT: POWER / PACK / GEOMETRY
const gP=group(railL,'POWER &amp; ENVIRONMENT',true);
slider(gP,'c','C-rate',0.2,6,0.05,v=>v.toFixed(2)+' C');
slider(gP,'tamb','Ambient',0,45,1,v=>v.toFixed(0)+' °C');
const gK=group(railL,'PACK &amp; CELLS',true);
segctl(gK,'fmt','Cell format',Object.keys(FMT).map(k=>[k,k]),
 ()=>{const fm=FMT[C.fmt];setCtl('cap',fm.cap);setCtl('rdc',fm.r);
  recalc();});
slider(gK,'ns','Series (Ns)',24,240,2,v=>v.toFixed(0)+'s');
slider(gK,'np','Parallel (Np)',2,40,1,v=>v.toFixed(0)+'p');
slider(gK,'cap','Cell capacity',2,30,0.5,v=>v.toFixed(1)+' Ah');
slider(gK,'rdc','Cell DCIR (25°)',3,60,0.5,v=>v.toFixed(1)+' mΩ');
slider(gK,'kdcir','DCIR fall',0,0.03,0.001,
 v=>(v*100).toFixed(1)+' %/°C');
const gG=group(railL,'GEOMETRY');
slider(gG,'pitch','Cell pitch',18,50,0.5,v=>v.toFixed(1)+' mm');
segctl(gG,'arr','Arrangement',
 [['Hexagonal','HEX'],['Square','SQR']]);
selctl(gG,'tplane','Tube plane',['Top of pack',
 'Interstitial (between rows)','Mid-height']);
const gF=group(railL,'FLUID');
selctl(gF,'fluid','Dielectric',P.fluids.map(f=>f.name));
const gL=group(railL,'LIMITS');
slider(gL,'tlim','Cell limit',35,60,1,v=>v.toFixed(0)+' °C');
switchctl(gL,'limc','Govern on core');
slider(gL,'hext','Casing h_ext',1,15,0.5,v=>v.toFixed(1)+' W/m²·K');
// RIGHT: WATER / TUBES / CIRCULATION
const gW=group(railR,'WATER LOOP',true);
slider(gW,'flow','Flow',0.5,60,0.5,v=>v.toFixed(1)+' L/min');
slider(gW,'twin','Inlet',2,40,1,v=>v.toFixed(0)+' °C');
slider(gW,'nt','Tubes',2,48,1,v=>v.toFixed(0));
selctl(gW,'loop','Loop fluid',Object.keys(WATERS));
const gT=group(railR,'TUBES &amp; FINS',true);
segctl(gT,'tshape','Section',
 [['round','RND'],['square','SQ'],['rect','RECT']],
 ()=>{vis();recalc();});
const odRow=segctl(gT,'tod','Tube OD',[[8,'8'],[10,'10'],[12,'12']],
 ()=>recalc());
const twRow=slider(gT,'tw','Width / side',4,25,0.5,
 v=>v.toFixed(1)+' mm');
const thRow=slider(gT,'th','Height',4,25,0.5,
 v=>v.toFixed(1)+' mm');
segctl(gT,'twall','Wall',[[0.5,'0.5'],[0.8,'0.8'],[1,'1.0']]);
selctl(gT,'tmat','Material',Object.keys(KT));
const fRow=switchctl(gT,'fins','Annular fins');
const gC=group(railR,'CIRCULATION',true);
segctl(gC,'circ','Mode',
 [['thermosiphon','THERM'],['stirred','STIR'],
  ['serpentine','SERP'],['extpump','EXT']],
 ()=>{vis();recalc();});
const uRow=slider(gC,'u','Velocity',0.005,0.15,0.005,
 v=>(v*1000).toFixed(0)+' mm/s');
const pRow=slider(gC,'plt','Plate thickness',1,2,0.5,
 v=>v.toFixed(1)+' mm');
const cRow=slider(gC,'plc','Plate contact',0.4,1,0.05,
 v=>v.toFixed(2));
const dRow=slider(gC,'pipeD','Ext pipe bore',6,50,1,
 v=>v.toFixed(0)+' mm');
const lRow=slider(gC,'pipeL','Ext pipe length',0.5,8,0.25,
 v=>v.toFixed(2)+' m');
function vis(){const ext=C.circ==='extpump',
  pl=C.circ==='serpentine'||ext,rnd=C.tshape==='round';
 uRow.style.display=C.circ==='thermosiphon'?'none':'block';
 pRow.style.display=pl?'block':'none';
 cRow.style.display=pl?'block':'none';
 dRow.style.display=ext?'block':'none';
 lRow.style.display=ext?'block':'none';
 odRow.style.display=rnd?'block':'none';
 twRow.style.display=rnd?'none':'block';
 thRow.style.display=C.tshape==='rect'?'block':'none';
 fRow.style.display=rnd?'block':'none';
 if(!rnd)C.fins=false;}
vis();
// ------------- PFD -------------
const pfd=document.getElementById('cp-pfd');
function tile(id,name,sub){return '<div class="pfd"><span>'+name+
 '</span><b id="p_'+id+'"></b><small id="s_'+id+'">'+(sub||'')+
 '</small></div>';}
pfd.innerHTML=tile('t','CAN / CORE')+
 '<div class="pfd"><span>MARGIN</span><b id="p_m"></b>'+
 '<div class="bar"><i id="p_tb"></i></div></div>'+
 tile('mc','MAX CONT C','A/T target')+
 tile('s','SPREAD','criterion 5 °C')+
 tile('q','HEAT')+tile('c2','CHILLER')+tile('p','PARASITICS')+
 tile('e','ENERGY')+tile('k','PACK MASS')+tile('o','OIL')+
 '<div class="ann">'+
 '<span class="lamp" id="l_lim">LIMIT</span>'+
 '<span class="lamp" id="l_lam">LAMINAR</span>'+
 '<span class="lamp" id="l_gap">GAP&lt;6</span>'+
 '<span class="lamp" id="l_spr">SPREAD&gt;5</span>'+
 '<span class="lamp" id="l_dry">DRY-OK</span></div>';
const $=id=>document.getElementById(id);
let snap=null,hist=[],mcVal=null,mcTimer=null;
function dstr(cur,ref,dec){if(!snap)return'';
 const d=cur-ref;if(Math.abs(d)<Math.pow(10,-dec)/2)return'';
 return ' <span class="d '+(d>0?'up':'dn')+'">('+(d>0?'+':'')+
  d.toFixed(dec)+')</span>';}
function recalc(){res=solve();
 const Tg=Tgov(res),m=C.tlim-Tg;
 $('p_t').innerHTML=res.Tb.toFixed(1)+' / '+res.Tcore.toFixed(1)+
  ' °C'+dstr(Tg,snap?Tgov(snap):0,1);
 $('s_t').textContent='governing: '+(C.limc?'core':'can');
 $('p_m').textContent=(m>=0?'+':'')+m.toFixed(1)+' °C';
 $('p_m').style.color=m<0?'#F87171':m<3?'#FCD34D':'#6EE7B7';
 const fr=Math.min(Math.max((Tg-15)/(C.tlim-15),0),1.15);
 $('p_tb').style.width=(fr*100).toFixed(0)+'%';
 $('p_tb').style.background=m<0?'#EF4444':m<3?'#F59E0B':'#34D399';
 $('p_mc').textContent=mcVal===null?'...':mcVal.toFixed(2)+' C';
 $('p_s').innerHTML=res.spread.toFixed(1)+' °C'+
  dstr(res.spread,snap?snap.spread:0,1);
 $('p_q').innerHTML=(res.Q/1000).toFixed(2)+' kW'+
  dstr(res.Q/1000,snap?snap.Q/1000:0,2);
 $('s_q').textContent='casing '+res.Qatm.toFixed(0)+' W free';
 $('p_c2').textContent=(res.chl.Pel/1000).toFixed(2)+' kW el';
 $('s_c2').textContent='duty '+((res.Qw+res.Pp)/1000).toFixed(2)+
  ' kW · COP '+res.chl.COP.toFixed(1);
 $('p_p').textContent=(res.Pp+res.Pc).toFixed(1)+' W';
 $('s_p').textContent=res.xl
  ?('oil '+(res.xl.dp/1000).toFixed(1)+' kPa · pipe '
    +res.xl.v.toFixed(1)+' m/s · '+res.xl.lpm.toFixed(0)+' L/min')
  :('pump '+res.Pp.toFixed(1)+
    (res.Pc>0?' + circ '+res.Pc.toFixed(1):''));
 $('p_e').textContent=res.ms.E.toFixed(1)+' kWh';
 $('s_e').textContent=res.g.N+' cells · '+
  (Math.round(C.ns)*D.v_nom).toFixed(0)+' V';
 $('p_k').innerHTML=res.ms.mp.toFixed(0)+' kg'+
  dstr(res.ms.mp,snap?snap.ms.mp:0,0);
 $('s_k').textContent=res.ms.whkg.toFixed(0)+' Wh/kg'+
  (res.pl.m>0?' · plates '+res.pl.m.toFixed(1)+' kg':'');
 $('p_o').textContent=(res.ms.Voil*1000).toFixed(0)+' L';
 $('s_o').textContent=res.ms.moil.toFixed(0)+' kg dielectric';
 $('l_lim').className='lamp'+(m<0?' red':'');
 $('l_lam').className='lamp'+(res.regime==='laminar'?' amb':'');
 $('l_gap').className='lamp'+(res.g.gap<6?' amb':'');
 $('l_spr').className='lamp'+(res.spread>5?' amb':'');
 $('l_dry').className='lamp'+(C.twin>=C.tamb+5?' grn':'');
 hist.push(Tg);if(hist.length>170)hist.shift();
 mcVal=null;
 if(mcTimer)clearTimeout(mcTimer);
 mcTimer=setTimeout(()=>{mcVal=maxC();
  $('p_mc').textContent=mcVal.toFixed(2)+' C';},160);
}
recalc();
// autopilots + snap + presets
$('cp-snap').addEventListener('click',()=>{snap=res;
 status.textContent='reference snapped';
 setTimeout(()=>status.textContent='live',2000);recalc();});
$('cp-maxc').addEventListener('click',()=>{const v=maxC();
 setCtl('c',Math.max(Math.floor((v-0.02)*100)/100,0.2));recalc();});
$('cp-turb').addEventListener('click',()=>{const W_=WATERS[C.loop];
 const g=res.g;
 const need=3100*W_.mu*g.ts.Pin*g.nt/4*60/W_.rho*1000;
 setCtl('flow',Math.min(Math.ceil(need*2)/2,60));recalc();});
$('cp-dry').addEventListener('click',()=>{
 setCtl('twin',Math.min(C.tamb+5,40));recalc();});
document.querySelectorAll('.preset').forEach(b=>{
 b.addEventListener('click',()=>{const p=b.dataset.p;
  if(p==='reset'){location.reload();return;}
  if(p==='base'){setCtl('c',2);C.circ='thermosiphon';
   setCtl('flow',10);setCtl('twin',20);}
  if(p==='serp'){C.circ='serpentine';setCtl('u',0.05);
   setCtl('plt',1.5);setCtl('plc',0.8);}
  if(p==='eco'){setCtl('twin',30);C.circ='thermosiphon';}
  if(p==='c4'){setCtl('c',4);C.circ='serpentine';
   setCtl('u',0.08);setCtl('flow',30);setCtl('twin',15);}
  if(p==='ext'){C.circ='extpump';setCtl('u',0.05);
   setCtl('plt',1.5);setCtl('plc',0.8);
   setCtl('pipeD',19);setCtl('pipeL',2.5);}
  document.querySelectorAll('.seg button').forEach(x=>{
   if(['THERM','STIR','SERP','EXT'].includes(x.textContent))
    x.classList.toggle('on',
     (x.textContent==='THERM'&&C.circ==='thermosiphon')||
     (x.textContent==='STIR'&&C.circ==='stirred')||
     (x.textContent==='SERP'&&C.circ==='serpentine')||
     (x.textContent==='EXT'&&C.circ==='extpump'));});
  vis();recalc();});});
// copy bridge
$('cp-copy').addEventListener('click',()=>{
 const out=JSON.stringify({c1:C.c,tamb:C.tamb,fluid:C.fluid,
  ns:Math.round(C.ns),np:Math.round(C.np),cap:C.cap,rdc:C.rdc,
  kdcir:C.kdcir,pitch:C.pitch,arr:C.arr,tplane:C.tplane,
  flow:C.flow,twin:C.twin,ntub:Math.round(C.nt),loop:C.loop,
  tod:C.tod,twall:C.twall,tmat:C.tmat,fins:C.fins,circ:C.circ,
  tshape:C.tshape,tw:C.tw,th:C.th,pipeD:C.pipeD,pipeL:C.pipeL,
  u:C.u,plt:C.plt,plc:C.plc,tlim:C.tlim,limc:C.limc,hext:C.hext});
 const done=()=>{status.textContent=
  'copied - paste in the apply box below';
  setTimeout(()=>status.textContent='live',3500);};
 if(navigator.clipboard&&navigator.clipboard.writeText)
  navigator.clipboard.writeText(out).then(done,
   ()=>prompt('Copy:',out));
 else prompt('Copy:',out);});
// ------------- canvas -------------
const cv=$('cp-cv'),cx=cv.getContext('2d');
const iv=$('cp-inst'),ix=iv.getContext('2d');
const tip=$('cp-tip'),insp=$('cp-insp'),inspB=$('cp-insp-b');
$('cp-x').addEventListener('click',()=>{sel=null;
 insp.style.display='none';});
let Wc=760,Hc=430,Wi=760,Hi=96;
function fit(){Wc=cv.clientWidth||760;Hc=cv.clientHeight||430;
 Wi=iv.clientWidth||760;Hi=iv.clientHeight||96;
 const d=window.devicePixelRatio||1;
 cv.width=Math.round(Wc*d);cv.height=Math.round(Hc*d);
 cx.setTransform(d,0,0,d,0,0);
 iv.width=Math.round(Wi*d);iv.height=Math.round(Hi*d);
 ix.setTransform(d,0,0,d,0,0);}
fit();if(window.ResizeObserver)new ResizeObserver(fit).observe(cv);
let R={},sel=null,hover=null,mx=-1,my=-1,tsec=0;
function geomDraw(){
 const mL=C.circ==='extpump'?86:20,mR=20,top=24,bot=16;
 const bx=mL,by=top,bw=Wc-mL-mR,bh=Hc-top-bot;
 const g=res.g;
 const oilTop=by+bh*(1-g.fillh/g.Lz);
 const cTop=by+bh*(1-(D.bottom_gap+g.H)/g.Lz),
  cBot=by+bh*(1-D.bottom_gap/g.Lz);
 const inter=C.tplane==='Interstitial (between rows)';
 const tubeY=inter?(cTop+cBot)/2:Math.max(oilTop+20,by+24);
 R={bx,by,bw,bh,oilTop,cTop,cBot,tubeY,
  n:Math.min(g.nrow,14),nt:Math.min(g.nt,14)};
 R.pitch=bw/(R.n+0.6);R.cw=R.pitch*(g.Dm/g.pm);
 R.tp=bw/(R.nt+1);
 R.serp=(C.circ==='serpentine'||C.circ==='extpump');
 R.unit=bw/R.n;R.cwS=R.unit*0.54;R.plZ=R.unit-R.cwS;
 R.plW=R.plZ*0.62;R.midY=(cTop+cBot)/2;
 R.tr=Math.max(3,Math.min(R.plW*0.46,9));}
function cellX(i){return R.serp?R.bx+R.unit*i+(R.unit-R.cwS-R.plZ)/2+2
 :R.bx+R.pitch*(0.4+i)+(R.pitch-R.cw)/2;}
function plateX(i){return cellX(i)+R.cwS+(R.plZ-R.plW)/2;}
function tcol(T,a){const lo=Math.min(C.twin,C.tamb),
 hi=Math.max(res.Tcore,res.Tb+1,lo+8);
 let f=Math.min(Math.max((T-lo)/(hi-lo),0),1);
 return 'rgba('+Math.round(80+f*161)+','+Math.round(150-f*68)+','+
  Math.round(240-f*158)+','+a+')';}
function rrect(x,y,w,h,r,f,s){cx.beginPath();cx.moveTo(x+r,y);
 cx.arcTo(x+w,y,x+w,y+h,r);cx.arcTo(x+w,y+h,x,y+h,r);
 cx.arcTo(x,y+h,x,y,r);cx.arcTo(x,y,x+w,y,r);cx.closePath();
 if(f)cx.fill();if(s)cx.stroke();}
function field(px,py){const u=res.ue*1000*4;
 if(R.serp){const l=Math.floor(py*R.n);
  return{vx:(l%2?-1:1)*u*2.4,vy:Math.sin(px*14+l*2)*2.5};}
 const k=C.circ==='stirred'?1:3,A=u*2.2;
 return{vx:A*Math.sin(k*Math.PI*px)*Math.cos(Math.PI*py),
  vy:-A*Math.cos(k*Math.PI*px)*Math.sin(Math.PI*py)};}
const oilP=[],watP=[];
for(let i=0;i<160;i++)oilP.push({x:Math.random(),y:Math.random()});
for(let i=0;i<70;i++)watP.push({x:Math.random(),lane:i%14});
function hitTest(x,y){
 if(R.serp){
  for(let i=0;i<R.n-1;i++){const tcx=plateX(i)+R.plW/2;
   if(Math.hypot(x-tcx,y-R.midY)<R.tr+6)return 'tubes';}
  for(let i=0;i<R.n;i++){const cl=cellX(i);
   if(x>=cl&&x<=cl+R.cwS&&y>=R.cTop&&y<=R.cBot)return 'cells';}
  if(x>R.bx&&x<R.bx+R.bw&&y>R.oilTop&&y<R.by+R.bh)return 'oil';
  if(x>R.bx&&x<R.bx+R.bw&&y>R.by&&y<R.oilTop)return 'head';
  return null;}
 for(let j=0;j<R.nt;j++){const tx=R.bx+R.tp*(j+1);
  if(Math.hypot(x-tx,y-R.tubeY)<15)return 'tubes';}
 if(y>=R.cTop&&y<=R.cBot&&x>R.bx&&x<R.bx+R.bw)return 'cells';
 if(x>R.bx&&x<R.bx+R.bw&&y>R.oilTop&&y<R.by+R.bh)return 'oil';
 if(x>R.bx&&x<R.bx+R.bw&&y>R.by&&y<R.oilTop)return 'head';
 return null;}
function stText(id){
 if(id==='cells')return{t:'Cells - '+res.g.N+' x '+C.fmt,
  h:'can '+res.Tb.toFixed(1)+' °C · core '+res.Tcore.toFixed(1)+' °C',
  b:'<p>Each cell: <b>'+(res.Q/res.g.N).toFixed(2)+' W</b> at '+
  C.c.toFixed(2)+'C. First film h = <b>'+res.hc.toFixed(0)+
  '</b> over '+res.g.Ac.toFixed(1)+' m² -> '+
  (res.Q*res.Rb).toFixed(1)+' °C. Gap '+res.g.gap.toFixed(1)+
  ' mm (penalty x'+gapf(res.g.gap).toFixed(2)+'). DCIR '+
  rOfT(res.Tb).toFixed(1)+' mΩ.</p>'};
 if(id==='oil')return{t:'Bulk oil',
  h:res.Til.toFixed(1)+' °C · '+(res.ue*1000).toFixed(1)+' mm/s',
  b:'<p>Oil <b>'+res.Til.toFixed(1)+' °C</b>, '+
  (res.ue*1000).toFixed(1)+' mm/s ('+C.circ+'). Spread <b>'+
  res.spread.toFixed(1)+' °C</b>. '+
  (res.ms.Voil*1000).toFixed(0)+' L / '+res.ms.moil.toFixed(0)+
  ' kg; casing sheds '+res.Qatm.toFixed(0)+' W.</p>'};
 if(id==='tubes')return{t:'Tubes, fins'+
  (C.circ==='serpentine'?', plates':''),
  h:'oil h '+res.ht.toFixed(0)+' · Re '+res.Re.toFixed(0)+
  ' ('+res.regime+')',
  b:'<p>'+res.g.nt+' x '+(C.tshape==='round'?C.tod+' mm':
   C.tshape==='square'?C.tw+' mm square':
   C.tw+'x'+C.th+' mm rect')+' '+C.tmat.toLowerCase()+
  ', L = '+res.g.Lt.toFixed(2)+' m each. Second film h = '+
  res.ht.toFixed(0)+', A = <b>'+res.Ao.toFixed(1)+' m²</b>'+
  (res.pl.A>0?' (plates +'+res.pl.A.toFixed(1)+' m², η '+
  res.pl.eta.toFixed(2)+')':'')+' -> '+
  (res.Q*res.Rot).toFixed(1)+' °C. Water: <b>'+res.regime+
  '</b>, Re '+res.Re.toFixed(0)+' -> '+
  (res.Q*res.Rin).toFixed(1)+' °C; rise '+res.dTw.toFixed(1)+
  ' °C.</p>'};
 return{t:'Headspace & enclosure',h:'gas blanket',
  b:'<p>'+(D.gas_gap*1000).toFixed(0)+' mm gas blanket. Box '+
  (res.g.Lx*1000).toFixed(0)+' x '+(res.g.Ly*1000).toFixed(0)+
  ' x '+(res.g.Lz*1000).toFixed(0)+' mm; structure '+
  res.ms.mst.toFixed(0)+' kg at the burst rating.</p>'};}
cv.addEventListener('mousemove',e=>{
 const r=cv.getBoundingClientRect();
 mx=e.clientX-r.left;my=e.clientY-r.top;hover=hitTest(mx,my);
 if(hover){const s=stText(hover);
  tip.innerHTML='<b>'+s.t+'</b><br>'+s.h+
   '<br><span style="color:#94A3B8">click to inspect</span>';
  tip.style.display='block';
  tip.style.left=Math.min(mx+14,Wc-250)+'px';
  tip.style.top=Math.min(my+12,Hc-70)+'px';
  cv.style.cursor='pointer';}
 else{tip.style.display='none';cv.style.cursor='crosshair';}});
cv.addEventListener('mouseleave',()=>{tip.style.display='none';
 hover=null;});
cv.addEventListener('click',()=>{if(!hover){sel=null;
  insp.style.display='none';return;}
 sel=hover;const s=stText(sel);
 inspB.innerHTML='<h4>'+s.t+'</h4>'+s.b;
 insp.style.display='block';});
// instruments strip: ladder (left) + FDR (right)
function instruments(){
 ix.fillStyle='#0B1220';ix.fillRect(0,0,Wi,Hi);
 const lw=Wi*0.55-30,lx0=16,ly0=34,lh=18;
 const segs=[['water rise',0.5*res.dTw,'#38BDF8'],
  ['water film',res.Q*res.Rin,'#0EA5E9'],
  ['wall',res.Q*res.Rw,'#94A3B8'],
  ['oil-tube film',res.Q*res.Rot,'#F59E0B'],
  ['can-oil film',res.Q*res.Rb,'#F97316'],
  ['core rise',res.Tcore-res.Tb,'#EF4444']];
 const tot=segs.reduce((a,s)=>a+s[1],0);
 ix.fillStyle='#8B9CF9';ix.font='600 9px Inter';
 ix.fillText('TEMPERATURE LADDER  '+C.twin.toFixed(0)+' °C -> '+
  res.Tcore.toFixed(1)+' °C  ('+tot.toFixed(1)+' °C total)',lx0,14);
 let x=lx0;
 segs.forEach(s=>{const w=Math.max(tot>0?s[1]/tot*lw:0,0);
  ix.fillStyle=s[2];ix.fillRect(x,ly0,Math.max(w-1.5,0.5),lh);
  if(w>34){ix.fillStyle='#0B1220';ix.font='8.5px Inter';
   ix.fillText(s[1].toFixed(1)+'°',x+3,ly0+12);}
  x+=w;});
 ix.fillStyle='#64748B';ix.font='8.5px Inter';
 let x2=lx0;segs.forEach(s=>{const w=tot>0?s[1]/tot*lw:0;
  if(w>52)ix.fillText(s[0],x2+3,ly0+lh+11);x2+=w;});
 const fx0=Wi*0.55+14,fw=Wi-fx0-16,fy0=22,fh=Hi-34;
 ix.strokeStyle='rgba(255,255,255,.12)';
 ix.strokeRect(fx0,fy0,fw,fh);
 ix.fillStyle='#8B9CF9';ix.font='600 9px Inter';
 ix.fillText('FDR - governing T',fx0,14);
 const lo=Math.min(...hist,C.twin),hi=Math.max(...hist,C.tlim+2);
 const yOf=v=>fy0+fh-(v-lo)/(hi-lo)*fh;
 ix.strokeStyle='#F87171';ix.setLineDash([4,3]);
 ix.beginPath();ix.moveTo(fx0,yOf(C.tlim));
 ix.lineTo(fx0+fw,yOf(C.tlim));ix.stroke();ix.setLineDash([]);
 ix.strokeStyle='#34D399';ix.lineWidth=1.5;ix.beginPath();
 hist.forEach((v,i)=>{const px=fx0+i/(170-1)*fw,py=yOf(v);
  i?ix.lineTo(px,py):ix.moveTo(px,py);});
 ix.stroke();ix.lineWidth=1;
 ix.fillStyle='#94A3B8';ix.font='8.5px Inter';
 ix.fillText(hi.toFixed(0),fx0+fw+2-14,fy0+8);
 ix.fillText(lo.toFixed(0),fx0+fw+2-14,fy0+fh);
}
let last=performance.now(),acc=0;
function frame(now){try{
 const dt=Math.min((now-last)/1000,0.05);last=now;tsec+=dt;acc+=dt;
 if(acc>0.5){hist.push(Tgov(res));
  if(hist.length>170)hist.shift();acc=0;}
 geomDraw();
 const{bx,by,bw,bh,oilTop,cTop,cBot,tubeY,n,nt,pitch,cw,tp}=R;
 cx.fillStyle='#0B1220';cx.fillRect(0,0,Wc,Hc);
 cx.strokeStyle='rgba(255,255,255,.45)';cx.lineWidth=2;
 cx.strokeRect(bx-6,by-6,bw+12,bh+12);
 cx.fillStyle='#0E1830';cx.fillRect(bx,by,bw,oilTop-by);
 cx.fillStyle='rgba(245,158,11,0.10)';
 cx.fillRect(bx,oilTop,bw,by+bh-oilTop);
 cx.strokeStyle='rgba(245,158,11,.6)';cx.lineWidth=1;
 cx.beginPath();cx.moveTo(bx,oilTop);cx.lineTo(bx+bw,oilTop);
 cx.stroke();
 const serp=R.serp,unit=R.unit,cwS=R.cwS,plZ=R.plZ,plW=R.plW,
  midY=R.midY,tr=R.tr;
 for(let i=0;i<n;i++){const x=cellX(i),w=serp?cwS:cw;
  const gcx=x+w/2;
  const gl=cx.createRadialGradient(gcx,midY,2,gcx,midY,w*1.7);
  const a=Math.min(0.12+res.Q/1000*0.05,0.4);
  gl.addColorStop(0,'rgba(241,82,82,'+a+')');
  gl.addColorStop(1,'rgba(241,82,82,0)');
  cx.fillStyle=gl;cx.fillRect(x-w,cTop-w,w*3,(cBot-cTop)+2*w);
  cx.fillStyle=tcol(res.Tb,0.96);
  cx.strokeStyle='rgba(255,255,255,.30)';cx.lineWidth=1;
  rrect(x,cTop,w,cBot-cTop,5,true,true);
  cx.fillStyle=tcol(res.Tcore,0.95);
  rrect(x+w*0.30,cTop+4,w*0.40,(cBot-cTop)-8,4,true,false);
  if(serp&&i<n-1){const plx=plateX(i);
   const mg=cx.createLinearGradient(plx,0,plx+plW,0);
   mg.addColorStop(0,'rgba(148,163,184,.5)');
   mg.addColorStop(.5,'rgba(216,224,235,.94)');
   mg.addColorStop(1,'rgba(148,163,184,.5)');
   cx.fillStyle=mg;cx.strokeStyle='rgba(226,232,240,.6)';cx.lineWidth=1;
   rrect(plx,cTop,plW,cBot-cTop,3,true,true);
   const tcx=plx+plW/2;
   cx.strokeStyle='rgba(120,134,158,.85)';cx.lineWidth=1.5;
   cx.beginPath();cx.moveTo(tcx,cTop+5);cx.lineTo(tcx,cBot-5);cx.stroke();
   cx.fillStyle='rgba(71,85,105,.96)';
   cx.beginPath();cx.arc(tcx,midY,tr+1.6,0,6.283);cx.fill();
   cx.fillStyle=tcol(C.twin+res.dTw*0.5,0.96);
   cx.beginPath();cx.arc(tcx,midY,tr,0,6.283);cx.fill();
   cx.fillStyle='rgba(186,230,253,.75)';
   cx.beginPath();cx.arc(tcx,midY,tr*0.42,0,6.283);cx.fill();}}
 if(!serp){for(let j=0;j<nt;j++){const tx=bx+tp*(j+1);
  cx.fillStyle='rgba(165,180,204,.4)';
  cx.beginPath();cx.arc(tx,tubeY,11,0,6.283);cx.fill();
  cx.fillStyle='#D97706';
  cx.beginPath();cx.arc(tx,tubeY,5.5,0,6.283);cx.fill();}}
 cx.font='10px Inter';
 if(serp){cx.fillStyle='#93C5FD';
  cx.fillText('tubes bonded in the plates, running INTO the page  ·  in '+
   C.twin.toFixed(0)+'°C (front) -> '+(C.twin+res.dTw).toFixed(1)+
   '°C (back)',bx+4,by+14);}
 else{cx.fillStyle='#7DD3FC';
  cx.fillText('in '+C.twin.toFixed(0)+'°C',bx+4,tubeY-16);
  cx.fillStyle=tcol(C.twin+res.dTw,1);
  cx.fillText('out '+(C.twin+res.dTw).toFixed(1)+'°C',bx+bw-64,tubeY-16);}
 for(const p of watP){p.x+=dt*Math.min(C.flow/20,2)*0.25;if(p.x>1)p.x-=1;
  const ang=p.x*6.283;let tx,ty,rr;
  if(serp){const ln=p.lane%Math.max(n-1,1);tx=plateX(ln)+plW/2;ty=midY;
   rr=tr*0.6;}
  else{tx=bx+tp*((p.lane%nt)+1);ty=tubeY;rr=3.4;}
  cx.fillStyle=tcol(C.twin+res.dTw*p.x,0.95);
  cx.beginPath();cx.arc(tx+Math.cos(ang)*rr,ty+Math.sin(ang)*rr,
   serp?1.4:1.8,0,6.283);cx.fill();}
 for(const p of oilP){const v=field(p.x,p.y);
  p.x+=v.vx*dt/bw*3.2;p.y+=v.vy*dt/(by+bh-oilTop)*3.2;
  if(p.x<0.005)p.x=0.005;if(p.x>0.995)p.x=0.995;
  if(p.y<0.005)p.y=0.005;if(p.y>0.995)p.y=0.995;
  const px=bx+p.x*bw,py=oilTop+p.y*(by+bh-oilTop);
  cx.fillStyle='rgba(251,191,36,.7)';
  cx.beginPath();cx.arc(px,py,1.7,0,6.283);cx.fill();}
 const pulse=0.5+0.5*Math.sin(tsec*4);
 cx.strokeStyle='rgba(248,113,113,'+(0.35+0.5*pulse)+')';
 cx.lineWidth=2.5;
 if(res.weakRegion==='tubes'){
  if(serp){const ln=Math.floor((n-1)/2),tcx=plateX(ln)+plW/2;
   cx.beginPath();cx.arc(tcx,midY,tr+7+3*pulse,0,6.283);cx.stroke();
   cx.fillStyle='rgba(248,113,113,.95)';cx.font='600 11px Inter';
   cx.fillText('weakest: '+res.weak,tcx+14,midY-tr-6);}
  else{const tx=bx+tp*(Math.floor(nt/2)+1);
   cx.beginPath();cx.arc(tx,tubeY,17+3*pulse,0,6.283);cx.stroke();
   cx.fillStyle='rgba(248,113,113,.95)';cx.font='600 11px Inter';
   cx.fillText('weakest: '+res.weak,tx+24,tubeY+4);}}
 else{const i=Math.floor(n/2),x=cellX(i),w=serp?cwS:cw;
  rrect(x-3-2*pulse,cTop-3-2*pulse,w+6+4*pulse,
   (cBot-cTop)+6+4*pulse,7,false,true);
  cx.fillStyle='rgba(248,113,113,.95)';cx.font='600 11px Inter';
  cx.fillText('weakest: '+res.weak,x+w+10,cTop+14);}
 if(C.circ==='extpump'){
  const yS=by+bh-14,yR=oilTop+16,xw=bx-6,xp2=34,
   pcx=(xw+xp2)/2-4,pcy=(yS+yR)/2,pr=13;
  cx.strokeStyle='rgba(245,158,11,.85)';cx.lineWidth=4;
  cx.lineCap='round';
  cx.beginPath();cx.moveTo(xw,yS);cx.lineTo(xp2,yS);
  cx.lineTo(xp2,pcy+pr);cx.stroke();
  cx.beginPath();cx.moveTo(xp2,pcy-pr);cx.lineTo(xp2,yR);
  cx.lineTo(xw,yR);cx.stroke();
  cx.fillStyle='#0B1220';cx.strokeStyle='#F59E0B';cx.lineWidth=2.5;
  cx.beginPath();cx.arc(xp2,pcy,pr,0,6.283);cx.fill();cx.stroke();
  const ra=tsec*(2+C.u*40);
  cx.strokeStyle='rgba(252,211,77,.9)';cx.lineWidth=2;
  for(let k=0;k<3;k++){const a=ra+k*2.094;
   cx.beginPath();cx.moveTo(xp2,pcy);
   cx.lineTo(xp2+Math.cos(a)*(pr-4),pcy+Math.sin(a)*(pr-4));
   cx.stroke();}
  const path=[[xw,yS],[xp2,yS],[xp2,pcy+pr],[xp2,pcy-pr],
   [xp2,yR],[xw,yR]];
  let Ltot=0;const seg=[];
  for(let k=0;k<path.length-1;k++){const L2=Math.hypot(
   path[k+1][0]-path[k][0],path[k+1][1]-path[k][1]);
   seg.push(L2);Ltot+=L2;}
  cx.fillStyle='rgba(251,191,36,.95)';
  for(let k=0;k<6;k++){
   let s=((tsec*C.u*260)+k*Ltot/6)%Ltot,j=0;
   while(s>seg[j]){s-=seg[j];j++;}
   const fx=path[j][0]+(path[j+1][0]-path[j][0])*s/seg[j],
    fy=path[j][1]+(path[j+1][1]-path[j][1])*s/seg[j];
   cx.beginPath();cx.arc(fx,fy,2.1,0,6.283);cx.fill();}
  cx.fillStyle='#FCD34D';cx.font='600 9px Inter';
  cx.fillText('PUMP',xp2-13,pcy+pr+12);
  if(res.xl){cx.fillStyle='#94A3B8';cx.font='8px Inter';
   cx.fillText((res.xl.dp/1000).toFixed(1)+' kPa',4,pcy-pr-8);
   cx.fillText(res.xl.P.toFixed(0)+' W',4,pcy-pr+2);}
  cx.fillStyle='#64748B';cx.font='8px Inter';
  cx.fillText('closed loop',2,yS+11);
  cx.fillText('no ext. HX',2,yR-8);}
 cx.fillStyle='#64748B';cx.font='9px Inter';
 cx.fillText(res.g.N+' cells · box '+(res.g.Lx*1000).toFixed(0)+
  ' x '+(res.g.Ly*1000).toFixed(0)+' x '+
  (res.g.Lz*1000).toFixed(0)+' mm',bx+4,by+bh+12);
 instruments();
}catch(e){fail(e.message);return;}
 requestAnimationFrame(frame);}
requestAnimationFrame(frame);
}catch(e){fail(e.message);}
})();
</script>
"""
    return tpl.replace("__PAYLOAD__", P)
