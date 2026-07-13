"""Cockpit: a pilot's control unit for the pack.

Every lever live, every readout linked, at animation frame rate. The
physics runs in the browser as a line-faithful port of the app's solver
(Churchill-Chu films with the gap penalty and mixed-convection blend,
Andrade viscosity at film temperature, Schmidt fins, plate fins, the
thermosiphon bisection, Hausen/Gnielinski water side, DCIR-coupled 60-pass
fixed point, pump laws, 45%-Carnot chiller). On load the surrogate is
audited against the injected full-Python solve and the deviation is shown
on the panel. A copy button exports the current settings for the Design
tab. Errors paint themselves on the canvas.
"""
import json


def cockpit_html(payload: dict) -> str:
    P = json.dumps(payload)
    tpl = r"""
<div id="cp-root" style="font-family:Inter,-apple-system,'Segoe UI',sans-serif">
<style>
 #cp-root{position:relative;background:#0B1220;border-radius:18px;
  overflow:hidden;color:#CBD5E1;box-shadow:0 1px 3px rgba(16,24,40,.15)}
 #cp-grid{display:grid;grid-template-columns:172px 1fr 172px;gap:0}
 .rail{padding:12px 10px;background:rgba(15,23,42,.55);min-height:520px}
 .rail h5{margin:2px 4px 8px;font-size:10px;letter-spacing:.14em;
  color:#8B9CF9}
 .ctl{margin:0 2px 11px;background:rgba(255,255,255,.04);
  border:1px solid rgba(255,255,255,.09);border-radius:10px;
  padding:7px 9px}
 .ctl label{display:flex;justify-content:space-between;font-size:10.5px;
  color:#94A3B8;margin-bottom:3px}
 .ctl label b{color:#fff;font-size:11px}
 .ctl input[type=range]{width:100%;accent-color:#8B9CF9;height:16px}
 .seg{display:flex;gap:3px}
 .seg button{flex:1;border:1px solid rgba(255,255,255,.12);
  background:transparent;color:#94A3B8;font-size:9.5px;padding:4px 2px;
  border-radius:7px;cursor:pointer}
 .seg button.on{background:rgba(139,156,249,.25);color:#fff;
  border-color:#8B9CF9}
 select{width:100%;background:#101A2E;color:#E2E8F0;
  border:1px solid rgba(255,255,255,.12);border-radius:7px;padding:4px;
  font-size:11px}
 .sw{display:flex;align-items:center;justify-content:space-between;
  font-size:11px;color:#CBD5E1}
 #cp-cv{width:100%;height:520px;display:block;cursor:crosshair}
 #cp-pfd{display:flex;gap:8px;flex-wrap:wrap;align-items:stretch;
  padding:10px 12px;background:rgba(15,23,42,.7);
  border-top:1px solid rgba(255,255,255,.1)}
 .pfd{background:rgba(255,255,255,.05);border:1px solid
  rgba(255,255,255,.10);border-radius:11px;padding:7px 12px;min-width:104px}
 .pfd span{display:block;font-size:9.5px;letter-spacing:.08em;
  color:#8B9CF9}
 .pfd b{font-size:16px;color:#fff}
 .pfd small{color:#94A3B8;font-size:10px}
 .bar{height:5px;border-radius:4px;background:rgba(255,255,255,.12);
  margin-top:4px;overflow:hidden}
 .bar i{display:block;height:100%;background:#34D399}
 .ann{display:flex;gap:6px;margin-left:auto;align-items:center;
  flex-wrap:wrap}
 .lamp{font-size:9.5px;font-weight:700;letter-spacing:.06em;
  padding:4px 9px;border-radius:7px;border:1px solid
  rgba(255,255,255,.12);color:#475569;background:rgba(255,255,255,.03)}
 .lamp.red{color:#fff;background:#B91C1C;border-color:#F87171;
  animation:bl .8s infinite alternate}
 .lamp.amb{color:#0B1220;background:#F59E0B;border-color:#FCD34D}
 .lamp.grn{color:#0B1220;background:#34D399;border-color:#6EE7B7}
 @keyframes bl{from{opacity:1}to{opacity:.55}}
 #cp-top{display:flex;align-items:center;gap:10px;padding:9px 14px;
  border-bottom:1px solid rgba(255,255,255,.1)}
 #cp-top h3{margin:0;font-size:13px;color:#fff;letter-spacing:.05em}
 #cp-chk{font-size:10.5px;color:#94A3B8}
 #cp-chk b{color:#34D399}
 #cp-copy{margin-left:auto;background:rgba(139,156,249,.2);color:#fff;
  border:1px solid #8B9CF9;border-radius:8px;padding:5px 12px;
  font-size:11px;cursor:pointer}
 #cp-copy:hover{background:rgba(139,156,249,.35)}
 #cp-tip{position:absolute;display:none;background:rgba(15,23,42,.94);
  border:1px solid rgba(255,255,255,.2);border-radius:8px;color:#E2E8F0;
  font-size:11px;padding:5px 9px;pointer-events:none;z-index:6;
  max-width:240px}
 #cp-insp{position:absolute;width:300px;max-height:320px;overflow-y:auto;
  background:rgba(15,23,42,.94);border:1px solid rgba(255,255,255,.2);
  border-radius:14px;color:#E2E8F0;font-size:12px;padding:12px 14px;
  display:none;z-index:5;right:186px;top:56px}
 #cp-insp h4{margin:0 0 6px;font-size:13px;color:#fff}
 #cp-insp .x{position:absolute;right:9px;top:7px;cursor:pointer;
  color:#94A3B8}
 #cp-insp p{margin:5px 0;line-height:1.45}
 #cp-insp b{color:#fff}
 #cp-status{font-size:10.5px;color:#94A3B8}
</style>
<div id="cp-top"><h3>PACK COCKPIT</h3>
 <span id="cp-chk">surrogate check: ...</span>
 <span id="cp-status">live</span>
 <button id="cp-copy">&#10697; Copy settings for Design</button></div>
<div id="cp-grid">
 <div class="rail" id="railL"><h5>POWER &amp; ENVIRONMENT</h5></div>
 <div style="position:relative"><canvas id="cp-cv"></canvas>
   <div id="cp-tip"></div>
   <div id="cp-insp"><span class="x" id="cp-x">&times;</span>
     <div id="cp-insp-b"></div></div></div>
 <div class="rail" id="railR"><h5>COOLING SYSTEM</h5></div>
</div>
<div id="cp-pfd"></div>
</div>
<script>
(function(){
"use strict";
const P = __PAYLOAD__;
const status = document.getElementById('cp-status');
function fail(m){status.textContent='error: '+m;status.style.color='#FCA5A5';
 try{const c=document.getElementById('cp-cv').getContext('2d');
  c.setTransform(1,0,0,1,0,0);c.fillStyle='#FCA5A5';
  c.font='13px Inter';c.fillText('Cockpit error: '+m,14,26);}catch(e){}}
window.addEventListener('error',e=>fail(e.message));
try{
const D=P.design,G=P.geom,B=P.base,W_=P.water,K=P.consts;
const FL={};P.fluids.forEach(f=>FL[f.name]=f);
// ---------------- controls state ----------------
const C={c:D.C1,tamb:D.T_amb,fluid:D.coolant,flow:D.flow_lpm,
 twin:D.T_water_in,nt:D.n_tubes,circ:D.circ0,u:D.u0,plt:D.plate_t*1000,
 plc:D.plate_contact,fins:D.fins_on,pitch:D.pitch*1000};
// ---------------- physics port ----------------
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
function hCell(f,Ts,Tb,gapmm,u){const p=props(f,0.5*(Ts+Tb));
 const Ra=Ra_(p,Ts-Tb,D.h_cell);
 const hn=nuVert(Ra,p.Pr)*gapf(gapmm)*p.k/D.h_cell;
 let hf=0;const Re=u*D.d_cell/p.nu;
 if(u>1e-6)hf=nuCB(Re,p.Pr)*p.k/D.d_cell;
 return{h:blend(hn,hf),hn,hf,Re};}
function hTube(f,Tb,Tw,u){const p=props(f,0.5*(Tb+Tw));
 const Ra=Ra_(p,Tb-Tw,D.tube_od);
 const hn=nuHor(Ra,p.Pr)*p.k/D.tube_od;
 let hf=0;if(u>1e-6)hf=nuCB(u*D.tube_od/p.nu,p.Pr)*p.k/D.tube_od;
 return{h:blend(hn,hf)};}
function hWater(md,di,L){const{mu,k,cp}=W_;const Pr=mu*cp/k;
 const Re=md>0?4*md/(Math.PI*mu*di):0;
 const lam=R=>{const gz=(di/L)*R*Pr;
  return 3.66+0.0668*gz/(1+0.04*Math.pow(gz,2/3));};
 const tur=R=>{const f=Math.pow(0.790*Math.log(R)-1.64,-2);
  return (f/8)*(R-1000)*Pr/(1+12.7*Math.sqrt(f/8)*(Math.pow(Pr,2/3)-1));};
 let Nu,reg;
 if(Re<=0){Nu=3.66;reg='no flow';}
 else if(Re<2300){Nu=lam(Re);reg='laminar';}
 else if(Re<3000){const w=(Re-2300)/700;
  Nu=(1-w)*lam(2300)+w*tur(3000);reg='transitional';}
 else{Nu=tur(Re);reg='turbulent';}
 return{h:Nu*k/di,Re,regime:reg};}
function finPack(h){const r1=D.tube_od/2,r2=r1+D.fin_h,r2c=r2+D.fin_t/2;
 const n=1/D.fin_p,Abare=Math.PI*D.tube_od*Math.max(0,1-D.fin_t/D.fin_p);
 const Afin=n*2*Math.PI*(r2c*r2c-r1*r1);
 const m=Math.sqrt(2*Math.max(h,1)/(D.k_fin*D.fin_t));
 const Lc=D.fin_h+D.fin_t/2,phi=1+0.35*Math.log(r2c/r1),x=m*Lc*phi;
 const eta=x>1e-9?Math.tanh(x)/x:1;
 return{Aeff:Abare+eta*Afin,eta,gain:(Abare+eta*Afin)/(Math.PI*D.tube_od)};}
function plateFin(h,pitch_m,on,t,contact){if(!on)return{A:0,eta:0,m:0};
 const k=205,L=D.h_cell,
 m=Math.sqrt(2*Math.max(h,5)/(k*t)),eta=Math.tanh(m*L)/Math.max(m*L,1e-9);
 const npl=Math.max(G.n_rows-1,1),len=G.plate_len;
 return{A:npl*2*L*len*eta*contact,eta,
  m:npl*L*len*t*2700};}
function thermo(f,Q,Toil){const p=props(f,Toil);
 const A=G.A_flow,Dh=G.D_h,Ll=2.2*G.fill_h,H=G.H_loop,Kl=K.K_loop;
 const resid=u=>p.beta*g0*H*Q/(u*A*p.cp)
  -(32*p.rho*p.nu*Ll*u/(Dh*Dh)+Kl*0.5*p.rho*u*u);
 let lo=1e-6,hi=0.08;
 for(let i=0;i<50;i++){const m=0.5*(lo+hi);
  if(resid(m)>0)lo=m;else hi=m;}
 const u=lo;return{u,dT:Math.min(Q/(p.rho*u*A*p.cp),60)};}
function rOfT(T){return D.r_dc*Math.exp(-D.k_dcir*(T-25));}
const rCore=1/(4*Math.PI*D.k_rad*D.h_cell);
function waterPump(flow,nt){const md=flow/60*W_.rho/1000,
 mdt=md/Math.max(nt,1),Ai=Math.PI*G.d_i*G.d_i/4,
 v=mdt/(W_.rho*Ai),Re=W_.rho*v*G.d_i/W_.mu,
 f=Re<2300?64/Math.max(Re,1):0.316*Math.pow(Re,-0.25),
 dp=(f*G.L_tube/G.d_i+6)*0.5*W_.rho*v*v;
 return dp*(md/W_.rho)/0.35;}
function stirP(f,u){if(u<=1e-6)return 0;const p=props(f,35);
 const dp=32*p.rho*p.nu*(2.2*G.fill_h)*u/(G.D_h*G.D_h)
  +K.K_loop*0.5*p.rho*u*u;
 return dp*(u*G.A_flow)/0.30;}
function serpP(f,u,pitch_m,plt){if(u<=1e-6)return 0;const p=props(f,35);
 const s=Math.max((pitch_m-D.d_cell-plt)/2,5e-4),
 L=G.plate_len,nch=Math.max(G.n_rows-1,1)*2,Ach=s*D.h_cell,
 dp=12*p.rho*p.nu*L*u/(s*s)+3*0.5*p.rho*u*u,
 md=p.rho*u*Ach*nch;
 return dp*(md/p.rho)/0.35;}
function chiller(Qw,Tin,Tamb){const Tc=Tin-5+273.15,Th=Tamb+10+273.15,
 lift=Math.max(Th-Tc,3),COP=Math.max(0.45*Tc/lift,0.4);
 return{COP,Pel:Qw/COP,lift};}
// -------- full steady solve (60-pass, DCIR-coupled) --------
function solve(){
 const f=FL[C.fluid],pitch_m=C.pitch/1000,
 gapmm=Math.max((pitch_m-D.d_cell)*1000,0.3),
 nt=Math.round(C.nt),
 md=C.flow/60*W_.rho/1000,mdt=md/Math.max(nt,1);
 const wat=hWater(mdt,G.d_i,G.L_tube);
 const Rin=1/Math.max(wat.h*Math.PI*G.d_i*G.L_tube*nt,1e-9);
 const Rw=Math.log(D.tube_od/G.d_i)/(2*Math.PI*D.k_tube*G.L_tube*nt);
 const Ratm=1/Math.max(D.h_ext*G.A_box_ext,1e-9);
 const serp=C.circ==='serpentine',stir=C.circ==='stirred';
 const uc=(serp||stir)?C.u:0;
 let Twin=C.twin,Til=Twin+8,Tb=Til+6,Twl=Twin+2,Q=1000,uts=0,dTl=0;
 let Rb=0,Rot=0,Ao=0,pl={A:0,eta:0,m:0},hc={h:0},ht={h:0};
 const I=C.c*D.cap_Ah;
 for(let it=0;it<60;it++){
  Q=G.N*I*I*rOfT(Tb)*1e-3;
  Q+=I*B.I0*B.Rbus0;                    // busbars, J held constant
  const ts=thermo(f,Math.max(Q,1),0.5*(Tb+Til));uts=ts.u;dTl=ts.dT;
  const ue=Math.max(uc,uts);
  hc=hCell(f,Tb,Til,gapmm,ue);
  ht=hTube(f,Til,Twl,ue);
  hc.h*=K.cal;ht.h*=K.cal;
  let Aeff;
  if(C.fins){const fp=finPack(ht.h);Aeff=fp.Aeff*G.L_tube*nt;}
  else Aeff=Math.PI*D.tube_od*G.L_tube*nt;
  pl=plateFin(K.cal*ht.h,pitch_m,serp,C.plt/1000,C.plc);
  Ao=Aeff+pl.A;
  Rb=1/Math.max(hc.h*G.A_cells,1e-9);
  Rot=1/Math.max(ht.h*Ao,1e-9);
  const Rc=Rot+Rw+Rin;
  let dTr=Q/Math.max(md*W_.cp,1e-9);
  const Ts=Twin+0.5*Math.min(dTr,60);
  const TilN=(Q+Ts/Rc+C.tamb/Ratm)/(1/Rc+1/Ratm);
  const Qw=(TilN-Ts)/Rc;
  const TwlN=Ts+Qw*(Rin+Rw);
  const TbN=TilN+Q*Rb;
  Til+=0.6*(TilN-Til);Tb+=0.6*(TbN-Tb);Twl+=0.6*(TwlN-Twl);
 }
 const dTw=Q/Math.max(md*W_.cp,1e-9);
 const Qw=(Til-(Twin+0.5*dTw))/(Rot+Rw+Rin);
 const Qatm=(Til-C.tamb)/(1/Math.max(D.h_ext*G.A_box_ext,1e-9));
 if(serp||D.interstitial)dTl*=0.35;
 const spread=dTw+dTl;
 const Tcore=Tb+(Q/G.N)*rCore;
 const Pp=waterPump(C.flow,nt);
 const Pc=serp?serpP(f,C.u,pitch_m,C.plt/1000):stir?stirP(f,C.u):0;
 const chl=chiller(Math.max(Qw,1)+Pp,Twin,C.tamb);
 const films=[['can-oil film',Rb],['oil-tube film',Rot],
  ['water film',Rin]];
 films.sort((a,b)=>b[1]-a[1]);
 return{Tb,Tcore,Til,Q,Qw,Qatm,dTw,dTl,spread,uts,
  ue:Math.max(uc,uts),Rb,Rot,Rin,Rw,hc:hc.h,ht:ht.h,hw:wat.h,
  Re:wat.Re,regime:wat.regime,Ao,pl,Pp,Pc,chl,gapmm,
  weak:films[0][0],
  weakRegion:films[0][0]==='can-oil film'?'cells':'tubes'};
}
// -------- surrogate audit at the design point --------
let res=solve();
const dchk=res.Tb-B.T_b;
document.getElementById('cp-chk').innerHTML=
 'surrogate check at design point: <b>'+
 (dchk>=0?'+':'')+dchk.toFixed(2)+' °C</b> vs full solver'+
 (Math.abs(dchk)>0.4?' <span style="color:#F59E0B">(review)</span>':'');
// ---------------- control rails ----------------
const railL=document.getElementById('railL'),
      railR=document.getElementById('railR');
function slider(rail,key,label,min,max,step,fmt,onch){
 const d=document.createElement('div');d.className='ctl';
 d.innerHTML='<label>'+label+'<b id="v_'+key+'"></b></label>'+
  '<input type="range" min="'+min+'" max="'+max+'" step="'+step+
  '" value="'+C[key]+'">';
 rail.appendChild(d);
 const inp=d.querySelector('input'),vv=d.querySelector('#v_'+key);
 const upd=()=>{vv.textContent=fmt(C[key]);};
 inp.addEventListener('input',()=>{C[key]=parseFloat(inp.value);
  upd();(onch||recalc)();});
 upd();return d;}
function segctl(rail,key,label,opts,onch){
 const d=document.createElement('div');d.className='ctl';
 d.innerHTML='<label>'+label+'</label><div class="seg"></div>';
 const seg=d.querySelector('.seg');
 opts.forEach(o=>{const b=document.createElement('button');
  b.textContent=o[1];if(C[key]===o[0])b.classList.add('on');
  b.addEventListener('click',()=>{C[key]=o[0];
   seg.querySelectorAll('button').forEach(x=>x.classList.remove('on'));
   b.classList.add('on');(onch||recalc)();});
  seg.appendChild(b);});
 rail.appendChild(d);return d;}
function switchctl(rail,key,label){
 const d=document.createElement('div');d.className='ctl sw';
 d.innerHTML='<span>'+label+'</span><input type="checkbox" '+
  (C[key]?'checked':'')+'>';
 d.querySelector('input').addEventListener('change',e=>{
  C[key]=e.target.checked;recalc();});
 rail.appendChild(d);}
slider(railL,'c','C-rate (continuous)',0.2,6,0.05,v=>v.toFixed(2)+' C');
slider(railL,'tamb','Ambient',0,45,1,v=>v.toFixed(0)+' °C');
(function(){const d=document.createElement('div');d.className='ctl';
 let o='';P.fluids.forEach(f=>{o+='<option'+(f.name===C.fluid?
  ' selected':'')+'>'+f.name+'</option>';});
 d.innerHTML='<label>Dielectric fluid</label><select>'+o+'</select>';
 d.querySelector('select').addEventListener('change',e=>{
  C.fluid=e.target.value;recalc();});
 railL.appendChild(d);})();
slider(railL,'pitch','Cell pitch',23,33,0.5,
 v=>v.toFixed(1)+' mm (gap '+Math.max(v-D.d_cell*1000,0.3).toFixed(1)+')');
slider(railR,'flow','Water flow',0.5,60,0.5,v=>v.toFixed(1)+' L/min');
slider(railR,'twin','Water inlet',2,40,1,v=>v.toFixed(0)+' °C');
slider(railR,'nt','Tubes',2,48,1,v=>v.toFixed(0));
segctl(railR,'circ','Circulation',
 [['thermosiphon','THERM'],['stirred','STIR'],['serpentine','SERP']],
 ()=>{uRow.style.display=C.circ==='thermosiphon'?'none':'block';
  pRow.style.display=C.circ==='serpentine'?'block':'none';
  cRow.style.display=C.circ==='serpentine'?'block':'none';recalc();});
const uRow=slider(railR,'u','Velocity',0.005,0.15,0.005,
 v=>(v*1000).toFixed(0)+' mm/s');
const pRow=slider(railR,'plt','Plate thickness',1.0,2.0,0.5,
 v=>v.toFixed(1)+' mm');
const cRow=slider(railR,'plc','Plate contact',0.4,1.0,0.05,
 v=>v.toFixed(2));
switchctl(railR,'fins','Tube fins');
uRow.style.display=C.circ==='thermosiphon'?'none':'block';
pRow.style.display=C.circ==='serpentine'?'block':'none';
cRow.style.display=C.circ==='serpentine'?'block':'none';
// ---------------- PFD + annunciators ----------------
const pfd=document.getElementById('cp-pfd');
pfd.innerHTML=
 '<div class="pfd"><span>CELL CAN / CORE</span><b id="p_t"></b>'+
 '<div class="bar"><i id="p_tb"></i></div></div>'+
 '<div class="pfd"><span>MARGIN TO '+D.T_limit.toFixed(0)+
 ' °C</span><b id="p_m"></b><small>governing: '+
 (D.limit_core?'core':'can')+'</small></div>'+
 '<div class="pfd"><span>SPREAD</span><b id="p_s"></b>'+
 '<small>criterion 5 °C</small></div>'+
 '<div class="pfd"><span>HEAT</span><b id="p_q"></b>'+
 '<small id="p_qs"></small></div>'+
 '<div class="pfd"><span>CHILLER</span><b id="p_c"></b>'+
 '<small id="p_cs"></small></div>'+
 '<div class="pfd"><span>PARASITICS</span><b id="p_p"></b>'+
 '<small id="p_ps"></small></div>'+
 '<div class="pfd"><span>ADDED MASS</span><b id="p_k"></b>'+
 '<small>plates vs baseline</small></div>'+
 '<div class="ann">'+
 '<span class="lamp" id="l_lim">LIMIT</span>'+
 '<span class="lamp" id="l_lam">LAMINAR</span>'+
 '<span class="lamp" id="l_gap">GAP&lt;6</span>'+
 '<span class="lamp" id="l_spr">SPREAD&gt;5</span>'+
 '<span class="lamp" id="l_dry">DRY-COOLER OK</span></div>';
const $=id=>document.getElementById(id);
function recalc(){res=solve();
 const Tg=D.limit_core?res.Tcore:res.Tb,m=D.T_limit-Tg;
 $('p_t').textContent=res.Tb.toFixed(1)+' / '+res.Tcore.toFixed(1)+' °C';
 const fr=Math.min(Math.max((Tg-15)/(D.T_limit-15),0),1.15);
 $('p_tb').style.width=(fr*100).toFixed(0)+'%';
 $('p_tb').style.background=m<0?'#EF4444':m<3?'#F59E0B':'#34D399';
 $('p_m').textContent=(m>=0?'+':'')+m.toFixed(1)+' °C';
 $('p_m').style.color=m<0?'#F87171':m<3?'#FCD34D':'#6EE7B7';
 $('p_s').textContent=res.spread.toFixed(1)+' °C';
 $('p_q').textContent=(res.Q/1000).toFixed(2)+' kW';
 $('p_qs').textContent='at '+C.c.toFixed(2)+'C · casing '+
  res.Qatm.toFixed(0)+' W';
 $('p_c').textContent=(res.chl.Pel/1000).toFixed(2)+' kW el';
 $('p_cs').textContent='duty '+((res.Qw+res.Pp)/1000).toFixed(2)+
  ' kW · COP '+res.chl.COP.toFixed(1);
 $('p_p').textContent=(res.Pp+res.Pc).toFixed(1)+' W';
 $('p_ps').textContent='pump '+res.Pp.toFixed(1)+
  (res.Pc>0?' + circ '+res.Pc.toFixed(1):'');
 $('p_k').textContent=(res.pl.m).toFixed(1)+' kg';
 $('l_lim').className='lamp'+(m<0?' red':'');
 $('l_lam').className='lamp'+(res.regime==='laminar'?' amb':'');
 $('l_gap').className='lamp'+(res.gapmm<6?' amb':'');
 $('l_spr').className='lamp'+(res.spread>5?' amb':'');
 $('l_dry').className='lamp'+(C.twin>=C.tamb+5?' grn':'');
}
recalc();
// ---------------- copy bridge ----------------
document.getElementById('cp-copy').addEventListener('click',()=>{
 const out=JSON.stringify({c1:C.c,tamb:C.tamb,fluid:C.fluid,
  flow:C.flow,twin:C.twin,ntub:Math.round(C.nt),circ:C.circ,
  u:C.u,plt:C.plt,plc:C.plc,fins:C.fins,pitch:C.pitch});
 const done=()=>{status.textContent='settings copied - paste below the '+
  'cockpit to apply';setTimeout(()=>status.textContent='live',3500);};
 if(navigator.clipboard&&navigator.clipboard.writeText)
  navigator.clipboard.writeText(out).then(done,()=>prompt('Copy:',out));
 else prompt('Copy:',out);});
// ---------------- canvas: live pack drawing ----------------
const cv=document.getElementById('cp-cv'),cx=cv.getContext('2d');
const tip=document.getElementById('cp-tip'),
 insp=document.getElementById('cp-insp'),
 inspB=document.getElementById('cp-insp-b');
document.getElementById('cp-x').addEventListener('click',
 ()=>{sel=null;insp.style.display='none';});
let Wc=760,Hc=520;
function fit(){Wc=cv.clientWidth||760;Hc=cv.clientHeight||520;
 const d=window.devicePixelRatio||1;
 cv.width=Math.round(Wc*d);cv.height=Math.round(Hc*d);
 cx.setTransform(d,0,0,d,0,0);}
fit();if(window.ResizeObserver)new ResizeObserver(fit).observe(cv);
let R={},sel=null,hover=null,mx=-1,my=-1,tsec=0;
function geom(){const mL=20,mR=20,top=26,bot=40;
 const bx=mL,by=top,bw=Wc-mL-mR,bh=Hc-top-bot;
 const oilTop=by+bh*(1-G.fill_frac);
 const cTop=by+bh*(1-G.cell_top_frac),cBot=by+bh*(1-G.cell_bot_frac);
 const tubeY=D.interstitial?(cTop+cBot)/2:Math.max(oilTop+22,by+26);
 R={bx,by,bw,bh,oilTop,cTop,cBot,tubeY,n:G.n_rows_draw,
  nt:Math.min(Math.round(C.nt),14)};
 R.pitch=bw/(R.n+0.6);R.cw=R.pitch*(D.d_cell/(C.pitch/1000));
 R.tp=bw/(R.nt+1);}
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
 if(C.circ==='serpentine'){const l=Math.floor(py*R.n);
  return{vx:(l%2?-1:1)*u*2.4,vy:Math.sin(px*14+l*2)*2.5};}
 const k=C.circ==='stirred'?1:3,A=u*2.2;
 return{vx:A*Math.sin(k*Math.PI*px)*Math.cos(Math.PI*py),
  vy:-A*Math.cos(k*Math.PI*px)*Math.sin(Math.PI*py)};}
const oilP=[],watP=[];
for(let i=0;i<160;i++)oilP.push({x:Math.random(),y:Math.random()});
for(let i=0;i<70;i++)watP.push({x:Math.random(),lane:i%12});
function hitTest(x,y){
 for(let j=0;j<R.nt;j++){const tx=R.bx+R.tp*(j+1);
  if(Math.hypot(x-tx,y-R.tubeY)<15)return 'tubes';}
 if(y>=R.cTop&&y<=R.cBot&&x>R.bx&&x<R.bx+R.bw)return 'cells';
 if(x>R.bx&&x<R.bx+R.bw&&y>R.oilTop&&y<R.by+R.bh)return 'oil';
 if(x>R.bx&&x<R.bx+R.bw&&y>R.by&&y<R.oilTop)return 'head';
 return null;}
function stText(id){
 if(id==='cells')return{t:'Cells',h:'can '+res.Tb.toFixed(1)+
  ' °C · core '+res.Tcore.toFixed(1)+' °C',
  b:'<p>Each cell: <b>'+(res.Q/G.N).toFixed(2)+' W</b> at '+
  C.c.toFixed(2)+'C. First film h = <b>'+res.hc.toFixed(0)+
  ' W/m²·K</b> -> '+(res.Q*res.Rb).toFixed(1)+' °C toll. DCIR '+
  rOfT(res.Tb).toFixed(1)+' mΩ ('+
  (100*(1-rOfT(res.Tb)/D.r_dc)).toFixed(0)+
  '% below 25 °C).</p>'};
 if(id==='oil')return{t:'Bulk oil',h:res.Til.toFixed(1)+
  ' °C · '+(res.ue*1000).toFixed(1)+' mm/s',
  b:'<p>Oil at <b>'+res.Til.toFixed(1)+' °C</b>, circulation '+
  (res.ue*1000).toFixed(1)+' mm/s ('+C.circ+'). Spread <b>'+
  res.spread.toFixed(1)+' °C</b>; casing sheds '+
  res.Qatm.toFixed(0)+' W free.</p>'};
 if(id==='tubes')return{t:'Tubes, fins'+
  (C.circ==='serpentine'?', plates':''),
  h:'oil h '+res.ht.toFixed(0)+' · Re '+res.Re.toFixed(0)+
  ' ('+res.regime+')',
  b:'<p>Second film: h = '+res.ht.toFixed(0)+', A = <b>'+
  res.Ao.toFixed(1)+' m²</b>'+(res.pl.A>0?' (plates +'+
  res.pl.A.toFixed(1)+' m², η '+res.pl.eta.toFixed(2)+')':'')+
  ' -> '+(res.Q*res.Rot).toFixed(1)+' °C. Water: <b>'+
  res.regime+'</b>, Re '+res.Re.toFixed(0)+' -> '+
  (res.Q*res.Rin).toFixed(1)+' °C; stream warms '+
  res.dTw.toFixed(1)+' °C.</p>'};
 return{t:'Headspace',h:'gas blanket',
  b:'<p>Nitrogen blanket over the oil; burst rating sizes the '+
  'enclosure. Change nothing here from the cockpit.</p>'};}
cv.addEventListener('mousemove',e=>{const r=cv.getBoundingClientRect();
 mx=e.clientX-r.left;my=e.clientY-r.top;hover=hitTest(mx,my);
 if(hover){const s=stText(hover);
  tip.innerHTML='<b>'+s.t+'</b><br>'+s.h+
   '<br><span style="color:#94A3B8">click to inspect</span>';
  tip.style.display='block';
  tip.style.left=Math.min(mx+14,Wc-250)+'px';
  tip.style.top=Math.min(my+12,Hc-70)+'px';cv.style.cursor='pointer';}
 else{tip.style.display='none';cv.style.cursor='crosshair';}});
cv.addEventListener('mouseleave',()=>{tip.style.display='none';
 hover=null;});
cv.addEventListener('click',()=>{if(!hover){sel=null;
  insp.style.display='none';return;}
 sel=hover;const s=stText(sel);
 inspB.innerHTML='<h4>'+s.t+'</h4>'+s.b;insp.style.display='block';});
let last=performance.now();
function frame(now){try{
 const dt=Math.min((now-last)/1000,0.05);last=now;tsec+=dt;
 geom();const{bx,by,bw,bh,oilTop,cTop,cBot,tubeY,n,nt,pitch,cw,tp}=R;
 cx.fillStyle='#0B1220';cx.fillRect(0,0,Wc,Hc);
 cx.strokeStyle='rgba(255,255,255,.45)';cx.lineWidth=2;
 cx.strokeRect(bx-6,by-6,bw+12,bh+12);
 cx.fillStyle='#0E1830';cx.fillRect(bx,by,bw,oilTop-by);
 cx.fillStyle='rgba(245,158,11,0.10)';
 cx.fillRect(bx,oilTop,bw,by+bh-oilTop);
 cx.strokeStyle='rgba(245,158,11,.6)';cx.lineWidth=1;
 cx.beginPath();cx.moveTo(bx,oilTop);cx.lineTo(bx+bw,oilTop);cx.stroke();
 for(let i=0;i<n;i++){const x=bx+pitch*(0.4+i)+(pitch-cw)/2;
  const gcx=x+cw/2,gcy=(cTop+cBot)/2;
  const gl=cx.createRadialGradient(gcx,gcy,2,gcx,gcy,cw*1.7);
  const a=Math.min(0.12+res.Q/1000*0.05,0.4);
  gl.addColorStop(0,'rgba(241,82,82,'+a+')');
  gl.addColorStop(1,'rgba(241,82,82,0)');
  cx.fillStyle=gl;
  cx.fillRect(x-cw,cTop-cw,cw*3,(cBot-cTop)+2*cw);
  cx.fillStyle=tcol(res.Tb,0.96);
  cx.strokeStyle='rgba(255,255,255,.30)';cx.lineWidth=1;
  rrect(x,cTop,cw,cBot-cTop,5,true,true);
  cx.fillStyle=tcol(res.Tcore,0.95);
  rrect(x+cw*0.30,cTop+4,cw*0.40,(cBot-cTop)-8,4,true,false);
  if(C.circ==='serpentine'&&i<n-1){
   cx.fillStyle='rgba(139,156,249,.9)';
   cx.fillRect(x+cw+(pitch-cw)/2-1.2,cTop,2.4,cBot-cTop);}}
 for(let j=0;j<nt;j++){const tx=bx+tp*(j+1);
  cx.fillStyle='rgba(165,180,204,.4)';
  cx.beginPath();cx.arc(tx,tubeY,11,0,6.283);cx.fill();
  cx.fillStyle='#D97706';
  cx.beginPath();cx.arc(tx,tubeY,5.5,0,6.283);cx.fill();}
 cx.fillStyle='#7DD3FC';cx.font='10px Inter';
 cx.fillText('in '+C.twin.toFixed(0)+'°C',bx+4,tubeY-16);
 cx.fillStyle=tcol(C.twin+res.dTw,1);
 cx.fillText('out '+(C.twin+res.dTw).toFixed(1)+'°C',
  bx+bw-64,tubeY-16);
 for(const p of watP){p.x+=dt*Math.min(C.flow/20,2)*0.25;
  if(p.x>1)p.x-=1;
  const tx=bx+tp*((p.lane%nt)+1),ang=p.x*6.283;
  cx.fillStyle=tcol(C.twin+res.dTw*p.x,0.95);
  cx.beginPath();
  cx.arc(tx+Math.cos(ang)*3.4,tubeY+Math.sin(ang)*3.4,1.8,0,6.283);
  cx.fill();}
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
  const tx=bx+tp*(Math.floor(nt/2)+1);
  cx.beginPath();cx.arc(tx,tubeY,17+3*pulse,0,6.283);cx.stroke();
  cx.fillStyle='rgba(248,113,113,.95)';cx.font='600 11px Inter';
  cx.fillText('weakest: '+res.weak,tx+24,tubeY+4);}
 else{const i=Math.floor(n/2),
  x=bx+pitch*(0.4+i)+(pitch-cw)/2;
  rrect(x-3-2*pulse,cTop-3-2*pulse,cw+6+4*pulse,
   (cBot-cTop)+6+4*pulse,7,false,true);
  cx.fillStyle='rgba(248,113,113,.95)';cx.font='600 11px Inter';
  cx.fillText('weakest: '+res.weak,x+cw+10,cTop+14);}
 const lx=bx+bw/2-70,ly=by+bh+14;
 const lo=Math.min(C.twin,C.tamb),
  hi=Math.max(res.Tcore,res.Tb+1,lo+8);
 for(let i=0;i<140;i++){cx.fillStyle=tcol(lo+(hi-lo)*i/140,1);
  cx.fillRect(lx+i,ly,1,8);}
 cx.fillStyle='#94A3B8';cx.font='9.5px Inter';
 cx.fillText(lo.toFixed(0)+'°C',lx-26,ly+8);
 cx.fillText(hi.toFixed(0)+'°C',lx+146,ly+8);
}catch(e){fail(e.message);return;}
 requestAnimationFrame(frame);}
requestAnimationFrame(frame);
}catch(e){fail(e.message);}
})();
</script>
"""
    return tpl.replace("__PAYLOAD__", P)
