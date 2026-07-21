"""render_cockpit.py - render the cockpit's canvas EXACTLY as a browser
would, by executing the component's JavaScript in node with node-canvas
supplying a real 2D context. Produces PNGs of the main pack view (and
the instrument strip) for any circulation mode and tube shape, so
drawing changes can be reviewed without a browser.

Usage:
    CP_DUMP=1 SMOKE=1 python app.py        # once, to dump the payload
    python render_cockpit.py serpentine round /tmp/serp.png
    python render_cockpit.py extpump rect /tmp/ext.png
Needs node plus the npm 'canvas' package on NODE_PATH.
"""
import json
import re
import subprocess
import sys

NODE_STUB = r"""
'use strict';
const {createCanvas} = require('canvas');
const REAL = {};
REAL['cp-cv']   = createCanvas(1120, 560);
REAL['cp-inst'] = createCanvas(1120, 96);
function _el(id){
 const e={_id:id,style:{},dataset:{},children:[],
  value:'0',checked:false,textContent:'',innerHTML:'',
  clientWidth:id==='cp-cv'?1120:(id==='cp-inst'?1120:300),
  clientHeight:id==='cp-cv'?560:(id==='cp-inst'?96:40),
  classList:{add(){},remove(){},toggle(){},contains(){return false;}},
  addEventListener(){},removeEventListener(){},
  appendChild(c){this.children.push(c);return c;},
  insertBefore(c){this.children.push(c);return c;},
  removeChild(){},setAttribute(){},getAttribute(){return null;},
  querySelector(){return _el();},querySelectorAll(){return [];},
  getBoundingClientRect(){return{left:0,top:0,width:this.clientWidth,
   height:this.clientHeight,right:this.clientWidth,
   bottom:this.clientHeight};},
  focus(){},blur(){},click(){}};
 if(REAL[id]){
  const cnv=REAL[id];
  Object.defineProperty(e,'width',{get:()=>cnv.width,
   set:v=>{cnv.width=v;}});
  Object.defineProperty(e,'height',{get:()=>cnv.height,
   set:v=>{cnv.height=v;}});
  e.getContext=()=>cnv.getContext('2d');
 } else {
  e.getContext=()=>new Proxy({},{
   get(t,k){
    if(k==='createLinearGradient'||k==='createRadialGradient')
     return ()=>({addColorStop(){}});
    if(k==='measureText')return ()=>({width:12});
    if(typeof k==='string')return ()=>{};
    return undefined;},
   set(){return true;}});
 }
 return e;}
const _els={};
const document={
 getElementById(id){return _els[id]||(_els[id]=_el(id));},
 createElement(){return _el();},
 querySelector(){return _el();},
 querySelectorAll(){return [];},
 addEventListener(){},
 body:_el('body'),documentElement:_el('html')};
let _raf=[];
const window={devicePixelRatio:1,addEventListener(){},
 requestAnimationFrame(cb){_raf.push(cb);},
 location:{reload(){}},innerWidth:1280,innerHeight:800};
const requestAnimationFrame=cb=>_raf.push(cb);
const navigator={clipboard:{writeText(){return Promise.resolve();}}};
const location={reload(){}};
const alert=()=>{},prompt=()=>'';
"""

NODE_DRIVE = r"""
;(function(){
 const fs=require('fs');
 let t=0;
 for(let i=0;i<__NFRAMES__;i++){
  const q=_raf.slice();_raf.length=0;t+=33.4;
  for(const cb of q)cb(t);
 }
 fs.writeFileSync('__OUT__', REAL['cp-cv'].toBuffer('image/png'));
 fs.writeFileSync('__OUT_INST__',
  REAL['cp-inst'].toBuffer('image/png'));
 console.log('RENDER_OK');
})();
"""


def render(circ="serpentine", shape="round", out="/tmp/cockpit.png",
           nframes=40, design_overrides=None):
    src = open("cockpit.py").read()
    js = re.search(r"<script>(.*)</script>", src, re.S).group(1)
    js = js.replace("catch(e){fail(e.message);return;}",
                    "catch(e){throw e;}")
    js = js.replace("catch(e){fail(e.message);}", "catch(e){throw e;}")
    pay = json.load(open("/tmp/cp_payload.json"))
    p = json.loads(json.dumps(pay))
    p["design"]["circ0"] = circ
    p["design"]["tshape"] = shape
    p["design"].setdefault("tw", 0.012)
    p["design"].setdefault("th", 0.008)
    p["design"].setdefault("pipe_id", 0.019)
    p["design"].setdefault("pipe_len", 2.5)
    p["design"].setdefault("ver", "dev")
    for k, v in (design_overrides or {}).items():
        p["design"][k] = v
    body = js.replace("__PAYLOAD__", json.dumps(p))
    drive = (NODE_DRIVE.replace("__NFRAMES__", str(nframes))
             .replace("__OUT_INST__", out.replace(".png", "_inst.png"))
             .replace("__OUT__", out))
    open("/tmp/_cp_render.js", "w").write(NODE_STUB + body + drive)
    r = subprocess.run(["node", "/tmp/_cp_render.js"],
                       capture_output=True, text=True, timeout=120,
                       env={"NODE_PATH": "/tmp/nodecanvas/node_modules",
                            "PATH": "/usr/bin:/bin:/usr/local/bin"})
    if r.returncode != 0 or "RENDER_OK" not in r.stdout:
        sys.exit("render failed:\n" + (r.stderr or r.stdout)[:1500])
    print(f"rendered {circ}/{shape} -> {out}")


if __name__ == "__main__":
    a = sys.argv[1:]
    render(a[0] if a else "serpentine",
           a[1] if len(a) > 1 else "round",
           a[2] if len(a) > 2 else "/tmp/cockpit.png")
