"""js_runtime_check.py - actually EXECUTE the cockpit JavaScript in
node with a stub DOM, once per circulation mode and tube shape, so
runtime-only faults (undefined functions, bad property chains) are
caught before shipping. `node --check` only parses; this runs init, a
recalc, and three animation frames per configuration, with the
component's own try/catch converted to re-throws so failures exit
non-zero.

Needs: node on PATH, and /tmp/cp_payload.json from
    CP_DUMP=1 SMOKE=1 python app.py
Run:  python js_runtime_check.py
"""
import json
import re
import subprocess
import sys

STUB = r"""
'use strict';
function _el(id){const e={_id:id,style:{},dataset:{},children:[],
 value:'0',checked:false,textContent:'',innerHTML:'',
 clientWidth:760,clientHeight:430,width:0,height:0,
 classList:{add(){},remove(){},toggle(){},contains(){return false;}},
 addEventListener(){},removeEventListener(){},
 appendChild(c){this.children.push(c);return c;},
 insertBefore(c){this.children.push(c);return c;},
 removeChild(){},setAttribute(){},getAttribute(){return null;},
 querySelector(){return _el();},querySelectorAll(){return [];},
 getBoundingClientRect(){return{left:0,top:0,width:760,height:430,
  right:760,bottom:430};},
 focus(){},blur(){},click(){}};
 e.getContext=function(){return CTX;};
 return e;}
const CTX=new Proxy({},{
 get(t,k){
  if(k==='createLinearGradient'||k==='createRadialGradient')
   return ()=>({addColorStop(){}});
  if(k==='measureText')return ()=>({width:12});
  if(k==='getTransform')return ()=>({a:1,b:0,c:0,d:1,e:0,f:0});
  if(typeof k==='string')return ()=>{};
  return undefined;},
 set(){return true;}});
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

DRIVE = r"""
;(function(){
 let t=0;
 for(let i=0;i<3;i++){
  const q=_raf.slice();_raf.length=0;t+=16.7;
  for(const cb of q)cb(t);
 }
 console.log('RUNTIME_OK');
})();
"""


def main():
    src = open("cockpit.py").read()
    js = re.search(r"<script>(.*)</script>", src, re.S).group(1)
    # make internal failures fatal so node's exit code reflects them
    js = js.replace("catch(e){fail(e.message);return;}",
                    "catch(e){throw e;}")
    js = js.replace("catch(e){fail(e.message);}", "catch(e){throw e;}")
    pay = json.load(open("/tmp/cp_payload.json"))
    runs = [("thermosiphon", "round"), ("stirred", "round"),
            ("prop", "round"), ("serpentine", "round"),
            ("serpentine", "square"), ("serpentine", "rect"),
            ("extpump", "round"), ("extpump", "rect")]
    fails = 0
    for circ, shape in runs:
        p = json.loads(json.dumps(pay))
        p["design"]["circ0"] = circ
        p["design"]["tshape"] = shape
        p["design"].setdefault("tw", 0.012)
        p["design"].setdefault("th", 0.008)
        p["design"].setdefault("pipe_id", 0.019)
        p["design"].setdefault("pipe_len", 2.5)
        body = js.replace("__PAYLOAD__", json.dumps(p))
        open("/tmp/_cp_run.js", "w").write(STUB + body + DRIVE)
        r = subprocess.run(["node", "/tmp/_cp_run.js"],
                           capture_output=True, text=True, timeout=60)
        ok = r.returncode == 0 and "RUNTIME_OK" in r.stdout
        print(f"{circ:12s} {shape:6s} "
              f"{'OK' if ok else 'FAIL'}")
        if not ok:
            fails += 1
            err = (r.stderr or r.stdout).strip().splitlines()
            print("   " + "\n   ".join(err[:8]))
    if fails:
        sys.exit(f"{fails} runtime configuration(s) FAILED")
    print("ALL RUNTIME CONFIGS PASS")


if __name__ == "__main__":
    main()
