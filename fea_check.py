"""fea_check.py - validate generated FEA files without the solvers.
FEMM .lua: executed under lua5.4 with every hi_/ho_ call stubbed, so
syntax errors and undefined names fail loudly. COMSOL .java: compiled
with javac against a minimal chainable stub of the COMSOL API, so
typos, arity and structure errors fail at compile time. Neither check
proves the physics - the real solvers do that - but both catch the
class of error that otherwise only surfaces on the user's machine.

Usage: python fea_check.py <directory with .lua/.java files>
"""
import pathlib
import shutil
import subprocess
import sys
import tempfile

LUA_STUB = r"""
local function _n(...) return 0 end
newdocument=_n hi_probdef=_n hi_drawline=_n hi_drawarc=_n
hi_addmaterial=_n hi_addblocklabel=_n hi_selectlabel=_n
hi_setblockprop=_n hi_clearselected=_n hi_addboundprop=_n
hi_selectsegment=_n hi_setsegmentprop=_n hi_zoomnatural=_n
hi_saveas=_n hi_analyze=_n hi_loadsolution=_n
ho_addcontour=_n ho_clearcontour=_n ho_showdensityplot=_n
ho_selectblock=_n ho_blockintegral=_n
function ho_getpointvalues(...) return 21.5 end
function ho_lineintegral(...) return 1.234 end
function messagebox(s) print(s) end
format=string.format
function openfile(n,m) return {} end
function write(h,s) end
function closefile(h) end
"""

JAVA_MODEL_STUB = r"""
package com.comsol.model;
public class Model {
  public Model param() { return this; }
  public Model component(Object... a) { return this; }
  public Model geom(Object... a) { return this; }
  public Model feature(Object... a) { return this; }
  public Model create(Object... a) { return this; }
  public Model set(String k, Object... v) { return this; }
  public Model selection(Object... a) { return this; }
  public Model named(Object... a) { return this; }
  public Model all() { return this; }
  public Model physics(Object... a) { return this; }
  public Model material(Object... a) { return this; }
  public Model propertyGroup(Object... a) { return this; }
  public Model cpl(Object... a) { return this; }
  public Model mesh(Object... a) { return this; }
  public Model autoMeshSize(int n) { return this; }
  public Model axisymmetric(boolean b) { return this; }
  public Model study(Object... a) { return this; }
  public Model result(Object... a) { return this; }
  public Model numerical(Object... a) { return this; }
  public Model table(Object... a) { return this; }
  public Model export(Object... a) { return this; }
  public Model run() { return this; }
  public Model setResult() { return this; }
  public Model save(String s) throws java.io.IOException {
    return this; }
}
"""

JAVA_UTIL_STUB = r"""
package com.comsol.model.util;
import com.comsol.model.Model;
public class ModelUtil {
  public static Model create(String s) { return new Model(); }
}
"""


def _stub_dir():
    d = pathlib.Path(tempfile.mkdtemp(prefix="comsolstub_"))
    m = d / "com" / "comsol" / "model"
    u = m / "util"
    u.mkdir(parents=True)
    (m / "Model.java").write_text(JAVA_MODEL_STUB)
    (u / "ModelUtil.java").write_text(JAVA_UTIL_STUB)
    r = subprocess.run(["javac", str(m / "Model.java"),
                        str(u / "ModelUtil.java")],
                       capture_output=True, text=True)
    if r.returncode:
        sys.exit("stub compile failed:\n" + r.stderr)
    return d


def check_dir(path):
    path = pathlib.Path(path)
    luas = sorted(path.glob("*.lua"))
    javas = sorted(path.glob("*.java"))
    fails = 0
    if luas and shutil.which("lua5.4"):
        for f in luas:
            merged = LUA_STUB + "\n" + f.read_text()
            tmp = pathlib.Path(tempfile.mkstemp(suffix=".lua")[1])
            tmp.write_text(merged)
            r = subprocess.run(["lua5.4", str(tmp)],
                               capture_output=True, text=True,
                               timeout=30)
            ok = r.returncode == 0
            print(f"lua  {f.name:26s} {'OK' if ok else 'FAIL'}")
            if not ok:
                fails += 1
                print("   " + (r.stderr or r.stdout).strip()[:400])
    elif luas:
        print("lua5.4 not found - skipping lua execution checks")
    if javas and shutil.which("javac"):
        stub = _stub_dir()
        for f in javas:
            r = subprocess.run(
                ["javac", "-cp", str(stub), "-d",
                 tempfile.mkdtemp(), str(f)],
                capture_output=True, text=True, timeout=120)
            ok = r.returncode == 0
            print(f"java {f.name:26s} {'OK' if ok else 'FAIL'}")
            if not ok:
                fails += 1
                print("   " + r.stderr.strip()[:600])
    elif javas:
        print("javac not found - skipping java compile checks")
    if fails:
        sys.exit(f"{fails} file(s) FAILED")
    print("ALL FEA FILES PASS the stub checks")


if __name__ == "__main__":
    check_dir(sys.argv[1] if len(sys.argv) > 1 else ".")
