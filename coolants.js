// Representative properties near 25 C, 1 atm (liquid metals at typical operating temp).
// Values are approximate and for teaching; real properties vary strongly with
// temperature, grade, and additives. Sources: standard heat-transfer references.
// cat: water | glycol | oil | dielectric | metal | gas
// rho kg/m3, cp J/kgK, k W/mK, mu Pa*s, bp/fp degC (null = not applicable)
window.COOLANTS = [
  {name:"Water",                       cat:"water",      rho:997,  cp:4180, k:0.606, mu:0.00089, bp:100, fp:0,    diel:false, use:"Best all-round performance. Conductive, corrosive, freezes."},
  {name:"Deionized water",             cat:"water",      rho:997,  cp:4180, k:0.606, mu:0.00089, bp:100, fp:0,    diel:"part",use:"Electronics loops; resistivity falls as it ages and leaches ions."},
  {name:"Ethylene glycol 50% (EGW)",   cat:"glycol",     rho:1070, cp:3300, k:0.37,  mu:0.0039,  bp:107, fp:-37,  diel:false, use:"Automotive and industrial antifreeze. Toxic."},
  {name:"Propylene glycol 50% (PGW)",  cat:"glycol",     rho:1040, cp:3400, k:0.36,  mu:0.006,   bp:106, fp:-33,  diel:false, use:"Non-toxic antifreeze for food, HVAC, and solar loops."},
  {name:"Ethylene glycol (neat)",      cat:"glycol",     rho:1110, cp:2400, k:0.252, mu:0.0161,  bp:197, fp:-13,  diel:false, use:"High boiling point; very poor without water dilution."},
  {name:"Mineral / transformer oil",   cat:"oil",        rho:860,  cp:1860, k:0.13,  mu:0.025,   bp:300, fp:-40,  diel:true,  use:"Immersion cooling and transformers. Cheap, stable, dielectric."},
  {name:"Silicone oil (5 cSt)",        cat:"oil",        rho:920,  cp:1600, k:0.117, mu:0.0046,  bp:140, fp:-65,  diel:true,  use:"Low-viscosity grade; wide temperature range, dielectric."},
  {name:"Silicone oil (50 cSt)",       cat:"oil",        rho:960,  cp:1550, k:0.16,  mu:0.048,   bp:200, fp:-55,  diel:true,  use:"Mid-viscosity grade. Pumping cost rises sharply."},
  {name:"Silicone oil (1000 cSt)",     cat:"oil",        rho:970,  cp:1550, k:0.16,  mu:0.97,    bp:250, fp:-50,  diel:true,  use:"Very viscous. Shows how viscosity destroys convective performance."},
  {name:"PAO synthetic oil",           cat:"oil",        rho:800,  cp:2200, k:0.14,  mu:0.012,   bp:260, fp:-50,  diel:true,  use:"Dielectric coolant for avionics and electronics; wide temp range."},
  {name:"Glycerol (glycerin)",         cat:"oil",        rho:1260, cp:2430, k:0.285, mu:0.95,    bp:290, fp:18,   diel:"part",use:"Extremely viscous and hygroscopic; mostly an instructive extreme."},
  {name:"Therminol / Dowtherm",        cat:"oil",        rho:1060, cp:1550, k:0.137, mu:0.004,   bp:257, fp:-30,  diel:"part",use:"Synthetic aromatic for high-temperature process heat transfer."},
  {name:"FC-72 (perfluorocarbon)",     cat:"dielectric", rho:1680, cp:1100, k:0.057, mu:0.00064, bp:56,  fp:-90,  diel:true,  use:"Two-phase immersion of electronics. Low CHF, high GWP."},
  {name:"HFE-7100 (Novec)",            cat:"dielectric", rho:1510, cp:1183, k:0.069, mu:0.00061, bp:61,  fp:-135, diel:true,  use:"Two-phase electronics; lower environmental impact than PFCs."},
  {name:"Galinstan (GaInSn)",          cat:"metal",      rho:6440, cp:370,  k:16.5,  mu:0.0024,  bp:1300,fp:-19,  diel:false, use:"Liquid-metal cooling for ultra-high heat flux. Conductive, wets metals."},
  {name:"NaK eutectic",                cat:"metal",      rho:866,  cp:982,  k:22,    mu:0.00047, bp:785, fp:-12,  diel:false, use:"Fast-reactor coolant. Extremely reactive with air and water."},
  {name:"Liquid sodium",               cat:"metal",      rho:927,  cp:1230, k:71,    mu:0.00069, bp:883, fp:98,   diel:false, use:"Fast-reactor coolant (props at ~100 C). Reactive; solid at room temp."},
  {name:"Air (reference)",             cat:"gas",        rho:1.18, cp:1005, k:0.026, mu:0.0000185,bp:null,fp:null, diel:true,  use:"Gas baseline. Included to show why liquids win by orders of magnitude."}
];

// Prandtl number: ratio of momentum to thermal diffusivity.
window.prandtl = c => c.mu * c.cp / c.k;

// Mouromtseff number: figure of merit for single-phase turbulent forced convection
// (higher is better). Mo = rho^0.8 * k^0.6 * cp^0.4 / mu^0.4
window.mouromtseff = c =>
  Math.pow(c.rho,0.8) * Math.pow(c.k,0.6) * Math.pow(c.cp,0.4) / Math.pow(c.mu,0.4);
