"""Builds explorer.html: a self-contained static viewer (no external requests)."""
import sys, os, json, glob
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
import cities

KEEP = ["fitness", "cost", "n_lines", "n_stations", "track_km", "mean_stations_per_line", "n_transfer_stations", "rings",
        "mode_share", "mean_headway", "mean_transfers_per_trip", "longest_share", "water_hops"]
runs = []
for p in sorted(glob.glob(os.path.join(ROOT, "results", "runs", "*.json"))):
    r = json.load(open(p))
    runs.append(dict(city=r["city"], budget=r["budget"], pen=r["pen"], mode=r["mode"], seed=r["seed"], units=round(r["budget_units"], 1),
                     lines=r["lines"], m={k: round(r["metrics"][k], 4) for k in KEEP}, phase=[round(x) for x in r["phase_fit"]]))
real = json.load(open(os.path.join(ROOT, "results", "real_eval.json")))
cityd = {}
for n in cities.CITY_ORDER:
    c = cities.get_city(n)
    dens = np.log1p((c.H + c.J) / c.cell ** 2)
    q = np.clip((dens / dens.max() * 9).round(), 0, 9).astype(int)
    cityd[n] = dict(name=c.label, R=c.R, C=c.C, cell=c.cell, budget=c.mid_budget, desc=c.desc,
                    water=["".join("1" if w else "0" for w in row) for row in c.water],
                    dens=["".join(str(v) for v in row) for row in q],
                    real={k: [l for l in v["lines"]] for k, v in (real.get(n, {}) or {}).items() if v["pen"] == "lowpen"})
    cityd[n]["real"] = {v["label"]: v["lines"] for k, v in (real.get(n, {}) or {}).items() if v["pen"] == "lowpen"}
    cityd[n]["realfit"] = {v["label"]: round(v["metrics"]["fitness"]) for k, v in (real.get(n, {}) or {}).items() if v["pen"] == "lowpen"}

HTML = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Subway network explorer</title>
<style>
:root{--bg:#fff;--fg:#1a1a1a;--mut:#666;--card:#f4f4f6;--bd:#d0d0d8}
@media(prefers-color-scheme:dark){:root{--bg:#16171a;--fg:#e8e8ea;--mut:#9a9aa3;--card:#202227;--bd:#3a3c44}}
body{margin:0;padding:16px;background:var(--bg);color:var(--fg);font:14px/1.4 system-ui,sans-serif;max-width:1200px;margin:auto}
h1{font-size:18px;margin:0 0 4px}.sub{color:var(--mut);margin-bottom:12px}
.ctl{display:flex;flex-wrap:wrap;gap:10px;margin-bottom:10px}
label{font-size:12px;color:var(--mut);display:flex;flex-direction:column;gap:2px}
select{padding:4px;background:var(--card);color:var(--fg);border:1px solid var(--bd);border-radius:4px}
.row{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:10px;margin-bottom:14px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:6px;padding:6px}
canvas{width:100%;height:auto;display:block;border-radius:4px}
table{border-collapse:collapse;font-size:11px;width:100%}td{padding:1px 3px;border-top:1px solid var(--bd)}td:last-child{text-align:right}
h2{font-size:14px;margin:6px 0 4px}.note{font-size:11px;color:var(--mut)}
</style></head><body>
<h1>Subway network explorer</h1>
<div class="sub">End-state evolved networks. Pick a condition (A) and optionally a second one (B) to compare across the 3 seeds. Real-city geography is a stylized approximation. White dots: transfer stations. Fitness = person-hours/day of generalized time saved vs driving.</div>
<div class="ctl" id="ctlA"></div>
<h2>Condition A</h2><div class="row" id="rowA"></div>
<div class="ctl" id="ctlB"></div>
<h2>Condition B</h2><div class="row" id="rowB"></div>
<div class="note" id="foot"></div>
<script>
const D=__DATA__;
const PAL=["#d62728","#1f77b4","#2ca02c","#ff7f0e","#9467bd","#8c564b","#e377c2","#17becf","#bcbd22","#7f7f7f"];
const cityKeys=Object.keys(D.cities);
function mk(id,lab,opts,val){const l=document.createElement('label');l.textContent=lab;const s=document.createElement('select');
 opts.forEach(o=>{const e=document.createElement('option');e.value=o[0];e.textContent=o[1];s.appendChild(e)});s.value=val;l.appendChild(s);document.getElementById(id).appendChild(l);return s}
function ctl(id,init,showReal){
 const c=mk(id,'City',cityKeys.map(k=>[k,D.cities[k].name]),init.city);
 const b=mk(id,'Budget',[['low','low (0.5x)'],['mid','mid (1x)'],['high','high (2x)']],init.budget);
 const p=mk(id,'Transfer penalty',[['lowpen','low (5 min)'],['highpen','high (20 min)']],init.pen);
 const m=mk(id,'Build mode',[['planned','planned'],['phased','phased (3 tranches)']],init.mode);
 const r=mk(id,'Overlay real network',[['none','none']],'none');
 const o={c,b,p,m,r};
 function refreshReal(){const cc=D.cities[c.value];r.innerHTML='';[['none','none']].concat(Object.keys(cc.real).map(k=>[k,k+' (fit '+cc.realfit[k]+')'])).forEach(x=>{const e=document.createElement('option');e.value=x[0];e.textContent=x[1];r.appendChild(e)})}
 c.addEventListener('change',refreshReal);refreshReal();return o}
const A=ctl('ctlA',{city:'baltimore',budget:'mid',pen:'lowpen',mode:'planned'});
const B=ctl('ctlB',{city:'baltimore',budget:'mid',pen:'lowpen',mode:'phased'});
function draw(cv,city,lines,realLines){
 const s=city.cell,W=city.C,H=city.R,k=6;cv.width=W*k*2;cv.height=H*k*2;const g=cv.getContext('2d');g.scale(2,2);
 for(let r=0;r<H;r++)for(let c=0;c<W;c++){
  if(city.water[r][c]=='1'){g.fillStyle='#8fb8dd'}else{const v=+city.dens[r][c];const x=Math.round(250-v*17);g.fillStyle='rgb('+x+','+x+','+x+')'}
  g.fillRect(c*k,r*k,k,k)}
 const P=(cell)=>[(cell[1]+.5)*k,(cell[0]+.5)*k];
 if(realLines){g.setLineDash([3,2]);g.lineWidth=1.4;g.strokeStyle='#000';realLines.forEach(l=>{g.beginPath();l.forEach((c,i)=>{const[p,q]=P(c);i?g.lineTo(p,q):g.moveTo(p,q)});g.stroke()});g.setLineDash([])}
 const cnt={};
 lines.forEach((l,i)=>{g.strokeStyle=PAL[i%PAL.length];g.lineWidth=2.2;g.lineCap='round';g.beginPath();l.forEach((c,j)=>{const[p,q]=P(c);j?g.lineTo(p,q):g.moveTo(p,q);cnt[c]=(cnt[c]||0)+1});g.stroke()});
 Object.keys(cnt).forEach(key=>{const c=key.split(',').map(Number);const[p,q]=P(c);g.beginPath();g.arc(p,q,cnt[key]>1?2.8:1.6,0,7);g.fillStyle=cnt[key]>1?'#fff':'#000';g.fill();if(cnt[key]>1){g.strokeStyle='#000';g.lineWidth=1;g.stroke()}})}
function find(o,seed){return D.runs.find(r=>r.city==o.c.value&&r.budget==o.b.value&&r.pen==o.p.value&&r.mode==o.m.value&&r.seed==seed)}
function render(o,rowId){
 const row=document.getElementById(rowId);row.innerHTML='';const city=D.cities[o.c.value];
 const real=o.r.value!='none'?city.real[o.r.value]:null;
 [0,1,2].forEach(sd=>{const r=find(o,sd);const card=document.createElement('div');card.className='card';
  if(!r){card.textContent='no run';row.appendChild(card);return}
  const cv=document.createElement('canvas');card.appendChild(cv);draw(cv,city,r.lines,real);
  const t=document.createElement('table');const m=r.m;
  const rows=[['seed',sd],['fitness (person-h/day)',Math.round(m.fitness)],['cost / budget',m.cost.toFixed(0)+' / '+r.units],['lines',m.n_lines],['stations',m.n_stations],['track km',m.track_km.toFixed(1)],
   ['stations per line',m.mean_stations_per_line.toFixed(1)],['transfer stations',m.n_transfer_stations],['rings',m.rings],['mode share',(m.mode_share*100).toFixed(1)+'%'],
   ['mean headway (min)',m.mean_headway.toFixed(1)],['transfers per trip',m.mean_transfers_per_trip.toFixed(2)],['longest line share',(m.longest_share*100).toFixed(0)+'%'],['water-crossing hops',m.water_hops]];
  if(r.mode=='phased')rows.push(['fitness after phase 1/2/3',r.phase.join(' / ')]);
  t.innerHTML=rows.map(x=>'<tr><td>'+x[0]+'</td><td>'+x[1]+'</td></tr>').join('');card.appendChild(t);row.appendChild(card)})}
function all(){render(A,'rowA');render(B,'rowB')}
[A,B].forEach(o=>Object.values(o).forEach(s=>s.addEventListener('change',all)));all();
document.getElementById('foot').textContent='Cities: '+cityKeys.map(k=>D.cities[k].name+' (mid budget '+D.cities[k].budget+' units)').join('; ')+'. Dashed black: selected real/proposed comparison network.';
</script></body></html>"""
html = HTML.replace("__DATA__", json.dumps(dict(runs=runs, cities=cityd), separators=(",", ":")))
open(os.path.join(ROOT, "explorer.html"), "w").write(html)
print("explorer.html", len(html) // 1024, "KB,", len(runs), "runs")
