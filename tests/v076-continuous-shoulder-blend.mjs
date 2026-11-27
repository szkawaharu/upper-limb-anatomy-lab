import assert from 'node:assert/strict';
import {readBlend} from '../scripts/calibration/continuous-shoulder-blend.mjs';
const f=readBlend(),s=f.sample,h=.0001;
for(let i=0;i<f.planes.length;i++)for(let j=0;j<f.nodes.length;j++)assert.equal(s(f.nodes[j],f.planes[i]),f.table[i][j]);
for(let e=0;e<=180;e++)for(let p=0;p<=90;p++){
 const a=s(e,p);assert(a>=0&&a<=.25+1e-12);
 if(e===0)assert.equal(a,0);
}
for(const p of [0,20,40,65,90])for(const e of f.nodes.slice(1,-1))assert(Math.abs((s(e+h,p)-2*s(e,p)+s(e-h,p))/h)<1e-5);
for(let e=0;e<=180;e+=3)assert(Math.abs((s(e,40+h)-2*s(e,40)+s(e,40-h))/h)<1e-5);
assert.throws(()=>s(-1,0));assert.throws(()=>s(0,91));
assert.equal(f.productAccepted,false);
console.log('PASS continuous candidate coefficient: evidence nodes, bounded C1 interpolation, neutral, no extrapolation; not bone collision acceptance');
