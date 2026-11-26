import assert from'node:assert/strict';
import '../anatomy/rest-attachments-v076.js';
const R=RestAttachmentsV076;
assert.equal(Object.keys(R.choices).length,17);
const refined=Float32Array.of(0,.125,.375,.75,1),old=[.1,.25,.5,.8,.9];
for(const [id,choice]of Object.entries(R.choices)){
 let calls=0;const smooth=(a,n)=>{calls++;assert.equal(n,6);return a.slice()},a=R.field(id,old,refined,smooth),b=R.field(id,old,refined,smooth);
 assert.deepEqual(a,b);assert.equal(a[0],0);assert.equal(a[4],1);assert(a.every(w=>w>=0&&w<=1));
 assert.equal(calls,choice==='smooth-midpoint'?2:0);
}
assert.throws(()=>R.field('new-unreviewed-muscle',old,refined,()=>{}),/Missing versioned/);
console.log('17 immutable rest-binding policies; original fixed masks retained; skeletal solver is not an input');
