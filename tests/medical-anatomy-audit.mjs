// Run from repository root: node tests/medical-anatomy-audit.mjs
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';

const html=readFileSync('index.html','utf8');
const script=html.match(/<script>([\s\S]*?)<\/script>/)?.[1];
assert.ok(script,'inline JS missing');
new Function(script);

const match=script.match(/const MEDICAL_ANATOMY=([\s\S]*?);\nlet medicalFocus=/);
assert.ok(match,'medical anatomy reference data missing');
const data=vm.runInNewContext('('+match[1]+')');

const keys=['bicepsLong','bicepsShort','brachialis','tricepsLong','tricepsLateral','tricepsMedial','deltoid','supraspinatus','infraspinatus','teresMinor','teresMajor','subscapularis','coracobrachialis'];
assert.deepEqual(Object.keys(data),keys);
for(const [key,row] of Object.entries(data)){
 for(const field of ['name','origin','insertion','nerve','action','note','source'])assert.ok(row[field]?.length, key+' missing '+field);
 assert.ok(Array.isArray(row.parts),key+' parts not an array');
 assert.ok(/^https:\/\/www\.ncbi\.nlm\.nih\.gov\/books\/NBK\d+\/$/.test(row.source),key+' reference URL invalid');
}
assert.equal(data.tricepsMedial.parts.length,0,'do not pretend a medial-head mesh exists');
assert.equal(data.coracobrachialis.parts.length,0,'do not pretend coracobrachialis mesh exists');
assert.deepEqual(data.subscapularis.parts,['v58Subscap'],'subscapularis has atlas loader ID');
assert.ok(script.includes("file:'FMA13414.stl'"));
const staticParts=[...script.matchAll(/const (?:MUSCLE_PARTS|TENDON_PARTS)=\[([\s\S]*?)\];/g)]
 .flatMap(m=>[...m[1].matchAll(/'([^']+)'/g)].map(v=>v[1]));
const scapSpecMatch=script.match(/const SCAPULAR_SPEC=(\[[\s\S]*?\]);\nconst SCAPULAR_CUFF=/);
assert.ok(scapSpecMatch,'scapular muscle specification missing');
const scapSpec=vm.runInNewContext(scapSpecMatch[1]);
assert.equal(scapSpec.length,4);
assert.equal(scapSpec.filter(s=>s.cuff).length,3);
const scapParts=scapSpec.flatMap(m=>[m.id,'tendon'+m.id[0].toUpperCase()+m.id.slice(1)]);
for(const row of Object.values(data))for(const part of row.parts)
 assert.ok(staticParts.includes(part)||scapParts.includes(part)||part==='v58Subscap','unknown model part '+part);
const htmlWithoutScript=html.replace(/<script>[\s\S]*?<\/script>/g,'');
const ids=[...htmlWithoutScript.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1]);
assert.equal(ids.length,new Set(ids).size,'duplicate HTML ids');
for(const id of ['medicalStudyPanel','medicalStudySelect','medicalStudyFields','medicalStudySource','medicalStudyFocus','medicalStudyReset','medicalStudyModelStatus']){
 assert.ok(ids.includes(id),'study element missing: '+id);
}
assert.ok(script.includes('medicalFocusPrevious'),'previous visibility state should be restored');
assert.ok(script.includes("if(medicalFocus&&MUSCLE_ALL.includes(name)"),'focus filter missing');
assert.ok(html.includes('医学部生向け教材としての形態・運動学的検証は未完了'),'validation limits missing');
assert.ok(/上肢3D解剖ラボ · v\d+\.\d+/.test(html),'v0.55 version missing');
assert.ok(script.includes('1+.26*contractionRate*exaggeration'),'bulge indicator must match biceps visual coefficient');
console.log('PASS v0.55: medical anatomy regression, 13 facts records + missing muscle honesty + source links + syntax');
