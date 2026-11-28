import assert from 'node:assert/strict';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
const read=p=>JSON.parse(readFileSync(p,'utf8'));
const base=read('test-results/shoulder-muscle-control/shoulder-muscle-report.json'),candidate=read('test-results/shoulder-muscle-preview/shoulder-muscle-report.json');
assert.equal(base.preview,false);assert.equal(candidate.preview,true);assert.equal(base.results.length,candidate.results.length);
const comparisons=[],regressions=[];
for(let i=0;i<base.results.length;i++){
 const a=base.results[i],b=candidate.results[i];assert.deepEqual(a.input,b.input);
 for(const group of ['active','other','neck']){
  const x=group==='neck'?a.neck:a[group].muscles,y=group==='neck'?b.neck:b[group].muscles;
  assert.equal(x.length,y.length);
  for(let j=0;j<x.length;j++){
   const c=x[j],d=y[j];assert.equal(c.id,d.id);assert.equal(c.vertices,d.vertices);assert.equal(c.triangles,d.triangles);assert.equal(c.attachmentSamples,d.attachmentSamples);
   const metrics={edgeMax:[c.edgeRatio[1],d.edgeRatio[1]],edgeMin:[c.edgeRatio[0],d.edgeRatio[0]],areaMax:[c.areaRatio[1],d.areaRatio[1]],areaMin:[c.areaRatio[0],d.areaRatio[0]],attachment:[c.attachmentMaxMm,d.attachmentMaxMm],degenerate:[c.newDegenerateFaces,d.newDegenerateFaces],invalidFrames:[c.invalidFrames,d.invalidFrames]};
   const row={pose:a.input.name,group,id:c.id,metrics,volumeRatio:[c.volumeRatio,d.volumeRatio],physicalVolumeVerified:false};
   comparisons.push(row);
   for(const [name,pair]of Object.entries(metrics)){
    if(pair.some(v=>v===null))continue;
    assert(pair.every(Number.isFinite),'Nonfinite metric '+name);
    const signed=name.endsWith('Min')?-1:1;
    const delta=signed*(pair[1]-pair[0]);
    if(delta>1e-6*Math.max(1,Math.abs(pair[0])))regressions.push({pose:row.pose,group,id:row.id,metric:name,before:pair[0],after:pair[1],delta});
   }
  }
 }
}
regressions.sort((a,b)=>b.delta-a.delta);
const report={cases:base.results.length,comparisons:comparisons.length,regressionCount:regressions.length,regressions,rows:comparisons,previewHash:candidate.previewHash,productAccepted:false,anatomicalAcceptance:false,reason:'Existing muscle deformation remains unvalidated; this is a paired surface/attachment regression audit, not a physical volume or collision certificate.'};
mkdirSync('test-results/shoulder-muscle-comparison',{recursive:true});
writeFileSync('test-results/shoulder-muscle-comparison/report.json',JSON.stringify(report,null,2)+'\n');
console.log('MUSCLE_COMPARISON '+JSON.stringify({cases:report.cases,comparisons:report.comparisons,regressionCount:report.regressionCount,worst:regressions.slice(0,12),productAccepted:false}));
