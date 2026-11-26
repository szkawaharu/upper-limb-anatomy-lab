import assert from'node:assert/strict';
for(const f of['engine/math','engine/fiber-path','engine/volume-constraint','engine/contact-constraint','engine/tendon-deformation','engine/muscle-deformation','scripts/experimental/differential-chart'])await import('../'+f+'.js');
const M=AnatomyMath,D=DifferentialChart,rest=Float32Array.of(0,0,0,1,0,0,.3,.5,.2,.7,.8,-.2),faces=Uint32Array.of(0,2,1,0,1,3,0,3,2,1,2,3),line=Array.from({length:4},(_,i)=>[i/3,0,0]);
for(const coincident of[false,true]){
 const B=D.prepare({rest,faces,weights:Float32Array.of(0,1,.3,.7),lane:Float32Array.of(0,0,.5,.8),material:Float32Array.of(0,1,.3,.7),field:{controls:[line,line.map(p=>M.add(p,[0,coincident?0:1,0]))]},reference:[0,1,0],sheet:true,blendRange:[0,1],definition:{source:{worldUnitsPerMm:1},tissues:{tendonStart:.86}}});
 for(const angle of[0,.7,Math.PI,5.5]){
  const bone={q:M.qaxis([.3,.4,.5],angle),p:[.3,-.4,.2],rest:[.2,.3,-.5]},s=D.pose(B,bone,bone),out=D.deform(s);
  for(let i=0;i<B.weights.length;i++)assert(M.len(M.sub(out.subarray(i*3,i*3+3),M.apply(bone,rest.subarray(i*3,i*3+3))))<2e-7,'Rigid motion changed source shape');
 }
 const a=M.identity(),b={q:M.qaxis([0,0,1],1.8),p:[1,.1,0],rest:[1,0,0]},state=D.pose(B,a,b),out=D.deform(state);
 for(const [i,bone]of[[0,a],[1,b]])assert(M.len(M.sub(D.evaluate(state,i).point,M.apply(bone,rest.subarray(i*3,i*3+3))))<1e-12);
 assert([...out].every(Number.isFinite));
}
console.log('Experimental differential chart: rigid shape invariance, coincident lanes, finite output, exact fixed boundaries. Actual source-mesh shape remains rejected by independent audits.');
