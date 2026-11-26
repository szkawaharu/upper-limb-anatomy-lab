/* Experimental source-embedded ARAP. No clinical material constants, force or
 * activation are inferred. Not enabled in the application until validation.
 * Sorkine/Alexa 2007 local/global energy; Sumner et al. 2007 node embedding.
 * Unlike the former, this graph uses positive inverse-edge weights.
 */
(function (root) {
  'use strict';
  const M=root.AnatomyMath;
  function matrix(q,out,offset=0){
    const [x,y,z,w]=q,xx=2*x*x,yy=2*y*y,zz=2*z*z,xy=2*x*y,xz=2*x*z,yz=2*y*z,wx=2*w*x,wy=2*w*y,wz=2*w*z;
    out.set([1-yy-zz,xy-wz,xz+wy,xy+wz,1-xx-zz,yz-wx,xz-wy,yz+wx,1-xx-yy],offset);
  }
  // Davenport/Horn symmetric eigensystem. Jacobi returns the largest algebraic
  // eigenvector, avoiding the ± eigenvalue ambiguity of unshifted power steps.
  function rotation(H,out=new Float64Array(9),K=new Float64Array(16),V=new Float64Array(16),qout=null){
    const [xx,xy,xz,yx,yy,yz,zx,zy,zz]=H;
    K.set([xx+yy+zz,yz-zy,zx-xz,xy-yx,yz-zy,xx-yy-zz,xy+yx,zx+xz,zx-xz,xy+yx,-xx+yy-zz,yz+zy,xy-yx,zx+xz,yz+zy,-xx-yy+zz]);
    V.fill(0);for(let i=0;i<4;i++)V[i*4+i]=1;
    for(let sweep=0;sweep<10;sweep++){
      let biggest=0,scale=0;for(let i=0;i<4;i++)scale=Math.max(scale,Math.abs(K[i*4+i]));
      for(let p=0;p<3;p++)for(let q=p+1;q<4;q++){
        const a=K[p*4+q];biggest=Math.max(biggest,Math.abs(a));if(Math.abs(a)<=Math.max(scale,1e-30)*1e-14)continue;
        const u=(K[q*4+q]-K[p*4+p])/(2*a),t=(u<0?-1:1)/(Math.abs(u)+Math.sqrt(1+u*u)),c=1/Math.sqrt(1+t*t),s=t*c;
        const pp=K[p*4+p],qq=K[q*4+q];K[p*4+p]=pp-t*a;K[q*4+q]=qq+t*a;K[p*4+q]=K[q*4+p]=0;
        for(let k=0;k<4;k++)if(k!==p&&k!==q){const ap=K[k*4+p],aq=K[k*4+q];K[k*4+p]=K[p*4+k]=c*ap-s*aq;K[k*4+q]=K[q*4+k]=s*ap+c*aq;}
        for(let k=0;k<4;k++){const vp=V[k*4+p],vq=V[k*4+q];V[k*4+p]=c*vp-s*vq;V[k*4+q]=s*vp+c*vq;}
      }
      if(biggest<=Math.max(scale,1e-30)*1e-13)break;
    }
    let j=0;for(let k=1;k<4;k++)if(K[k*4+k]>K[j*4+j])j=k;
    const q=[V[4+j],V[8+j],V[12+j],V[j]];if(qout){const sign=q.reduce((s,x,k)=>s+x*qout[k],0)<0?-1:1;for(let k=0;k<4;k++)qout[k]=q[k]*sign;}matrix(q,out);return out;
  }
  function rotate(R,o,x,y,z,out,k=0){out[k]=R[o]*x+R[o+1]*y+R[o+2]*z;out[k+1]=R[o+3]*x+R[o+4]*y+R[o+5]*z;out[k+2]=R[o+6]*x+R[o+7]*y+R[o+8]*z;}
  function prepare(G,embedding){
    const n=G.nodeWeights.length,free=Int32Array.from(G.free),index=new Int32Array(n);index.fill(-1);free.forEach((x,i)=>index[x]=i);
    const factor={};for(const [k,v]of Object.entries(G.factor))factor[k]=k.includes('Indices')||k.includes('Offsets')?Int32Array.from(v):Float64Array.from(v);
    return {...G,n,nodes:Float64Array.from(G.nodes),nodeWeights:Float64Array.from(G.nodeWeights),nodeClass:Int8Array.from(G.nodeClass),edges:Int32Array.from(G.edges),edgeWeights:Float64Array.from(G.edgeWeights),free,index,factor,embedding,
      x:new Float64Array(n*3),R:new Float64Array(n*9),q:new Float64Array(n*4),covariance:new Float64Array(n*9),rhs:new Float64Array(n*3),solveBuffer:new Float64Array(free.length*3),K:new Float64Array(16),V:new Float64Array(16),rotationBuffer:new Float64Array(9)};
  }
  function solveFactor(G){
    const F=G.factor,N=G.free.length,z=G.solveBuffer;
    for(let i=0;i<N;i++)for(let k=0;k<3;k++){
      let r=G.rhs[G.free[i]*3+k];for(let j=F.lowerOffsets[i];j<F.lowerOffsets[i+1];j++)r-=F.lowerValues[j]*z[F.lowerIndices[j]*3+k];z[i*3+k]=r/F.diagonal[i];
    }
    for(let i=N-1;i>=0;i--)for(let k=0;k<3;k++){
      let r=z[i*3+k];for(let j=F.upperOffsets[i];j<F.upperOffsets[i+1];j++)r-=F.upperValues[j]*z[F.upperIndices[j]*3+k];z[i*3+k]=r/F.diagonal[i];G.x[G.free[i]*3+k]=z[i*3+k];
    }
  }
  function dualBlend(a,b,w,p){
    const dual=bone=>M.qmul([...M.sub(bone.p,M.qrot(bone.q,bone.rest)),0],bone.q).map(x=>x*.5),da=dual(a),db=dual(b),q=a.q.map((x,i)=>(1-w)*x+w*b.q[i]),d=da.map((x,i)=>(1-w)*x+w*db[i]),n=Math.hypot(...q);
    if(n<1e-8)return null;
    const real=q.map(x=>x/n),t=M.qmul(d.map(x=>x/n),[-real[0],-real[1],-real[2],real[3]]);
    return M.add(M.qrot(real,p),t.slice(0,3).map(x=>2*x));
  }
  function seed(G,a,b,originAt,initializer){
    G.initializationFallbacks=0;
    for(let i=0;i<G.n;i++){
      const C=G.nodes.subarray(i*3,i*3+3),A=originAt?originAt(C):a,w=G.nodeWeights[i],cls=G.nodeClass[i],ap=M.apply(A,C),bp=M.apply(b,C);
      // A geometric initializer only; the free graph then solves ARAP. The
      // continuous bone quaternion lift avoids a changing shortest-arc branch.
      const q=cls===0?A.q:cls===1?b.q:M.qblendLift(A.q,b.q,w);
      let x=cls===0?ap:cls===1?bp:null;
      if(!x&&(cls===3||initializer==='dqs')){x=dualBlend(A,b,w,C);if(!x)G.initializationFallbacks++;}
      matrix(q,G.R,i*9);G.q.set(q,i*4);G.x.set(x||M.mix(ap,bp,w),i*3);
    }
  }
  function local(G){
    const C=G.nodes,X=G.x,H=G.covariance;H.fill(0);
    for(let e=0;e<G.edgeWeights.length;e++){
      const i=G.edges[e*2],j=G.edges[e*2+1],w=G.edgeWeights[e];
      for(let a=0;a<3;a++)for(let b=0;b<3;b++){const h=w*(C[i*3+a]-C[j*3+a])*(X[i*3+b]-X[j*3+b]);H[i*9+a*3+b]+=h;H[j*9+a*3+b]+=h;}
    }
    for(const i of G.free){rotation(H.subarray(i*9,i*9+9),G.rotationBuffer,G.K,G.V,G.q.subarray(i*4,i*4+4));G.R.set(G.rotationBuffer,i*9);}
  }
  function globalStep(G){
    const C=G.nodes,X=G.x,R=G.R,Y=G.rhs;Y.fill(0);
    for(let e=0;e<G.edgeWeights.length;e++){
      const i=G.edges[e*2],j=G.edges[e*2+1],w=G.edgeWeights[e],v0=C[i*3]-C[j*3],v1=C[i*3+1]-C[j*3+1],v2=C[i*3+2]-C[j*3+2];
      for(let k=0;k<3;k++){
        const a=i*9+k*3,b=j*9+k*3,value=w*.5*((R[a]+R[b])*v0+(R[a+1]+R[b+1])*v1+(R[a+2]+R[b+2])*v2);
        Y[i*3+k]+=value;Y[j*3+k]-=value;
        if(G.index[j]<0)Y[i*3+k]+=w*X[j*3+k];if(G.index[i]<0)Y[j*3+k]+=w*X[i*3+k];
      }
    }
    solveFactor(G);
  }
  function pose(binding,a,b,options={}){
    const G=binding.graph,iterations=options.iterations??48,start=performance.now();seed(G,a,b,options.originAt,binding.graphInitializer);
    const state={binding,a,b,originAt:options.originAt,contact:options.contact||null,graph:G,iterations,solveMs:null,initializationFallbacks:G.initializationFallbacks,activation:null,force:null,vertexDiagnostics:()=>({stretch:null,determinant:null})};
    for(let k=0;k<iterations;k++){
      local(G);globalStep(G);
      if(options.obstacle)for(const i of G.free){const C=G.nodes.subarray(i*3,i*3+3),p=G.x.subarray(i*3,i*3+3),r=M.len(M.sub(C,options.obstacle.restCenter)),v=M.sub(p,options.obstacle.center),d=M.len(v),minimum=Math.min(r,options.obstacle.radius);if(d<minimum&&d>1e-10)G.x.set(M.add(options.obstacle.center,M.mul(v,minimum/d)),i*3);}
      if(options.volume)projectVolume(state,3);
    }
    local(G);
    if(options.volume)projectVolume(state,8);state.solveMs=performance.now()-start;
    // Cached poses must not alias the next solve, including the other arm.
    state.graph={...G,x:G.x.slice(),R:G.R.slice(),q:G.q.slice()};return state;
  }
  function point(state,index){
    const B=state.binding,G=state.graph,p=B.rest.subarray(index*3,index*3+3),w=B.weights[index];
    if(w===0||w===1)return M.apply(w?state.b:state.originAt?state.originAt(p):state.a,p);
    const result=[0,0,0],tmp=[0,0,0],E=G.embedding,o=index*8,center=[0,0,0],q=[0,0,0,0];
    for(let k=0;k<4;k++){
      const weight=E[o+4+k];if(!weight)continue;
      const i=E[o+k],c=i*3;
      if(B.graphEmbedding==='orthogonal'){
        for(let j=0;j<3;j++){center[j]+=weight*G.nodes[c+j];result[j]+=weight*G.x[c+j];}for(let j=0;j<4;j++)q[j]+=weight*G.q[i*4+j];
      }else{
        rotate(G.R,i*9,p[0]-G.nodes[c],p[1]-G.nodes[c+1],p[2]-G.nodes[c+2],tmp);for(let j=0;j<3;j++)result[j]+=weight*(tmp[j]+G.x[c+j]);
      }
    }
    if(B.graphEmbedding==='orthogonal'){const n=Math.hypot(...q);if(n<1e-8)throw Error('Ambiguous embedded rotation');const offset=M.qrot(q.map(x=>x/n),M.sub(p,center));for(let j=0;j<3;j++)result[j]+=offset[j];}
    // A centroid constraint alone does not fix an entire source footprint.
    // This normalized quartic bubble is zero with zero derivative at both
    // Dirichlet ends; its maximum is one. It is a common basis function, not
    // a muscle/pose fitted coefficient. The original source masks stay exact.
    if(B.graphBoundary==='quartic-dirichlet-bubble'){
      const a=state.originAt?state.originAt(p):state.a,base=M.mix(M.apply(a,p),M.apply(state.b,p),w),gate=16*w*w*(1-w)*(1-w);
      for(let j=0;j<3;j++)result[j]=base[j]+gate*(result[j]-base[j]);
    }
    return state.contact?root.ContactConstraint.project(p,result,w,state.contact):result;
  }
  function deform(state,out=new Float32Array(state.binding.rest.length)){
    for(let i=0;i<out.length/3;i++)out.set(point(state,i),i*3);return out;
  }
  function surfaceAudit(state,now){
    const B=state.binding,G=state.graph,E=G.embedding;let opposed=0,ambiguous=0,minCos=1,worst=null;
    for(let k=0;k<B.faces.length;k+=3){
      const ids=B.faces.subarray(k,k+3),p=Array.from(ids,i=>B.rest.subarray(i*3,i*3+3)),x=Array.from(ids,i=>now.subarray(i*3,i*3+3)),N=M.cross(M.sub(p[1],p[0]),M.sub(p[2],p[0])),posed=M.cross(M.sub(x[1],x[0]),M.sub(x[2],x[0])),expected=[0,0,0],tmp=[0,0,0];
      for(const i of ids){const w=B.weights[i];if(w===0||w===1){const a=w?state.b:state.originAt?state.originAt(B.rest.subarray(i*3,i*3+3)):state.a,v=M.qrot(a.q,N);for(let j=0;j<3;j++)expected[j]+=v[j]/3;}else for(let j=0;j<4;j++){const node=E[i*8+j],weight=E[i*8+4+j]/3;rotate(G.R,node*9,...N,tmp);for(let c=0;c<3;c++)expected[c]+=weight*tmp[c];}}
      if(M.len(expected)<M.len(N)*1e-9||M.len(posed)<1e-15){ambiguous++;continue;}
      const cos=M.dot(expected,posed)/(M.len(expected)*M.len(posed));if(cos<0)opposed++;if(cos<minCos){minCos=cos;worst={triangle:k/3,vertices:[...ids],sourceVertexIds:B.sourceVertexIds?Array.from(ids,i=>B.sourceVertexIds[i]):null,cos};}
    }
    return {method:'posed triangle normal versus blended source-face rotation; orientation warning only, not a volumetric-Jacobian certificate',opposedFaces:opposed,ambiguousFaces:ambiguous,minCos,worst};
  }
  function volumeData(state,gradient=false){
    const B=state.binding,Q=B.graphQuadrature,G=state.graph;if(!Q)return null;
    const origin=M.apply(state.a,Q.sourceCenter),pts=Q.sourceVertices.map(i=>M.sub(point(state,i),origin)),grad=gradient?new Float64Array(G.n*3):null;
    let volume=Q.referenceVolume+M.dot(M.sub(M.apply(state.b,Q.sourceCenter),origin),M.qrot(state.b.q,Q.pureBVector))/6;
    for(let i=0;i<Q.multipliers.length;i++){
      const face=Q.faces.slice(i*3,i*3+3),p=face.map(j=>pts[j]),mult=Q.multipliers[i]/6;
      volume+=mult*(M.dot(p[0],M.cross(p[1],p[2]))-Q.controlDet[i]);
      if(gradient)for(let j=0;j<3;j++){
        const id=Q.sourceVertices[face[j]],w=B.weights[id];if(w===0||w===1)continue;
        const v=M.cross(p[(j+1)%3],p[(j+2)%3]),o=id*8,gate=B.graphBoundary==='quartic-dirichlet-bubble'?16*w*w*(1-w)*(1-w):1;
        for(let k=0;k<4;k++){const n=G.embedding[o+k],s=G.embedding[o+4+k]*mult*gate;if(G.index[n]<0)continue;for(let c=0;c<3;c++)grad[n*3+c]+=s*v[c];}
      }
    }
    return {volume,grad};
  }
  function projectVolume(state,steps){
    const G=state.graph,Q=state.binding.graphQuadrature;if(!Q)return;
    const F=G.factor,N=G.free.length,z=G.solveBuffer,saved=new Float64Array(G.x.length);
    for(let s=0;s<steps;s++){
      const {volume,grad}=volumeData(state,true),error=Q.referenceVolume-volume;
      state.volumeProjection={ratio:volume/Q.referenceVolume,method:'global source-integral constraint in graph energy metric; local incompressibility unproven'};
      if(Math.abs(error)<Q.referenceVolume*1e-6)return;
      for(let i=0;i<N;i++)for(let k=0;k<3;k++){let r=grad[G.free[i]*3+k];for(let j=F.lowerOffsets[i];j<F.lowerOffsets[i+1];j++)r-=F.lowerValues[j]*z[F.lowerIndices[j]*3+k];z[i*3+k]=r/F.diagonal[i];}
      for(let i=N-1;i>=0;i--)for(let k=0;k<3;k++){let r=z[i*3+k];for(let j=F.upperOffsets[i];j<F.upperOffsets[i+1];j++)r-=F.upperValues[j]*z[F.upperIndices[j]*3+k];z[i*3+k]=r/F.diagonal[i];}
      let denominator=0,max=0;for(let i=0;i<N;i++){for(let k=0;k<3;k++)denominator+=grad[G.free[i]*3+k]*z[i*3+k];max=Math.max(max,Math.hypot(...z.subarray(i*3,i*3+3)));}
      if(!(denominator>1e-20)){state.volumeProjection.reason='singular-volume-gradient';return;}
      const multiplier=error/denominator,trust=Math.min(1,G.cellWorld/Math.max(Math.abs(multiplier)*max,1e-20));saved.set(G.x);let accepted=false;
      for(let ls=0;ls<12;ls++){
        const alpha=multiplier*trust*2**-ls;
        for(let i=0;i<N;i++)for(let k=0;k<3;k++)G.x[G.free[i]*3+k]=saved[G.free[i]*3+k]+alpha*z[i*3+k];
        const next=volumeData(state).volume;if(next>0&&Math.abs(Q.referenceVolume-next)<Math.abs(error)){accepted=true;break;}
      }
      if(!accepted){G.x.set(saved);state.volumeProjection.reason='volume-step-rejected';return;}
    }
  }
  root.DeformationGraph=Object.freeze({prepare,rotation,pose,point,deform,volumeData,surfaceAudit});
})(globalThis);
