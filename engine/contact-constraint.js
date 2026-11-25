(function(root){'use strict';const {sub,add,mul,len,clamp,smooth}=root.AnatomyMath;
// Rest-relative ellipsoid envelope. This is a conservative radial proxy, not
// exact signed distance, collision response, or an inter-muscle solver.
function project(rest,p,t,env){if(!env||t<=0||t>=1)return p;const r=sub(rest,env.restCenter||env.center),v=sub(p,env.center),radii=env.radii,l0=len(r.map((x,i)=>x/radii[i])),l=len(v.map((x,i)=>x/radii[i]));if(l<1e-7)return p;const scale=env.scale||.00555,avg=(radii[0]+radii[1]+radii[2])/3,gate=smooth(t/.12)*smooth((1-t)/.12)*(env.proximal?1-smooth((t-.45)/.4):1),target=clamp(l,l0-(env.inMm??2)*scale/avg,l0+(env.outMm??6)*scale/avg),d=mul(v,target/l-1),cap=(env.maxMm??20)*scale,n=len(d);return add(p,mul(d,gate*Math.min(1,cap/Math.max(n,1e-12))))}
root.ContactConstraint=Object.freeze({project,kind:'rest-relative-ellipsoid-proxy'});
})(globalThis);
