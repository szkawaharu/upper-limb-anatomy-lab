/* WebGL1 implementation of the same material-frame map as the CPU reference.
 * Up to nine slips; fixed array indexing for mobile GLSL ES 1.00 compilers. */
(function (root) {
  "use strict";
  function vertexShader() {
    function sample(name, array) {
      let code = `void ${name}(float l,float t,out vec3 c,out vec3 T){vec3 a=vec3(0.),b=a,d=a,e=a;`;
      for (let i = 0; i < 9; i++)
        code += `float w${i}=max(0.,1.-abs(l-${i}.));a+=${array}[${i * 4}]*w${i};b+=${array}[${i * 4 + 1}]*w${i};d+=${array}[${i * 4 + 2}]*w${i};e+=${array}[${i * 4 + 3}]*w${i};`;
      return (
        code +
        "float u=1.-t;c=a*u*u*u+b*3.*u*u*t+d*3.*u*t*t+e*t*t*t;T=3.*(b-a)*u*u+6.*(d-b)*u*t+3.*(e-d)*t*t;}"
      );
    }
    function frame(name, sample) {
      return `void ${name}(float l,float t,vec3 ref,out vec3 c,out vec3 T,out vec3 U,out vec3 N,out float tl,out float ul){vec3 v,x,y,z;${sample}(l,t,c,v);tl=length(v);T=v/max(tl,0.0000001);float l0=clamp(l-.25,0.,lanes-1.),l1=clamp(l+.25,0.,lanes-1.);${sample}(l0,t,x,z);${sample}(l1,t,y,z);v=(y-x)/max(l1-l0,0.000001);v-=T*dot(v,T);ul=length(v);if(ul<0.0000001)v=ref-T*dot(ref,T);U=v/max(length(v),0.0000001);N=normalize(cross(T,U));}`;
    }
    return `precision highp float;
attribute vec3 pos;attribute vec3 normal;attribute float weight;attribute float lane;attribute float materialT;attribute float originBlend;
uniform mat4 mvp,model,anchorA,anchorB,originAnchor;uniform vec4 qa,qb,originQA;uniform float originMixOn;uniform vec2 blendRange;uniform vec3 restControls[36],nowControls[36];uniform float lanes,sheet,tendonStart;uniform vec3 reference;
uniform float contactOn,contactProximal,contactUnilateral,unitScale,contactIn,contactOut,contactMax;uniform vec3 contactCenter,contactRestCenter,contactRadii;
varying vec3 N;varying vec3 W;
vec3 qrot(vec4 q,vec3 v){vec3 t=2.*cross(q.xyz,v);return v+q.w*t+cross(q.xyz,t);}
vec3 align(vec3 u,vec3 v,vec3 p,vec3 axis){float d=clamp(dot(u,v),-1.,1.);if(d<-.9999)return 2.*axis*dot(axis,p)-p;vec3 k=cross(u,v),a=cross(k,p);return p+a+cross(k,a)/max(1.+d,.000001);}
${sample("restCurve", "restControls")}${sample("nowCurve", "nowControls")}
${frame("restFrame", "restCurve")}${frame("nowFrame", "nowCurve")}
void main(){float t=weight;float k=originMixOn>.5?originBlend:0.;vec4 aq=normalize(mix(qa,originQA*(dot(qa,originQA)<0.?-1.:1.),k));vec4 ap=mix(anchorA*vec4(pos,1.),originAnchor*vec4(pos,1.),k);vec4 p=mix(ap,anchorB*vec4(pos,1.),t);vec3 n=mix(qrot(aq,normal),mat3(anchorB)*normal,t);
if(t>0.&&t<1.){vec3 rc,rt,ru,rn,sc,st,su,sn;float rtl,rul,stl,sul;restFrame(lane,materialT,reference,rc,rt,ru,rn,rtl,rul);nowFrame(lane,materialT,qrot(aq,reference),sc,st,su,sn,stl,sul);
if(rtl>.0000001&&stl>.0000001&&length(ru)>.5&&length(su)>.5){vec4 q=normalize(mix(aq,qb*(dot(aq,qb)<0.?-1.:1.),smoothstep(0.,1.,mix(blendRange.x,blendRange.y,materialT))));vec3 refT=qrot(q,rt),refU=qrot(q,ru);su=normalize(align(refT,st,refU,refU));sn=normalize(cross(st,su));float lt=clamp(stl/rtl,.35,2.8),lu=sheet>.5?clamp(sul/max(rul,.0000001),.4,2.5):inversesqrt(lt),ln=clamp(sheet>.5?1./(lt*lu):lu,.4,1.6),tw=smoothstep(tendonStart,1.,t);lu=mix(lu,1.,tw);ln=mix(ln,1.,tw);vec3 o=pos-rc,candidate=sc+st*(dot(o,rt)*lt)+su*(dot(o,ru)*lu)+sn*(dot(o,rn)*ln);float gate=smoothstep(0.,.12,t)*smoothstep(0.,.12,1.-t);p.xyz=mix(p.xyz,candidate,gate);n=mix(n,st*(dot(normal,rt)/lt)+su*(dot(normal,ru)/lu)+sn*(dot(normal,rn)/ln),gate);}
if(contactOn>.5){vec3 r=(pos-contactRestCenter)/contactRadii,v=p.xyz-contactCenter;float l0=length(r),l=length(v/contactRadii),avg=(contactRadii.x+contactRadii.y+contactRadii.z)/3.,gate=smoothstep(0.,.12,t)*smoothstep(0.,.12,1.-t)*(contactProximal>.5?1.-smoothstep(.45,.85,t):1.);if(l>.0000001){float target=contactUnilateral>.5?max(l,min(l0,1.)-contactIn*unitScale/avg):clamp(l,l0-contactIn*unitScale/avg,l0+contactOut*unitScale/avg);vec3 d=v*(target/l-1.);p.xyz+=d*gate*min(1.,contactMax*unitScale/max(length(d),.000000000001));}}
}vec4 wp=model*p;W=wp.xyz;N=normalize(mat3(model)*n);gl_Position=mvp*p;}`;
  }
  root.MuscleGPU = Object.freeze({ vertexShader, maxLanes: 9 });
})(globalThis);
