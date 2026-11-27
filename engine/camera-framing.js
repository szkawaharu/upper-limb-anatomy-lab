/* Pure orthographic framing shared by initial fit and clipping guard.
 * All dimensions are CSS pixels; input points must already be posed in world
 * coordinates, including mirrored arms. A fit is a display operation only. */
(function(root){
"use strict";
const M=root.AnatomyMath,{dot,cross,norm,add,sub,mul}=M;
function basis(yaw,pitch){
 const outward=[Math.sin(yaw)*Math.cos(pitch),Math.sin(pitch),Math.cos(yaw)*Math.cos(pitch)];
 const right=norm(cross([0,1,0],outward)),up=cross(outward,right);
 return {right,up};
}
function layout(w,h,pads){
 w=Math.max(1,w);h=Math.max(1,h);
 // Very small or transitional viewport: preserve a positive usable rectangle.
 let [l,r,t,b]=pads.map(x=>Math.max(0,x));
 const kx=Math.min(1,Math.max(0,w-1)/Math.max(1,l+r)),ky=Math.min(1,Math.max(0,h-1)/Math.max(1,t+b));
 l*=kx;r*=kx;t*=ky;b*=ky;
 return {w,h,l,r,t,b,areaW:w-l-r,areaH:h-t-b,cx:(l+w-r)/2,cy:(t+h-b)/2};
}
function screenBounds(points,camera,w,h){
 const {right,up}=basis(camera.yaw,camera.pitch);
 let minX=Infinity,maxX=-Infinity,minY=Infinity,maxY=-Infinity;
 for(const p of points){
  const d=sub(p,camera.target),x=w/2+dot(d,right)*h/camera.zoom,y=h/2-dot(d,up)*h/camera.zoom;
  minX=Math.min(minX,x);maxX=Math.max(maxX,x);minY=Math.min(minY,y);maxY=Math.max(maxY,y);
 }
 return {minX,maxX,minY,maxY};
}
function clipped(points,camera,w,h,pads){
 if(!points.length)return false;
 const a=layout(w,h,pads),s=screenBounds(points,camera,a.w,a.h);
 // Small tolerance for projection roundoff, not a relaxation of safe bounds.
 return s.minX<a.l-1e-6||s.maxX>a.w-a.r+1e-6||s.minY<a.t-1e-6||s.maxY>a.h-a.b+1e-6;
}
function fit(points,camera,w,h,pads,{minSpan=0,zoomScale=1}={}){
 if(!points.length)return null;
 const a=layout(w,h,pads),{right,up}=basis(camera.yaw,camera.pitch);
 let minX=Infinity,maxX=-Infinity,minY=Infinity,maxY=-Infinity;
 for(const p of points){const x=dot(p,right),y=dot(p,up);minX=Math.min(minX,x);maxX=Math.max(maxX,x);minY=Math.min(minY,y);maxY=Math.max(maxY,y)}
 const natural=Math.max(1.85,(maxX-minX)*a.h/a.areaW,(maxY-minY)*a.h/a.areaH)*1.04;
 // Use final zoom for centering. Never center for one span and apply another.
 const span=Math.max(minSpan,natural/zoomScale),base=mul(add(points[0],points.at(-1)),.5);
 const target=add(base,add(mul(right,(minX+maxX)/2-dot(base,right)+(a.w/2-a.cx)*span/a.h),
                          mul(up,(minY+maxY)/2-dot(base,up)+(a.cy-a.h/2)*span/a.h)));
 if(!Number.isFinite(span)||!target.every(Number.isFinite))return null;
 return {target,zoom:span,yaw:camera.yaw,pitch:camera.pitch};
}
function protect(points,camera,w,h,pads){
 // No recenter or zoom when all points already fit; never zoom in on gesture end.
 return clipped(points,camera,w,h,pads)?fit(points,camera,w,h,pads,{minSpan:camera.zoom}):null;
}
root.CameraFraming=Object.freeze({fit,protect,clipped,screenBounds,layout});
})(globalThis);
