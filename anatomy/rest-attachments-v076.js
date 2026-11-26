/* Frozen rest-binding choices captured at 57e3530. These are the existing
 * estimated fields, NOT newly measured attachment footprints or pose-specific
 * corrections. A future anatomical revision must version these data explicitly.
 * No skeletal solver, current pose, camera, or deformation score is consulted.
 */
(function(root){
  'use strict';
  const choices=Object.freeze({
    v58RhomboidMajor:'original',v58RhomboidMinor:'midpoint',v58Levator:'original',
    v58Serratus:'smooth-midpoint',v58Subscap:'refined',v59TrapUpper:'original',
    v59TrapMiddle:'original',v59TrapLower:'smooth-midpoint',v59Lat:'refined',
    v60PecClav:'original',v60PecStern:'smooth-midpoint',v60PecAbd:'refined',
    v60PecMinor:'original',v60Subclavius:'smooth-midpoint',v60DeltAnt:'midpoint',
    v60DeltMid:'smooth-midpoint',v60DeltPost:'smooth-midpoint'
  });
  function field(id,original,refined,smooth){
    const choice=choices[id];
    if(!choice)throw Error('Missing versioned rest attachment field '+id);
    if(original.length!==refined.length)throw Error('Attachment field length');
    const fix=values=>Float32Array.from(values,(w,i)=>refined[i]<=1e-8||refined[i]>=1-1e-8?refined[i]:Math.max(0,Math.min(1,w)));
    if(choice==='refined')return Float32Array.from(refined);
    if(choice==='original')return fix(original);
    // Preserve Float32 rounding in the captured construction, including its
    // exact Dirichlet masks. This is a baseline-data separation, not retuning.
    const mid=fix(Float32Array.from(refined,(w,i)=>(w+original[i])*.5));
    return choice==='midpoint'?mid:smooth(mid,6);
  }
  root.RestAttachmentsV076=Object.freeze({choices,field,sourceCommit:'57e3530f077781d1502f15d478ef4e708495597a',confidence:'legacy-estimated-rest-fields'});
})(globalThis);
