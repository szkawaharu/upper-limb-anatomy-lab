/* BodyParts3D 3.0 pinned surfaces, mm. These are calibrated estimates, not
 * measured subject-specific functional joint centres. See calibration audit. */
(function(root){
'use strict';
root.ShoulderRigData=Object.freeze({
 version:'0.76.0',atlasCommit:'2c83c698788cff57f6317bddcf712e7acfbea7d8',license:'CC BY-SA 2.1 Japan',
 landmarks:{
  IJ:[.28306,-141.64795,1342.02948],
  C7:[-.85246,-24.25698,1384.21708],
  T8:[-.92073,-2.09918,1182.50001],
  PX:[-.64925,-197.22623,1170.09708],
  SC:[-14.91024,-129.66467,1341.37453], // medial clavicle terminal-cap centroid
  AC:[-147.16389,-62.15942,1357.28268], // paired clavicle/acromion surface patch
  GH:[-162.59912,-72.14044,1325.58960], // robust proximal-medial humeral sphere
  RH:[-232,-68,1039], DU:[-240,-110,825],
  EL_LAT:[-230.01776,-62.78176,1058.07301], EL_MED:[-213.38389,-65.33045,1056.05264],
  EL:[-221.70083,-64.05610,1057.06283], // midpoint of fitted distal condylar centres
  AA:[-155.36700,-51.18210,1356.91003],
  TS:[-63.87020,-6.27689,1331.10999],
  AI:[-91.16020,.02043,1201.5],
  PC:[-136.90500,-95.10980,1332.53003],
 },
 elbowAxis:[[-230.01776,-62.78176,1058.07301],[-213.38389,-65.33045,1056.05264]],
 headSphere:{radiusMm:21.51355,fitRmsMm:1.35566,method:'proximal-medial source cap, robust sphere'},
 confidence:{thorax:'IJ/PX source surface extrema; C7/T8 posterior process caps, approximate ISB landmarks',SC:'surface-cap estimate; counterpart distance 2.556mm',AC:'paired-surface estimate; patch vertex gaps 0.216–1.715mm',GH:'sphere fit estimate; radius21.514mm, residual RMS1.356mm',AA:'manual seed snapped to source surface',TS:'manual seed snapped to source surface',AI:'source inferior extremum',PC:'manual seed snapped to source surface',EL:'distal condylar cylinder-section estimate; fit RMS0.603/0.633mm; not measured functional axis'},
 // Approximate mean points read from Ludewig et al. 2009, Fig5A/B, pp383–384.
 // Columns: HT elevation, SC retraction/elevation/posterior roll,
 // AC internal rotation/upward rotation/posterior tilt. Degrees, positive names.
 // 0° is TableI. End-ROM points are group endpoint means, not matched samples.
 curves:{
  abduction:[[0,19.2,5.9,.1,60,2.5,-8.4],[15,27,13,1,55.6,3.1,-4],[30,31,14.5,3,58.4,6.3,-.5],[60,37,16.7,9,60.8,12,7],[90,40,18,15.8,61.6,16.3,13.1],[120,43.7,19.5,25,62.3,17.1,14.7],[136,45.1,19.4,31,63.5,16.9,15.5]],
  scapular:[[0,19.2,5.9,.1,60,2.5,-8.4],[15,23.4,10.6,.5,57.3,4.2,-5],[30,25.6,11.5,2,59.2,8.3,-2.1],[60,29.7,14.3,7,61.4,13.5,4.3],[90,32.2,15.6,15,63,15.2,6],[120,36.5,16.9,24,63.6,15.1,10.5],[139,38.9,17.7,32,64.3,15.8,13.8]],
  flexion:[[0,19.2,5.9,.1,60,2.5,-8.4],[15,18.8,10.5,-1.3,58.5,7.7,-3.5],[30,21,10,0,60.8,10,-1.2],[60,23,12.5,7,64,14.2,4.7],[90,24.7,13.9,15,65.5,15.3,7],[120,31.5,14.7,24.8,66.2,15.4,11.5],[131,34.7,14.6,30,68.5,15.5,16.1]],
 },
 curveSource:'Ludewig2009 DOI10.2106/JBJS.G.01483 Fig5A/B; hand digitization approx ±1–2deg, n12; SC data n11',
 beyondData:'C1 plateau of last group endpoint; high elevation is operable but unvalidated',
 priorStdDeg:[4,2,3,4,3,5], // numerical prior scales informed by figure uncertainty; NOT ROM limits
 clinicalValidation:false
});
})(globalThis);
