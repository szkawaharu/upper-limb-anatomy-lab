"""Extract geometry-only reference paths. No force/activation parameters are
adopted. Coordinates remain in the cited model's body-local metre frames;
registration to BodyParts3D is a separate, uncompleted operation.
"""
import sys,json,hashlib
from pathlib import Path
import xml.etree.ElementTree as ET
src=Path(sys.argv[1]);out=Path(sys.argv[2]);root=ET.parse(src).getroot();vector=lambda s:[float(x)for x in s.split()]if s else None
muscles=[]
for m in root.findall('.//Millard2012EquilibriumMuscle'):
 points=[]
 for p in m.findall('./GeometryPath/PathPointSet/objects/*'):
  points.append({'type':p.tag,'name':p.attrib.get('name'),'bodyFrame':p.findtext('socket_parent_frame'),'locationMetres':vector(p.findtext('location'))})
 wraps=[{c.tag:(c.text or'').strip()for c in w}for w in m.findall('./GeometryPath/PathWrapSet/objects/*')]
 muscles.append({'name':m.attrib['name'],'points':points,'wraps':wraps})
markers=[{'name':m.attrib['name'],'bodyFrame':m.findtext('socket_parent_frame'),'locationMetres':vector(m.findtext('location'))}for m in root.findall('.//MarkerSet/objects/Marker')]
obstacles=[]
for b in root.findall('.//BodySet/objects/Body'):
 for w in b.findall('./WrapObjectSet/objects/*'):
  obstacles.append({'body':b.attrib['name'],'type':w.tag,'name':w.attrib['name'],'parameters':{c.tag:(c.text or'').strip()for c in w}})
report={'adoptedIntoApplication':False,'registeredToBodyParts3D':False,'units':'metres; source body-local frames','sourceRepository':'https://github.com/ComputationalBiomechanicsLab/rmr-solver','sourceCommit':'b47b06cd23c872e69f92687270e47d4c560bb14b','sourceTag':'v1.1','sourcePath':'OpenSim Models/for RMR solver/TSM_subject_noWeight.osim','modelSha256':hashlib.sha256(src.read_bytes()).hexdigest(),'license':'CC BY 4.0','copyright':'2023 Italo Belli, Sagar Joshi, J. Micah Prendergast, Irene Beck, Cosimo dell Santina, Luka Peternel, Ajay Seth','attribution':'Modified thoracoscapular shoulder model; Belli et al. 2023 DOI10.1371/journal.pone.0295003; based on Seth et al. 2019 DOI10.3389/fnbot.2019.00090','forceModelAdopted':False,'muscles':muscles,'markers':markers,'wrappingObjects':obstacles}
out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'muscles':len(muscles),'markers':len(markers),'wrapObjects':len(obstacles),'registered':False}))
