// Offline verification against git blob IDs read from the pinned upstream tree.
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
const source=JSON.parse(readFileSync(new URL('../docs/validation/v076/atlas-sources.json',import.meta.url))),root=process.argv[2]||'/tmp/upper-limb-atlas';
for(const f of source.files){const b=readFileSync(join(root,f.file)),sha=createHash('sha1').update(Buffer.from('blob '+b.length+'\0')).update(b).digest('hex');if(sha!==f.gitBlobSha||b.length!==f.bytes)throw Error('Pinned atlas mismatch: '+f.file)}
console.log(JSON.stringify({atlasCommit:source.atlasCommit,verified:source.files.length,landmarkAccuracyVerified:false}));
