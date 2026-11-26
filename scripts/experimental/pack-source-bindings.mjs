import{readFileSync,writeFileSync,mkdirSync}from'node:fs';import{gzipSync}from'node:zlib';import{createHash}from'node:crypto';
const bindings=JSON.parse(readFileSync('test-results/v076-source-clean/bindings.json')),out='test-results/v076-source-clean/packets';mkdirSync(out,{recursive:true});const manifest=[];
for(const B of bindings){const packet=gzipSync(Buffer.from(JSON.stringify(B)));writeFileSync(out+'/'+B.spec.id+'.json.gz',packet);manifest.push({id:B.spec.id,bytes:packet.length,sha256:createHash('sha256').update(packet).digest('hex'),meshPreparation:B.meshPreparation});}
writeFileSync(out+'/manifest.json',JSON.stringify(manifest,null,2));console.log('Compressed source bindings bytes',manifest.reduce((s,x)=>s+x.bytes,0));
