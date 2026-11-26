// Rejected production candidate: independent closest projection worsens atlas edge/area outliers.
// Diagnostic helper only. See docs/VALIDATION_V076.md.
export function closestMaterial(c, p, preferred = .5) {
    const a=c[0].map((x,k)=>c[3][k]-3*c[2][k]+3*c[1][k]-x),
      b=c[0].map((x,k)=>3*(c[2][k]-2*c[1][k]+x)),
      d=c[0].map((x,k)=>3*(c[1][k]-x));
    const loss=t=>{let s=0;for(let k=0;k<3;k++){const x=((a[k]*t+b[k])*t+d[k])*t+c[0][k]-p[k];s+=x*x}return s};
    let best=0,value=loss(0);
    const take=t=>{const v=loss(t);if(v<value||(v===value&&Math.abs(t-preferred)<Math.abs(best-preferred))){best=t;value=v}};
    const cells=32,values=Array.from({length:cells+1},(_,i)=>loss(i/cells));
    take(1);
    for(let i=0;i<=cells;i++)if(values[i]<=(values[i-1]??Infinity)&&values[i]<=(values[i+1]??Infinity)){
      let lo=Math.max(0,i-1)/cells,hi=Math.min(cells,i+1)/cells;
      for(let n=0;n<28;n++){const x=lo+(hi-lo)*.3819660112501051,y=lo+(hi-lo)*.6180339887498949;if(loss(x)<loss(y))hi=y;else lo=x}
      take((lo+hi)/2);
    }
    return best;
  }
