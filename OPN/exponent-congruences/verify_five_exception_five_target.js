"use strict";
// Node.js, standard library only. Exact five non-Euler exceptions with 3 absent.
// Requires the mathematical reductions in five-exception-fixed-core.md.
function coarseCases() {
{
const X=1e6,P=500;
const is=new Uint8Array(X+1).fill(1,2),rem=new Float64Array(X+1),fs=new Map(),ps=[];
for(let p=2;p*p<=X;p++)if(is[p])for(let j=p*p;j<=X;j+=p)is[j]=0;
for(let p=5;p<=X;p++)if(is[p]){ps.push(p);rem[p]=p*p+p+1;fs.set(p,[]);}
function pow(a,e,p){let b=1;while(e){if(e%2)b=b*a%p;a=a*a%p;e=Math.floor(e/2);}return b;}
for(const p of [3,...ps]){let a=[];if(p===3)a=[1];else if(p%3===1){for(let k=2;k<p;k++){let z=pow(k,(p-1)/3,p);if(z!==1){a=[z,p-1-z];break;}}}for(const z of a)for(let j=z;j<=X;j+=p)if(rem[j]){fs.get(j).push(p);while(rem[j]%p===0)rem[j]/=p;}}
const caps=new Map();
for(const p of ps){let f=fs.get(p);if(rem[p]>1)f.push(rem[p]);let l=p%3===1?p:f.at(-1);if(p>P)(caps.get(l)??caps.set(l,[]).get(l)).push(p);}

const F=(a,b=1)=>[BigInt(a),BigInt(b)], mul=(a,b)=>[a[0]*b[0],a[1]*b[1]], add=(a,b)=>[a[0]*b[1]+a[1]*b[0],a[1]*b[1]], cmp=(a,b)=>a[0]*b[1]-a[1]*b[0], prod=as=>as.reduce(mul,F(1)), hp=p=>F(p,p-1);
function req(v,msg){if(!v)throw Error(msg);}
function froot(p){let l=1000,h=p*1000;while(h-l>1){let k=Math.floor((l+h)/2);if(BigInt(k)**6n<=BigInt(p)*10n**18n)l=k;else h=k;}return F(l,l-1000);}

const weightBound=F(1006,1000), W=new Map([...caps].map(([l,a])=>[l,prod(a.map(hp))]));
for(const weight of W.values())req(cmp(weight,weightBound)<0,"tail group below 1.006");
const small=ps.filter(p=>p<=P), other=small.filter(p=>p%3===2&&p!==5&&p!==11);
const cand=[...new Set([...small.filter(p=>p%3===1),...small.filter(p=>p%3===2).flatMap(p=>fs.get(p))])].sort((a,b)=>a-b);

const specs=[["13 present",[7,19,31,13],0,[]]];
const candSet=new Set(cand);
const stats=[];
const topWeights=[...W.entries()].sort((a,b)=>{let d=cmp(b[1],a[1]);return d>0n?1:d<0n?-1:0;}).slice(0,12);
for(const [name,fixedT,e,banned] of specs){
 const setup=cand.filter(p=>!fixedT.includes(p)&&!banned.includes(p));const free=6-e-fixedT.length;
 const anonymous=ps.filter(p=>p%3===1&&!candSet.has(p)&&!banned.includes(p)).slice(0,6);
 let n=0,surv=0,best=F(0),worst,remain=[];
 function rec(i,t){
  const T=[...fixedT,...t], S=new Set(T), L=new Set([7,19,31,...T]);
  const reg=other.filter(p=>fs.get(p).every(r=>S.has(r))), opt=other.filter(p=>!reg.includes(p));
  const active=[5,7,11,...reg,...T.filter(p=>p!==7&&p<=P)];
  for(let u=0;u<=free-t.length;u++){
   if(T.length+u<5||(!T.some(p=>p%12===1)&&!u))continue;
   const used=opt.slice(0,free-t.length-u);
   const ub=prod([...active.map(hp),...used.map(hp),...[...L].map(r=>W.get(r)??F(1)),...topWeights.filter(([r])=>!L.has(r)).slice(0,6-L.size).map(([,w])=>w)]);
   const d=add(F(6,X),mul(F(X+2,X*10000),prod([...T,...anonymous.slice(0,u)].map(froot))));
   const threshold=add(ub,mul(F(2),d));
   n++;if(cmp(threshold,F(2))>=0){surv++;remain.push([T,u]);}
   if(cmp(threshold,best)>0){best=threshold;worst=[T,u,active,used,Number(ub[0])/Number(ub[1]),Number(d[0])/Number(d[1])];}
  }
  if(t.length<free)for(let j=i;j<setup.length;j++)rec(j+1,[...t,setup[j]]);
 }
 rec(0,[]);
 stats.push({name,n,surv,threshold:Number(best[0])/Number(best[1]),worst,remain});
} 



const cases=stats[0].remain.filter(([T,u])=>T.length===5&&!u).map(([T])=>T);
if(cases.length!==73)throw Error("five-target count");
return cases;
}
}
function screen(cases) {
// Run with Node.js, taking the exact candidate-screen JSON from stdin.
const X=10n**16n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
const ok=Array.from({length:576},(_,i)=>i).reduce((s,i)=>(s.add(i*i%576),s),new Set());
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n;if(!ok.has(Number(d%576n)))return;let r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
function ints(p,role,optional){
 p=BigInt(p);const v=optional?[[F(1),F(1)]]:[];
 for(let e=p===7n?22:role==="Euler"?1:2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n)]);break;}
  v.push([val,val]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");
const hiThreshold=F(2);

function req(condition,message){if(!condition)throw Error(message);}


function order(p,R){
 const r=BigInt(R), x=p%r; if(!x)throw Error("nonunit");
 const pow=(e)=>{let a=x,b=1n;for(let n=BigInt(e);n;n>>=1n,a=a*a%r)if(n&1n)b=b*a%r;return b};
 let ord=R-1,k=ord;
 for(let d=2;d*d<=k;d++)if(k%d===0){while(k%d===0)k/=d;while(ord%d===0&&pow(ord/d)===1n)ord/=d;}
 if(k>1)while(ord%k===0&&pow(ord/k)===1n)ord/=k;
 return ord;
}

function regular(T,p,optional){
 const es=[...new Set(T.map(r=>order(p,r)).filter(n=>n%6===3))].sort((a,b)=>a-b).map(n=>n-1);
 const choices=optional?[[F(1),F(1)]]:[];
 let cap=0,maxpow=1n;
 while((p-1n)*maxpow<=10n**10n){cap++;maxpow*=p;}
 for(const e of es){
  let power=p**BigInt(Math.min(e,cap)), h=F(p*power-1n,(p-1n)*power);
  if(e>=cap){choices.push([h,F(p,p-1n)]);break;}
  choices.push([h,h]);
 }
 return choices;
}
function prime(p){for(let j=2;j*j<=p;j++)if(p%j===0)return false;return true;}
function boxes(T,pi,R,q,B){
 const count=T.length;
 const fine=count===5?140n:360n;
 if((X+2n)*8n**BigInt(count)>=fine*X*3n**BigInt(count))throw Error("fine tail bound");
 const loThreshold=F(2n*(10n**8n-fine),10n**8n);
 const domains=[
  ...T.map(p=>ints(p,p===pi?"Euler":"ex",false)),
  ...(q?[ints(q,"ex",false)]:[]),
  ...R.filter(p=>BigInt(q)!==p).map(p=>regular(T,p,p!==5n&&p!==11n))
 ];
 if(domains.some(x=>!x.length))return null;
 const lower=Array(domains.length+1);lower[domains.length]=F(1);
 for(let i=domains.length-1;i>=0;i--)lower[i]=mul(lower[i+1],domains[i][0][0]);
 const suffix=Array(domains.length+1);suffix[domains.length]=q||T.length===6?F(1):F(B+1,B);
 for(let i=domains.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],domains[i].at(-1)[1]);
 let best=F(0),num=0;
 function go(i,lo,hi){if(num)return;
  if(cmp(mul(lo,lower[i]),hiThreshold)>0||cmp(mul(hi,suffix[i]),loThreshold)<0)return;
  if(i===domains.length){num++;if(cmp(hi,best)>0)best=hi;return;}
  for(const [a,b] of domains[i])go(i+1,mul(lo,a),mul(hi,b));
 }
 go(0,F(1),F(1));if(!num)return null;
 const gap=loThreshold[0]*best[1]-loThreshold[1]*best[0];
 return [q,num,gap>0n?String((loThreshold[0]*best[1])/gap):"unbounded"];
}



const B=500;
const qs=Array.from({length:B},(_,i)=>i+1).filter(p=>p>11&&p%3===2&&prime(p));
let n=0,out=[];
for(const T of cases){
 const R=smooth(T).filter(p=>T.every(t=>p%BigInt(t)!==0n));
 if(!R.includes(5n)||!R.includes(11n))throw Error("mandatory preimage");
 for(const pi of T.filter(p=>p%12===1)){
  n++;
  const remaining=[...qs,0].filter(q=>boxes(T,pi,R,q,B));
  if(remaining.length)out.push([T,pi,R.map(String),remaining]);
 }
}
if(n!==107||out.length!==1||out[0][0].join(",")!=="7,19,31,13,79"||out[0][1]!==13||out[0][3].join(",")!=="311,0")throw Error("unexpected explicit target screen");
return {n,out};
}
function lastBranch() {
// Run with Node.js, taking the exact candidate-screen JSON from stdin.
const X=10n**16n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
const ok=Array.from({length:576},(_,i)=>i).reduce((s,i)=>(s.add(i*i%576),s),new Set());
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n;if(!ok.has(Number(d%576n)))return;let r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
let tightened=false;
function ints(p,role,optional){
 p=BigInt(p);const v=optional?[[F(1),F(1),0,false]]:[];
 for(let e=p===7n?22:role==="Euler"?(tightened?9:1):2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n),e,true]);break;}
  v.push([val,val,e,false]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");
const hiThreshold=F(2);



function order(p,R){
 const r=BigInt(R), x=p%r; if(!x)throw Error("nonunit");
 const pow=(e)=>{let a=x,b=1n;for(let n=BigInt(e);n;n>>=1n,a=a*a%r)if(n&1n)b=b*a%r;return b};
 let ord=R-1,k=ord;
 for(let d=2;d*d<=k;d++)if(k%d===0){while(k%d===0)k/=d;while(ord%d===0&&pow(ord/d)===1n)ord/=d;}
 if(k>1)while(ord%k===0&&pow(ord/k)===1n)ord/=k;
 return ord;
}

function regular(T,p,optional){
 if(tightened){
  if(p===653n)return [[F(1),F(1),0,false]];
  if([5n,11n,23n,191n].includes(p)){const t=F(p*p+p+1n,p*p);return [[t,t,2,false]];}
  throw Error("unrecognized regular preimage");
 }
 const es=[...new Set(T.map(r=>order(p,r)).filter(n=>n%6===3))].sort((a,b)=>a-b).map(n=>n-1);
 const choices=optional?[[F(1),F(1),0,false]]:[];
 let cap=0,maxpow=1n;
 while((p-1n)*maxpow<=10n**10n){cap++;maxpow*=p;}
 for(const e of es){
  let power=p**BigInt(Math.min(e,cap)), h=F(p*power-1n,(p-1n)*power);
  if(e>=cap){choices.push([h,F(p,p-1n),e,true]);break;}
  choices.push([h,h,e,false]);
 }
 return choices;
}
function prime(p){for(let j=2;j*j<=p;j++)if(p%j===0)return false;return true;}
function boxes(T,pi,R,q,B){
 const count=T.length;
 const fine=count===5?140n:360n;
 if((X+2n)*8n**BigInt(count)>=fine*X*3n**BigInt(count))throw Error("fine tail bound");
 const loThreshold=F(2n*(10n**8n-fine),10n**8n);
 const bases=[...T,...(q?[q]:[]),...R.filter(p=>BigInt(q)!==p).map(Number)];
 const domains=[
  ...T.map(p=>ints(p,p===pi?"Euler":"ex",false)),
  ...(q?[ints(q,"ex",false)]:[]),
  ...R.filter(p=>BigInt(q)!==p).map(p=>regular(T,p,p!==5n&&p!==11n))
 ];
 if(domains.some(x=>!x.length))return null;
 const lower=Array(domains.length+1);lower[domains.length]=F(1);
 for(let i=domains.length-1;i>=0;i--)lower[i]=mul(lower[i+1],domains[i][0][0]);
 const suffix=Array(domains.length+1);suffix[domains.length]=q||T.length===6?F(1):F(B+1,B);
 for(let i=domains.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],domains[i].at(-1)[1]);
 let best=F(0),num=0,seen=new Map();
 function go(i,lo,hi,es){
  if(cmp(mul(lo,lower[i]),hiThreshold)>0||cmp(mul(hi,suffix[i]),loThreshold)<0)return;
  if(i===domains.length){if(es[bases.indexOf(13)]?.[0]===1&&[191,653].some(p=>es[bases.indexOf(p)]?.[0]))return;
 num++;if(cmp(hi,best)>0)best=hi;for(let j=0;j<bases.length;j++)(seen.get(bases[j])??seen.set(bases[j],new Set()).get(bases[j])).add(es[j][0]+(es[j][1]?'+':''));return;}
  for(const [a,b,e,t] of domains[i])go(i+1,mul(lo,a),mul(hi,b),[...es,[e,t]]);
 }
 go(0,F(1),F(1),[]);if(!num)return null;
 const gap=loThreshold[0]*best[1]-loThreshold[1]*best[0];
 return [q,num,best.map(String),Object.fromEntries([...seen].map(([p,e])=>[p,[...e]])),gap>0n?String((loThreshold[0]*best[1])/gap):"unbounded"];
}




function req(condition,message){if(!condition)throw Error(message);}
const T=[7,19,31,13,79], R=smooth(T).sort((a,b)=>a<b?-1:a>b?1:0);
req(R.join(",")==="5,11,23,191,653","regular preimages");
for(const p of [191n,653n])req((p*p+p+1n)%169n===0n,"13-adic bound");
const fixed=boxes(T,13,R,311,500), unknown=boxes(T,13,R,0,500);
req(fixed[1]===24&&fixed[3][5].join(",")==="8,38+","a5 bounds for 311");
const lengths=[...new Set(T.map(r=>order(5n,r)).filter(m=>m%6===3))].sort((a,b)=>a-b);
req(lengths.join(",")==="3,9,39","admissible lengths at 5");
req(5n**6n+5n**3n+1n===19n*829n&&prime(829)&&829%3===1,"Phi9(5)");
const numerator=(5n**39n-1n)*4n,den=(5n**13n-1n)*(5n**3n-1n),phi39=numerator/den;
const gcd=(a,b)=>{while(b)[a,b]=[b,a%b];return a;};
req(numerator%den===0n&&gcd(phi39,5n**3n-1n)===1n&&gcd(phi39,5n**13n-1n)===1n,"all Phi39 primes have order 39");
let residue=phi39;for(const p of T.map(BigInt))while(residue%p===0n)residue/=p;
req(residue>1n,"Phi39 needs another target");
req(unknown[1]===12&&unknown[4]==="14798","unknown exception bound");
for(const [p,values] of [[5,"2"],[11,"2"],[23,"2"],[191,"2"],[653,"0"],[13,"9+"]])req(unknown[3][p].join(",")===values,"unknown exception profile");
tightened=true;
const qList=Array.from({length:14798},(_,i)=>i+1).filter(p=>p>500&&p%3===2&&prime(p));
const candidates=qList.filter(q=>boxes(T,13,R,q,500));
req(qList.length===828&&candidates.length===25,"bounded prime intervals");
return candidates;
}
function closeLastBranch(candidates) {
// Run with Node.js, taking the exact candidate-screen JSON from stdin.
const X=10n**30n, max=X*X+X+1n;
const G=([a,b])=>{let x=a,y=b;while(y)[x,y]=[y,x%y];return [a/x,b/x];};
const F=(a,b=1n)=>[BigInt(a),BigInt(b)], mul=(a,b)=>G([a[0]*b[0],a[1]*b[1]]), cmp=(a,b)=>a[0]*b[1]-a[1]*b[0];
function isqrt(n){if(n<2n)return n;let g=1n<<BigInt(Math.ceil(n.toString(2).length/2));let h=(g+n/g)>>1n;while(h<g){g=h;h=(g+n/g)>>1n;}return g;}
const ok=Array.from({length:576},(_,i)=>i).reduce((s,i)=>(s.add(i*i%576),s),new Set());
function smooth(T){
 let ans=[];
 function go(i,v){
  if(i===T.length){let d=4n*v-3n;if(!ok.has(Number(d%576n)))return;let r=isqrt(d),q=(r-1n)/2n;if(r*r===d&&q%6n===5n&&q<=X)ans.push(q);return;}
  for(let w=v;w<=max;w*=BigInt(T[i]))go(i+1,w);
 }
 go(0,1n);return ans;
}
const tightened=true;
function ints(p,role,optional){
 p=BigInt(p);const v=optional?[[F(1),F(1),0,false]]:[];
 for(let e=p===7n?22:role==="Euler"?(tightened?9:1):2;;e+=role==="Euler"?4:2){
  if(((e+1)%3===0)!==(role==="reg"))continue;
  let pow=p**BigInt(e),val=F(p*pow-1n,(p-1n)*pow);
  if((p-1n)*pow>10n**10n){v.push([val,F(p,p-1n),e,true]);break;}
  v.push([val,val,e,false]);
 }
 return v;
}
if(8n**4n>=7n*5n**4n||(X+2n)*8n**5n>=140n*X*3n**5n)throw Error("fine tail bound");
const hiThreshold=F(2);



function order(p,R){
 const r=BigInt(R), x=p%r; if(!x)throw Error("nonunit");
 const pow=(e)=>{let a=x,b=1n;for(let n=BigInt(e);n;n>>=1n,a=a*a%r)if(n&1n)b=b*a%r;return b};
 let ord=R-1,k=ord;
 for(let d=2;d*d<=k;d++)if(k%d===0){while(k%d===0)k/=d;while(ord%d===0&&pow(ord/d)===1n)ord/=d;}
 if(k>1)while(ord%k===0&&pow(ord/k)===1n)ord/=k;
 return ord;
}

function regular(T,p,optional){
 if(tightened){
  if(p===653n)return [[F(1),F(1),0,false]];
  if([5n,11n,23n,191n].includes(p)){const t=F(p*p+p+1n,p*p);return [[t,t,2,false]];}
  throw Error("unrecognized regular preimage");
 }
 const es=[...new Set(T.map(r=>order(p,r)).filter(n=>n%6===3))].sort((a,b)=>a-b).map(n=>n-1);
 const choices=optional?[[F(1),F(1),0,false]]:[];
 let cap=0,maxpow=1n;
 while((p-1n)*maxpow<=10n**10n){cap++;maxpow*=p;}
 for(const e of es){
  let power=p**BigInt(Math.min(e,cap)), h=F(p*power-1n,(p-1n)*power);
  if(e>=cap){choices.push([h,F(p,p-1n),e,true]);break;}
  choices.push([h,h,e,false]);
 }
 return choices;
}
function prime(p){for(let j=2;j*j<=p;j++)if(p%j===0)return false;return true;}
function boxes(T,pi,R,q,B){
 const count=T.length;
 const fine=count===5?140n:360n;
 if((X+2n)*8n**BigInt(count)>=fine*X*3n**BigInt(count))throw Error("fine tail bound");
 const loThreshold=F(2n*(10n**15n-fine),10n**15n);
 const bases=[...T,...(q?[q]:[]),...R.filter(p=>BigInt(q)!==p).map(Number)];
 const domains=[
  ...T.map(p=>ints(p,p===pi?"Euler":"ex",false)),
  ...(q?[ints(q,"ex",false)]:[]),
  ...R.filter(p=>BigInt(q)!==p).map(p=>regular(T,p,p!==5n&&p!==11n))
 ];
 if(domains.some(x=>!x.length))return null;
 const lower=Array(domains.length+1);lower[domains.length]=F(1);
 for(let i=domains.length-1;i>=0;i--)lower[i]=mul(lower[i+1],domains[i][0][0]);
 const suffix=Array(domains.length+1);suffix[domains.length]=q||T.length===6?F(1):F(B+1,B);
 for(let i=domains.length-1;i>=0;i--)suffix[i]=mul(suffix[i+1],domains[i].at(-1)[1]);
 let best=F(0),num=0,seen=new Map();
 function go(i,lo,hi,es){
  if(cmp(mul(lo,lower[i]),hiThreshold)>0||cmp(mul(hi,suffix[i]),loThreshold)<0)return;
  if(i===domains.length){if(es[bases.indexOf(13)]?.[0]===1&&[191,653].some(p=>es[bases.indexOf(p)]?.[0]))return;
 num++;if(cmp(hi,best)>0)best=hi;for(let j=0;j<bases.length;j++)(seen.get(bases[j])??seen.set(bases[j],new Set()).get(bases[j])).add(es[j][0]+(es[j][1]?'+':''));return;}
  for(const [a,b,e,t] of domains[i])go(i+1,mul(lo,a),mul(hi,b),[...es,[e,t]]);
 }
 go(0,F(1),F(1),[]);if(!num)return null;
 const gap=loThreshold[0]*best[1]-loThreshold[1]*best[0];
 return [q,num,best.map(String),Object.fromEntries([...seen].map(([p,e])=>[p,[...e]])),gap>0n?String((loThreshold[0]*best[1])/gap):"unbounded"];
}




function req(condition,message){if(!condition)throw Error(message);}
const T=[7,19,31,13,79], R=smooth(T).sort((a,b)=>a<b?-1:a>b?1:0);
req(R.join(",")==="5,11,23,191,653","regular preimages");
for(const p of [191n,653n])req((p*p+p+1n)%169n===0n,"13-adic bound");

const remains=candidates.filter(q=>boxes(T,13,R,q,500));
req(remains.length===0,"all five-target cases excluded at 10^30");
return R.map(String);
}
const cases=coarseCases();
console.log("coarse five-target sets:", cases.length);
const screened=screen(cases);
console.log("explicit interval screen:", JSON.stringify(screened));
const candidates=lastBranch();
console.log("small-tail survivors:", candidates.length);
console.log("relaxed regular preimages through 10^30:", closeLastBranch(candidates).join(", "));
console.log("five-target cases remaining: 0");
